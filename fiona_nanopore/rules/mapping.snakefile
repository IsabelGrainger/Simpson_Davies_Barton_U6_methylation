def get_reference_fn(wc):
    if wc.ref_type == 'transcriptome':
        return config['transcriptome_fasta_fn']
    else:
        return config['genome_fasta_fn']


def get_junc_bed(wc):
    if wc.ref_type == 'genome':
        return 'juncs/merged.filtered.bed'
    else:
        return []


def minimap2_general_parameters(wc):
    params = ['-a', '-L', '--cs=long']
    if wc.ref_type == 'transcriptome':
        params += ['-k14', '--for-only', '--secondary=no']
    else:
        intron_size = config['minimap2_parameters'].get(
            'max_intron_size', 200_000
        )
        params += ['-k14', '-uf', '-w5', '--splice',
                   '-g2000', f'-G{intron_size}',
                   '--end-seed-pen=15',
                   '-A1', '-B2', '-O2,32', '-E1,0', '-C9',
                   '--splice-flank=yes', '-z200']
    return ' '.join(params)


def minimap2_junc_parameters(wc, input):
    if wc.ref_type == 'genome':
        junc_bonus = config['minimap2_parameters'].get(
            'junction_bonus', 10
        )
        params = f'--junc-bed {input.bed} --junc-bonus={junc_bonus}'
    else:
        params = ''
    return params

            
rule map_with_minimap:
    input:
        fastq='basecalled_data/{sample_name}.dna.fastq',
        bed=get_junc_bed,
        reference=get_reference_fn,
    output:
        'aligned_data/{sample_name}.{ref_type}.bam'
    threads: 
        8
    params:
        mm2_params=minimap2_general_parameters,
        junc_params=minimap2_junc_parameters,
    conda:
        'env_yamls/minimap2.yml'
    shell:
        '''
        ../scripts/minimap2/minimap2 -t{threads} \
          {params.mm2_params} \
          {params.junc_params} \
          {input.reference} \
          {input.fastq} |
        samtools view -bS - |
        samtools sort -m 1G -@ {threads} -o - - > {output}
        samtools index {output}
        '''


rule run_2passtools_score:
    input:
        bam='aligned_data/{sample_name}.genome_1pass.bam',
        reference=config['genome_fasta_fn']
    output:
        bed='juncs/{sample_name}.all.bed'
    threads:
        8
    conda:
        'env_yamls/2passtools.yml'
    shell:
        '''
        2passtools score \
          -j 4 -d 20 -m "GTAG|GCAG|ATAG" \
          -w 128 -k 6 -lt 0.1 -ht 0.6 \
          --stranded \
          -f {input.reference} \
          -o {output.bed} \
          {input.bam}
        '''


def merge_input(wc):
    sample_names = glob_wildcards(
        'raw_data/{sample_name}'
    ).sample_name
    return expand(
        'juncs/{sample_name}.all.bed',
        sample_name=sample_names
    )
        

rule run_2passtools_merge:
    input:
        bams=merge_input,
        reference=config['genome_fasta_fn']
    output:
        all_='juncs/merged.all.bed',
        filtered='juncs/merged.filtered.bed',
    threads:
        8
    conda:
        'env_yamls/2passtools.yml'
    shell:
        '''
        2passtools merge \
          -j 4 -d 20 -m "GTAG|GCAG|ATAG" \
          -w 128 -k 6 -lt 0.1 -ht 0.6 \
          -f {input.reference} \
          -o {output.all_} \
          {input.bams}
        2passtools filter -o {output.filtered} {output.all_}
        '''


def sample_name_subset(cond):
    sample_names = glob_wildcards(
        'raw_data/{sample_name}'
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.rsplit('_', 1) == cond]
    return cond_sample_names


def pool_input(wc):
    return expand(
        'aligned_data/{sample_name}.{ref_type}.bam',
        sample_name=sample_name_subset(wc.cond),
        ref_type=wc.ref_type
    )


rule pool_bams:
    input:
        pool_input
    output:
        'aligned_data/pooled/{cond}.{ref_type}.bam'
    threads: 12
    conda:
        'env_yamls/minimap2.yml'
    shell:
        '''
        samtools merge -r -@ {threads} {output} {input}
        samtools index {output}
        '''


def get_gene_track_input(wc):
    sample_names = glob_wildcards('raw_data/{sample_name}').sample_name
    conds = set([sn.rsplit('_', 1)[0] for sn in sample_names])
    bams = expand('aligned_data/pooled/{cond}.genome.bam', cond=conds)
    return {
        'illumina_psi_fit': illumina(f'splicing/denovo/psi_fit.csv'),
        'bams': bams,
        'fasta': ancient(config['genome_fasta_fn']),
    }


rule generate_splicing_gene_tracks:
    input:
        unpack(get_gene_track_input)
    output:
        gene_tracks=directory('figures/splicing/gene_tracks/{gene_id}_splicing_gene_tracks')
    conda:
        'env_yamls/nb_seqlogos.yaml'
    log:
        notebook='notebook_processed/{gene_id}_splicing_gene_tracks.py.ipynb'
    notebook:
        'notebook_templates/splicing_gene_tracks.py.ipynb'