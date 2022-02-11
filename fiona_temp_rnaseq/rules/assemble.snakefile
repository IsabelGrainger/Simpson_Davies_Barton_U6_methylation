def match_nanopore_input(wc):
    # nanopore data is single temp so only use geno
    cond = wc.cond.rsplit('_', 1)[0]
    return nanopore(f'aligned_data/pooled/{cond}.genome.bam')


rule run_stringtie:
    input:
        illumina='aligned_data/pooled/{cond}.sorted.bam',
        nanopore=match_nanopore_input,
    output:
        gtf='assemblies/{cond}.gtf'
    params:
        annot=config['annot_gtf_fn'],
    threads:
        24
    conda:
        'env_yamls/stringtie.yaml'
    shell:
        '''
        stringtie -p {threads} --mix --rf -c 1 -s 1 -g 0 -M 10 \
          -G {params.annot} \
          -o {output.gtf} \
          {input.illumina} {input.nanopore}
        '''


def get_assembly_input(wc):
    sample_names = glob_wildcards('raw_data/{sample_name}.1.fastq.gz').sample_name
    conds = set([sn.rsplit('_', 1)[0] for sn in sample_names])
    return expand(
        'assemblies/{conds}.gtf',
        conds=conds
    )


rule stringtie_merge:
    input:
        get_assembly_input
    output:
        'assemblies/merged_stringtie_assembly.gtf'
    params:
        annot=config['annot_gtf_fn'],
    threads: 12
    conda:
        'env_yamls/stringtie.yaml'
    shell:
        '''
        stringtie --merge -g 0 -F 0 -T 0 -f 0.001 -i \
          -G {params.annot} \
          -o {output}.tmp.gtf {input}
        python ../scripts/relabel_stringtie_merge_gene_ids.py {output}.tmp.gtf {output}
        rm {output}.tmp.gtf
        '''


rule gtf_to_fasta:
    input:
        'assemblies/merged_stringtie_assembly.gtf'
    output:
        'assemblies/merged_stringtie_assembly.fa'
    params:
        genome_fasta_fn=config['genome_fasta_fn']
    conda:
        'env_yamls/stringtie.yaml'
    shell:
        '''
        gffread -w {output} -g {params.genome_fasta_fn} {input}
        '''


rule annotate_orfs:
    input:
        gtf='assemblies/merged_stringtie_assembly.gtf',
        fasta='assemblies/merged_stringtie_assembly.fa',
    output:
        gtf='assemblies/merged_stringtie_assembly.orf.gtf',
        orf_annot='assemblies/merged_stringtie_assembly.orf_annot.csv'
    conda:
        'env_yamls/transuite.yaml'
    shell:
        '''
        python ../scripts/TranSuite/transuite.py Auto \
          --cds 50 \
          --gtf {input.gtf} \
          --fasta {input.fasta} \
          --outpath assemblies \
          --outname orf \
          --ptc 90
        cp assemblies/orf_TranSuite_output/orf_transfix/orf_transfix.gtf {output.gtf}
        cp assemblies/orf_TranSuite_output/orf_transfeat/orf_transfeat.csv {output.orf_annot}
        '''