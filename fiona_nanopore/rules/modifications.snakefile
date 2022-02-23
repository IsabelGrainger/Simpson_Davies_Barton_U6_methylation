import itertools as it

def sample_name_subset(cond):
    sample_names = glob_wildcards('raw_data/{sample_name}').sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    assert len(cond_sample_names), cond
    return cond_sample_names


def yanocomp_input(wildcards):
    conds = wildcards.comp.split('_vs_')
    return [
        expand(
            'nanopolish/{sample_name}.collapsed.h5',
            sample_name=sample_name_subset(c)
        ) for c in conds
    ]


def get_comp_flag(wildcards, input):
    sample_names = glob_wildcards('nanopolish/{sample_name}.collapsed.h5', files=input).sample_name
    conds = [sn.rsplit('_', 1)[0] for sn in sample_names]
    comp_flag = ' '.join(f'-i {fn} {c}' for fn, c in zip(input, conds))
    return comp_flag


rule yanocomp_gmmtest:
    input:
        unpack(yanocomp_input),
    output:
        bed='yanocomp/{comp,\w+}.bed'
    params:
        comp_flag=get_comp_flag
    threads: 12
    conda:
        'env_yamls/yanocomp.yaml'
    resources:
        job_class='long',
        hostname='c6[43]*'
    shell:
        '''
        CURRDIR=$(pwd)
        cp --parents -L {input} $TMPDIR
        cd $TMPDIR
        yanocomp gmmtest \
          {params.comp_flag} \
          -p {threads} \
          -o "${{CURRDIR}}/{output.bed}"
        '''


def yanocomp_upsetplot_input(wc):
    sites = expand(
        'yanocomp/{comp}.bed',
        comp=config.get('multicomp', [])
    )
    return {
        'yanocomp_sites': sites,
        'miclip_peaks': ancient(config['miclip_peaks']),
        'gtf': ancient(config['gtf_fn']),
        'fasta': ancient(config['genome_fasta_fn'])
    }


def get_multicomp_combinations():
    multicomp = config.get('multicomp', [])
    for mc in multicomp:
        conds = mc.split('_vs_')
        # assume last position is overall control
        cntrl = conds[-1]
        comps = [f'{c1}_vs_{c2}' for c1, c2 in it.combinations(conds, r=2)]
        comps_intersection = [f'{comp1}__and__{comp2}' for comp1, comp2 in it.permutations(comps, r=2)]
        comps_difference = [f'{comp1}__not__{comp2}' for comp1, comp2 in it.permutations(comps, r=2)]
        comps_difference_miclip = [f'{comp1}__not__{comp2}__miclip' for comp1, comp2 in it.permutations(comps, r=2)]
        yield from comps + comps_intersection + comps_difference + comps_difference_miclip


rule generate_yanocomp_upsetplots:
    input:
        unpack(yanocomp_upsetplot_input)
    output:
        sites=[f'yanocomp/{comp}.posthoc.bed' for comp in get_multicomp_combinations()],
        site_upset='figures/yanocomp/yanocomp_site_upset.svg',
        site_effect_size='figures/yanocomp/yanocomp_site_effect_size.svg',
    conda:
        'env_yamls/nb_seqlogos.yaml'
    log:
        notebook='notebook_processed/yanocomp_upsetplots.py.ipynb'
    notebook:
        'notebook_templates/yanocomp_upsetplots.py.ipynb'


rule motif_analysis:
    input:
        'yanocomp/{cond}.posthoc.bed'
    output:
        meme=directory('yanocomp/motif_detection/{cond}.meme'),
        seqs='yanocomp/motif_detection/{cond}.seqs.fa'
    params:
        fasta=config['genome_fasta_fn'],
        fai=config['genome_fasta_fn'] + '.fai',
        expected_motif=config.get('expected_motif', 'NNANN'),
        minw=len(config.get('expected_motif', 'NNANN')) - 2,
        maxw=len(config.get('expected_motif', 'NNANN')) + 2,
    conda:
        'env_yamls/meme.yml'
    shell:
        '''
        bedtools slop -b 5 \
          -i {input} \
          -g <(cut -f-2 {params.fai}) |
        bedtools merge -s -c 4,5,6 -o collapse,max,distinct \
          -i stdin |
        bedtools getfasta -s \
          -fi {params.fasta} \
          -fo {output.seqs} \
          -bed stdin
        meme -oc {output.meme} \
          -cons {params.expected_motif} \
          -dna -nmotifs 1 \
          -minw {params.minw} -maxw {params.maxw} \
          -mod zoops \
          {output.seqs}
        '''


def build_regex(wc):
    iupac_to_regex = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'S': '[GC]', 'W': '[AT]',
        'R': '[AG]', 'Y': '[CT]',
        'K': '[GT]', 'M': '[AC]',
        'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]',
        'N': '[ACGT]'
    }
    motif = config.get('expected_motif', 'NNANN').upper()
    regex = []
    for iupac in motif:
        regex.append(iupac_to_regex[iupac])
    return ''.join(regex)


rule find_motifs:
    input:
        seqs='yanocomp/motif_detection/{cond}.seqs.fa',
        meme='yanocomp/motif_detection/{cond}.meme'
    output:
        fimo=directory('yanocomp/motif_detection/{cond}.fimo'),
        motifs='yanocomp/motif_detection/{cond,\w+}.motifs.bed',
    params:
        regex=build_regex
    conda:
        'env_yamls/meme.yml'
    shell:
        '''
        fimo --norc --thresh 0.1 \
            -oc {output.fimo} \
            {input.meme}/meme.txt \
            {input.seqs}
        sort -k1,1 -k6,6nr {output.fimo}/fimo.gff |
        awk '!seen[$1]++' |
        awk -v OFS='\\t' '
            $0 !~ "^##" {{
                split($1,a,":");
                split(a[2],b,"(");
                split(b[1],pos,"-");
                split($9,attr,";");
                split(attr[6],seq_attr,"=");
                seq=seq_attr[2];
                if (seq ~ /{params.regex}/) {{
                    strand=substr(b[2], 1, 1);
                    if (strand == "+") {{
                        start=$4+pos[1];
                        end=$5+pos[1];
                    }}
                    else {{
                        start=pos[2]-($5-1);
                        end=pos[2]-($4-1);
                    }}
                    print a[1], start - 1, end, seq, $6, strand
                }}
            }}' > {output.motifs}
            
        '''
    

rule generate_yanocomp_logos:
    input:
        sites='yanocomp/{comp}.posthoc.bed',
        motifs='yanocomp/motif_detection/{comp}.motifs.bed',
        gtf=ancient(config['gtf_fn']),
        fasta=ancient(config['genome_fasta_fn'])
    output:
        logo_plot='figures/yanocomp/{comp}_yanocomp_motif_logo.svg',
        distrib_plot='figures/yanocomp/{comp}_yanocomp_distrib_barplot.svg',
    conda:
        'env_yamls/nb_seqlogos.yaml'
    log:
        notebook='notebook_processed/{comp}_yanocomp_logos.py.ipynb'
    notebook:
        'notebook_templates/yanocomp_logos.py.ipynb'


rule generate_m6a_gene_tracks:
    input:
        yanocomp=expand('yanocomp/{comp}.bed', comp=config['multicomp']),
        yanocomp_posthoc=expand('yanocomp/{comp}.posthoc.bed', comp=get_multicomp_combinations()),
        miclip_cov=ancient(config['miclip_coverage']),
        miclip_peaks=ancient(config['miclip_peaks']),
        der_sites=ancient(config['der_sites']),
        gtf=ancient(config['gtf_fn']),
    output:
        gene_track='figures/yanocomp/gene_tracks/{gene_id}_m6a_gene_track.svg'
    conda:
        'env_yamls/nb_seqlogos.yaml'
    log:
        notebook='notebook_processed/{gene_id}_m6a_gene_tracks.py.ipynb'
    notebook:
        'notebook_templates/m6a_gene_tracks.py.ipynb'