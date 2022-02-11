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
        bed='yanocomp/{comp,\w+}.bed',
        smpreds='yanocomp/{comp}.smpreds.json.gz'
    params:
        comp_flag=get_comp_flag
    threads: 12
    conda:
        'env_yamls/yanocomp.yaml'
    resources:
        job_class='long'
    shell:
        '''
        CURRDIR=$(pwd)
        cp --parents -L {input} $TMPDIR
        cd $TMPDIR
        yanocomp gmmtest \
          {params.comp_flag} \
          -p {threads} \
          -o "${{CURRDIR}}/{output.bed}" \
          -s "${{CURRDIR}}/{output.smpreds}"
        '''


def differr_input(wildcards):
    return {
        'cntrl_events': expand(
            'aligned_data/{sample_name}.genome.bam',
            sample_name=sample_name_subset(wildcards.cntrl)
        ),
        'treat_events': expand(
            'aligned_data/{sample_name}.genome.bam',
            sample_name=sample_name_subset(wildcards.treat)
        )
    }


rule differr:
    input:
        unpack(differr_input),
    output:
        bed='differr/{treat,\w+}_vs_{cntrl,\w+}.bed'
    params:
        cntrl_flag=lambda wc, input: '-b ' + ' -b '.join(input.cntrl_events),
        treat_flag=lambda wc, input: '-a ' + ' -a '.join(input.treat_events),
        reference=os.path.abspath(config['genome_fasta_fn']),
    threads: 48
    conda:
        'env_yamls/differr.yaml'
    resources:
        job_class='short'
    shell:
        '''
        CURRDIR=$(pwd)
        for BAM in {input}; do
           cp --parents -L $BAM $TMPDIR
           cp --parents -L "${{BAM}}.bai" $TMPDIR
        done
        cd $TMPDIR
        differr \
          {params.cntrl_flag} \
          {params.treat_flag} \
          -r {params.reference} \
          -p {threads} \
          -o "${{CURRDIR}}/{output.bed}"
        '''


rule motif_analysis:
    input:
        '{method}/{cond}.bed'
    output:
        meme=directory('{method}/motif_detection/{cond}.meme'),
        seqs='{method}/motif_detection/{cond}.seqs.fa'
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
        seqs='{method}/motif_detection/{cond}.seqs.fa',
        meme='{method}/motif_detection/{cond}.meme'
    output:
        fimo=directory('{method}/motif_detection/{cond}.fimo'),
        motifs='{method}/motif_detection/{cond,\w+}.motifs.bed',
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