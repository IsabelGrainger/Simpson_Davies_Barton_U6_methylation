rule build_STAR_index:
    '''Create the index required for alignment with STAR'''
    output:
        directory('annotations/STAR_index')
    log:
        'logs/STAR_idx.log'
    threads: 24
    params:
        fasta_fn=config['genome_fasta_fn'],
        gtf_fn=config['annot_gtf_fn'],
        overhang=149,
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        mkdir {output};
        STAR \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles {params.fasta_fn} \
          --sjdbGTFfile {params.gtf_fn} \
          --sjdbOverhang {params.overhang}
        mv Log.out {log}
        '''


rule map_with_STAR:
    '''map reads with STAR spliced aligner'''
    input:
        read='raw_data/{sample_name}.1.fastq.gz',
        mate='raw_data/{sample_name}.2.fastq.gz',
        index='annotations/STAR_index'
    output:
        bam='aligned_data/{sample_name}.sorted.bam',
        bai='aligned_data/{sample_name}.sorted.bam.bai',
        stats='aligned_data/{sample_name}.sorted.bamstats',
        sjdb='aligned_data/{sample_name}.sjdb.tsv',
        chimeric='aligned_data/{sample_name}.chimeric_aln.tsv'
    log:
        progress='logs/{sample_name}.STAR_progress.log',
        final='logs/{sample_name}.STAR_final.log',
        main='logs/{sample_name}.STAR.log'
    threads: 28
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        TOPDIR=$(pwd)
        STAR_TMP_DIR="aligned_data/{wildcards.sample_name}.tmpdir"
        mkdir -p $STAR_TMP_DIR
        cd $STAR_TMP_DIR
        STAR \
          --runThreadN {threads} \
          --genomeDir $TOPDIR/{input.index} \
          --readFilesIn $TOPDIR/{input.read} $TOPDIR/{input.mate} \
          --readFilesCommand "zcat" \
          --outFilterMultimapNmax 3 \
          --alignSJoverhangMin 8 \
          --alignSJDBoverhangMin 3 \
          --outFilterMismatchNmax 4 \
          --alignIntronMin 60 \
          --alignIntronMax 20000 \
          --outSAMtype BAM Unsorted \
          --chimOutType Junctions \
          --chimSegmentMin 15 \
          --chimScoreJunctionNonGTAG 0 \
          --chimSegmentReadGapMax 20000

        cd $TOPDIR
        mv ${{STAR_TMP_DIR}}/Log.progress.out {log.progress}
        mv ${{STAR_TMP_DIR}}/Log.final.out {log.final}
        mv ${{STAR_TMP_DIR}}/Log.out {log.main}

        mv ${{STAR_TMP_DIR}}/SJ.out.tab {output.sjdb}
        mv ${{STAR_TMP_DIR}}/Chimeric.out.junction {output.chimeric}
        samtools sort \
          -m 2G -@ {threads} \
          -o {output.bam} \
          "${{STAR_TMP_DIR}}/Aligned.out.bam"
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.stats}
        rm -rf $STAR_TMP_DIR
        '''


rule split_strand:
    input:
        bam='aligned_data/{sample_name}.sorted.bam',
        bai='aligned_data/{sample_name}.sorted.bam.bai'
    output:
        bam=temp('aligned_data/{sample_name}.sorted.{strand}.bam'),
        bai=temp('aligned_data/{sample_name}.sorted.{strand}.bam.bai')
    params:
        samflags_1=lambda wc: '-f 128 -F 16' if wc.strand == 'fwd' else '-f 144',
        samflags_2=lambda wc: '-f 80' if wc.strand == 'fwd' else '-f 64 -F 16'
    threads: 4
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        samtools view -@ {threads} -b {params.samflags_1} {input.bam} > {output.bam}.1.bam
        samtools index -@ {threads} {output.bam}.1.bam
        samtools view -@ {threads} -b {params.samflags_2} {input.bam} > {output.bam}.2.bam
        samtools index -@ {threads} {output.bam}.2.bam
        samtools merge -@ {threads} {output.bam} {output.bam}.1.bam {output.bam}.2.bam
        samtools index -@ {threads} {output.bam}
        rm {output.bam}.[12].bam
        rm {output.bam}.[12].bam.bai
        '''


def sample_name_subset(cond):
    sample_names = glob_wildcards(
        'raw_data/{sample_name}.1.fastq.gz'
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


rule pooled_genome_coverage:
    input:
        bams=lambda wc: expand(
            'aligned_data/{sample_name}.sorted.{strand}.bam',
            sample_name=sample_name_subset(wc.cond),
            strand=wc.strand,
        )
    output:
        bw='coverage_tracks/{cond}.{strand}.bw'
    params:
        fasta_fn=config['genome_fasta_fn']
    conda:
        'env_yamls/star.yaml'
    threads: 12
    shell:
        '''
        samtools merge -f -@ {threads} {output.bw}.tmp.bam {input.bams}
        samtools index {output.bw}.tmp.bam
        samtools depth -d0 {output.bw}.tmp.bam | 
          awk -v OFS='\t' '{{print $1, $2-1, $2, $3}}' > {output}.tmp.bdg
        bedGraphToBigWig {output}.tmp.bdg <(cut -f-2 {params.fasta_fn}.fai) {output}
        rm {output}.tmp.bdg {output.bw}.tmp.bam {output.bw}.tmp.bam.bai
        '''


rule pooled_alignments:
    input:
        bams=lambda wc: expand(
            'aligned_data/{sample_name}.sorted.bam',
            sample_name=sample_name_subset(wc.cond),
        )
    output:
        bam='aligned_data/pooled/{cond}.sorted.bam',
        bai='aligned_data/pooled/{cond}.sorted.bam.bai'
    conda:
        'env_yamls/star.yaml'
    threads: 12
    shell:
        '''
        samtools merge -f -@ {threads} {output.bam} {input.bams}
        samtools index {output.bam}
        '''