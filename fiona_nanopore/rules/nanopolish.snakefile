import os


rule nanopolish_index:
    input:
        fastq='basecalled_data/{sample_name}.dna.fastq',
        fast5='raw_data/{sample_name}',
        seq_summary='sequencing_summaries/{sample_name}_sequencing_summary.txt'
    output:
        "basecalled_data/{sample_name}.dna.fastq.index",
        "basecalled_data/{sample_name}.dna.fastq.index.fai",
        "basecalled_data/{sample_name}.dna.fastq.index.gzi",
        "basecalled_data/{sample_name}.dna.fastq.index.readdb"
    shell:
        '''
        ../scripts/nanopolish/nanopolish index \
        -s {input.seq_summary} \
        -d {input.fast5} \
        {input.fastq}
        '''


rule nanopolish_eventalign:
    input:
        fastq='basecalled_data/{sample_name}.dna.fastq',
        readdb="basecalled_data/{sample_name}.dna.fastq.index.readdb",
        bam='aligned_data/{sample_name}.transcriptome.bam',
        reference=os.path.abspath(config['transcriptome_fasta_fn']),
        gtf=os.path.abspath(config["gtf_fn"])
    output:
        events='nanopolish/{sample_name}.collapsed.h5',
    threads: 12
    conda:
        'env_yamls/yanocomp.yaml'
    resources:
        job_class='long',
        local_free='500G',
    shell:
        '''
        NANOPOLISH=$(readlink -f ../scripts/nanopolish/nanopolish)
        OUTPUT=$(readlink -f {output.events})
        #cp -L --parents `find basecalled_data/ -name '{wildcards.sample_name}.dna.fastq*'` $TMPDIR
        #cp -L --parents `find aligned_data/ -name '{wildcards.sample_name}.bam*'` $TMPDIR
        #cp -L --parents `find raw_data/{wildcards.sample_name}/ -name '*.fast5'` $TMPDIR
        #cd $TMPDIR
        $NANOPOLISH eventalign -t 8 \
          --scale-events --signal-index --print-read-names \
          -r {input.fastq} \
          -b {input.bam} \
          -g {input.reference} |
        yanocomp prep \
          {input.gtf} \
          -h $OUTPUT -p 4
        '''