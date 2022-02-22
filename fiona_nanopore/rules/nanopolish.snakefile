import os


rule f5c_index:
    input:
        fastq='basecalled_data/{sample_name}.dna.fastq',
        fast5='raw_data/{sample_name}',
    output:
        "basecalled_data/{sample_name}.dna.fastq.index",
        "basecalled_data/{sample_name}.dna.fastq.index.fai",
        "basecalled_data/{sample_name}.dna.fastq.index.gzi",
        "basecalled_data/{sample_name}.dna.fastq.index.readdb"
    threads: 12
    shell:
        '''
        ../scripts/f5c-v0.7/f5c_x86_64_linux index --iop 12 \
        -d {input.fast5} \
        {input.fastq}
        '''


rule f5c_eventalign:
    input:
        fastq='basecalled_data/{sample_name}.dna.fastq',
        readdb="basecalled_data/{sample_name}.dna.fastq.index.readdb",
        bam='aligned_data/{sample_name}.transcriptome.bam',
        reference=ancient(os.path.abspath(config['transcriptome_fasta_fn'])),
        gtf=ancient(os.path.abspath(config['gtf_fn'])),
    output:
        events='nanopolish/{sample_name}.collapsed.h5',
    threads: 32
    conda:
        'env_yamls/yanocomp.yaml'
    resources:
        hostname="c6[34]20*",
        job_class="long",
    shell:
        '''
        F5C=$(readlink -f ../scripts/f5c-v0.7/f5c_x86_64_linux)
        OUTPUT=$(readlink -f {output.events})
        cp -L --parents `find basecalled_data/ -name '{wildcards.sample_name}.dna.fastq*'` $TMPDIR
        cp -L --parents `find aligned_data/ -name '{wildcards.sample_name}.transcriptome.bam*'` $TMPDIR
        cp -L --parents `find raw_data/{wildcards.sample_name}/ -name '*.fast5'` $TMPDIR
        cd $TMPDIR
        $F5C eventalign --rna -t 14 --iop 14 \
          --scale-events --signal-index --print-read-names \
          -r {input.fastq} \
          -b {input.bam} \
          -g {input.reference} |
        yanocomp prep \
          -g {input.gtf} \
          -h $OUTPUT -p 4
        '''