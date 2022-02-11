import os
import re


rule fastqc:
    input:
        'raw_data/{sample_name}.{read}.fastq.gz'
    output:
        'qc/{sample_name}.{read}_fastqc.html'
    conda:
        'env_yamls/fastqc.yaml'
    shell:
        '''
        fastqc -o qc {input}
        '''


rule bwa_index:
    input:
        lambda wc: config['fastq_screen_genomes'][wc.genome]
    output:
        'annotations/fastq_screen/{genome}.bwt'
    conda:
        'env_yamls/fastq_screen.yaml'
    shell:
        '''
        bwa index -p annotations/fastq_screen/{wildcards.genome} {input}
        '''


rule fastq_screen_config:
    input:
        expand(
            'annotations/fastq_screen/{genome}.bwt',
            genome=config['fastq_screen_genomes']
        )
    output:
        'annotations/fastq_screen/fastq_screen.config'
    run:
        with open(output[0], 'w') as f:
            for fasta in input:
                genome = re.search('^annotations/fastq_screen/(\w+).bwt$', fasta).group(1)
                prefix = os.path.splitext(fasta)[0]
                f.write(f'DATABASE\t{genome}\t{prefix}\tBWA\n')
            


rule fastq_screen:
    input:
        read='raw_data/{sample_name}.{read}.fastq.gz',
        config='annotations/fastq_screen/fastq_screen.config'
    output:
        txt='qc/{sample_name}.{read}_screen.txt',
        png='qc/{sample_name}.{read}_screen.html'
    params:
        subset=100000,
        aligner='bwa'
    conda:
        'env_yamls/fastq_screen.yaml'
    threads: 8
    shell:
        '''
        fastq_screen --outdir qc \
          --force \
          --aligner {params.aligner} \
          --conf {input.config} \
          --subset {params.subset} \
          --threads {threads} \
          {input.read}
        '''


def multiqc_input(wildcards):

    sample_names = glob_wildcards(
        'raw_data/{sample_name}.1.fastq.gz'
    ).sample_name

    return expand(
        ['qc/{sample_name}.{read}_fastqc.html',
         'qc/{sample_name}.{read}_screen.txt',
         'quantification/ref/{sample_name}/quant.sf'],
        sample_name=sample_names,
        read=[1, 2]
    )
        


rule multiqc:
    input:
        multiqc_input
    output:
        'qc/multiqc_report.html'
    conda:
        'env_yamls/multiqc.yaml'
    shell:
        '''
        multiqc -f -dd 2 -o qc qc quantification/ref
        '''