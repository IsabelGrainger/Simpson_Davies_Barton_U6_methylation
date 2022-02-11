def salmon_index_input(wc):
    input_ = {
        'decoys_fasta_fn': ancient(config['genome_fasta_fn'])
    }
    if wc.ref_type == 'ref':
        input_['transcriptome_fasta_fn'] = ancient(config['transcriptome_fasta_fn'])
    elif wc.ref_type == 'denovo':
        input_['transcriptome_fasta_fn'] = 'assemblies/merged_stringtie_assembly.fa'
    return input_


rule build_salmon_index:
    input:
        unpack(salmon_index_input)
    output:
        directory('annotations/{ref_type}_salmon_index')
    conda:
        'env_yamls/salmon.yaml'
    threads: 8
    shell:
        '''
        cat {input.decoys_fasta_fn} | grep "^>" | cut -c2- | cut -d" " -f1 > decoys.txt
        cat {input.transcriptome_fasta_fn} {input.decoys_fasta_fn} > tmp.fa
        salmon index -p {threads} -i {output} -t tmp.fa -d decoys.txt
        rm tmp.fa decoys.txt
        '''


rule pseudoalign_with_salmon:
    input:
        read='raw_data/{sample_name}.1.fastq.gz',
        mate='raw_data/{sample_name}.2.fastq.gz',
        index='annotations/{ref_type}_salmon_index'
    output:
        'quantification/{ref_type}/{sample_name}/quant.sf'
    params:
        prefix=lambda wc, output: os.path.split(output[0])[0]
    conda:
        'env_yamls/salmon.yaml'
    threads: 8
    shell:
        '''
        salmon quant -l A -p {threads} \
          --validateMappings \
          -i {input.index} \
          -1 {input.read} -2 {input.mate} \
          -o {params.prefix}
        '''


rule run_edger:
    input:
        salmon_quants=expand(
            'quantification/{{ref_type}}/{sample_name}/quant.sf',
            sample_name=glob_wildcards('raw_data/{sample_name}.1.fastq.gz').sample_name
        )
    output:
        cpm='quantification/{ref_type}/cpm.csv',
        dge='quantification/{ref_type}/dge.csv',
        pca_plot='figures/dge/{ref_type}/pca.svg',
        fio1_plot='figures/dge/{ref_type}/fio1_expression.svg',
        mafs_plot='figures/dge/{ref_type}/mafs_expression.svg',
        ft_plot='figures/dge/{ref_type}/ft_genes_expression.svg',
    conda:
        'env_yamls/nb_rpy2_edger.yaml'
    log:
        notebook='notebook_processed/edger_{ref_type}.py.ipynb'
    notebook:
        'notebook_templates/edger.py.ipynb'