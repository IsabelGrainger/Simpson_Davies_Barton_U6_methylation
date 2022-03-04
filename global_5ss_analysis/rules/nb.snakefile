rule suppa_index:
    input:
        gtf=lambda wc: ancient(config['genomes'][wc.organism]['gtf_fn']),
    output:
        'annotations/{organism}_A5_strict.gtf',
    params:
        output_prefix=lambda wc: f'annotations/{wc.organism}'
    conda:
        'env_yamls/suppa.yaml'
    shell:
        '''
        suppa.py generateEvents \
          -f ioe --pool-genes \
          -e SS \
          -i {input} \
          -o {params.output_prefix}
        '''


rule generate_splice_site_logos:
    input:
        fasta=lambda wc: ancient(config['genomes'][wc.organism]['genome_fasta_fn']),
        gtf=lambda wc: ancient(config['genomes'][wc.organism]['gtf_fn']),
        a5_gtf='annotations/{organism}_A5_strict.gtf',
    output:
        gurag_logo='figures/{organism}/{organism}_gurag_seqlogos.svg',
        aggu_logo='figures/{organism}/{organism}_aggu_seqlogos.svg',
        u5_u6_ecdfs='figures/{organism}/{organism}_u5_u6_ecdfs.svg',
        u5_u6_altsplice_scatter='figures/{organism}/{organism}_u5_u6_altsplice_scatter.svg',
    conda:
        'env_yamls/nb_seqlogos.yaml'
    log:
        notebook='notebook_processed/{organism}_u5_u6_splice_site_logos.py.ipynb'
    notebook:
        'notebook_templates/u5_u6_splice_site_logos.py.ipynb'