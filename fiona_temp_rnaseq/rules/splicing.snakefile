import os

def suppa_index_input(wc):
    if wc.ref_type == 'ref':
        return ancient(config['annot_gtf_fn'])
    elif wc.ref_type == 'denovo':
        return 'assemblies/merged_stringtie_assembly.gtf'
    return input_


rule suppa_index:
    input:
        suppa_index_input
    output:
        expand(
            'annotations/suppa/{{ref_type}}/stringtie_assembly_{event_type}_strict.{file_type}',
            event_type=['A3', 'A5', 'MX', 'RI', 'SE', 'AF', 'AL'],
            file_type=['ioe', 'gtf'],
        )
    params:
        output_prefix=lambda wc: f'annotations/suppa/{wc.ref_type}/stringtie_assembly'
    conda:
        'env_yamls/suppa.yaml'
    shell:
        '''
        suppa.py generateEvents \
          -f ioe --pool-genes \
          -e SE SS MX RI FL \
          -i {input} \
          -o {params.output_prefix}
        '''


rule suppa_nmd_index:
    input:
        'assemblies/merged_stringtie_assembly.orf_annot.csv'
    output:
        'annotations/suppa/denovo/stringtie_assembly_NMD_strict.ioe'
    params:
        genome_fasta_fn=config['genome_fasta_fn']
    conda:
        'env_yamls/suppa.yaml'
    shell:
        '''
        python ../scripts/transfeat_to_ioe.py \
          {input} {output}
        '''


rule rename_suppa_gtf:
    input:
        'annotations/suppa/{ref_type}/stringtie_assembly_{event_type}_strict.gtf'
    output:
        gtf='annotations/suppa/{ref_type}/stringtie_assembly_{event_type}_friendly.gtf',
        tsv='annotations/suppa/{ref_type}/stringtie_assembly_{event_type}_friendly.tsv',
    params:
        output_prefix=lambda wc, output: os.path.splitext(output.gtf)[0]
    conda:
        'env_yamls/suppa.yaml'
    shell:
        '''
        python ../scripts/rename_suppa_gtf.py -g {input} -o {params.output_prefix}
        '''


def sample_name_subset(cond):
    sample_names = glob_wildcards('raw_data/{sample_name}.1.fastq.gz').sample_name
    return [sn for sn in sample_names if sn.rsplit('_', 1)[0] == cond]


rule suppa_agg_counts:
    input:
        lambda wc: expand(
            'quantification/{{ref_type}}/{sample_name}/quant.sf',
            sample_name=sample_name_subset(wc.cond)
        )
    output:
        'quantification/{ref_type}/{cond}.quant.sf'
    conda:
        'env_yamls/suppa.yaml'
    shell:
        '''
        multipleFieldSelection.py -k 1 -f 4 \
          -i {input} \
          -o {output}
        '''


rule suppa_psi:
    input:
        annot='annotations/suppa/{ref_type}/stringtie_assembly_{event_type}_strict.ioe',
        counts='quantification/{ref_type}/{cond}.quant.sf',
    output:
        psi='splicing/{ref_type}/{cond}_{event_type}.psi'
    params:
        output_prefix=lambda wc, output: os.path.splitext(output.psi)[0]
    conda:
        'env_yamls/suppa.yaml'
    shell:
        '''
        suppa.py psiPerEvent --total-filter 1 \
          -i {input.annot} -e {input.counts} \
          -o {params.output_prefix}
        '''


rule suppa_diffsplice:
    input:
        annot='annotations/suppa/{ref_type}/stringtie_assembly_{event_type}_strict.ioe',
        cntrl_counts='quantification/{ref_type}/{cntrl}.quant.sf',
        treat_counts='quantification/{ref_type}/{treat}.quant.sf',
        cntrl_psi='splicing/{ref_type}/{cntrl}_{event_type}.psi',
        treat_psi='splicing/{ref_type}/{treat}_{event_type}.psi'
    output:
        dpsi='splicing/{ref_type}/{treat}_vs_{cntrl}_{event_type}.dpsi'
    params:
        output_prefix=lambda wc, output: os.path.splitext(output.dpsi)[0]
    conda:
        'env_yamls/suppa.yaml'
    shell:
        '''
        suppa.py diffSplice -m empirical -gc \
          -i {input.annot} \
          -p {input.treat_psi} {input.cntrl_psi} \
          -e {input.treat_counts} {input.cntrl_counts} \
          -o {params.output_prefix}
        '''


def linear_model_psi_input(wc):
    sample_names = glob_wildcards('raw_data/{sample_name}.1.fastq.gz').sample_name
    conds = set([sn.rsplit('_', 1)[0] for sn in sample_names])
    event_types = ['SE', 'A3', 'A5', 'MX', 'RI', 'AF', 'AL']
    psi = expand(
        'splicing/{{ref_type}}/{cond}_{event_type}.psi',
        cond=conds, event_type=event_types
    )
    friendly_codes = expand(
        'annotations/suppa/{ref_type}/stringtie_assembly_{event_type}_friendly.tsv',
        ref_type=wc.ref_type, event_type=event_types
    )
    return {'psi': psi, 'friendly_codes': friendly_codes}
    


rule linear_model_psi:
    input:
        unpack(linear_model_psi_input)
    output:
        psi='splicing/{ref_type}/psi.csv',
        psi_fit='splicing/{ref_type}/psi_fit.csv',
        event_type_temp_plot='figures/splicing/{ref_type}/temp_event_types.svg',
        event_type_geno_plot='figures/splicing/{ref_type}/geno_event_types.svg',
        gxt_ddpsi_histogram='figures/splicing/{ref_type}/gxt_ddpsi_histogram.svg',
        splice_event_upset='figures/splicing/{ref_type}/splice_event_upset.pdf',
    conda:
        'env_yamls/nb_statsmodels.yaml'
    threads: 16
    log:
        notebook='notebook_processed/psi_fit_{ref_type}.py.ipynb'
    notebook:
        'notebook_templates/psi_fit.py.ipynb'


rule generate_logos_and_heatmaps:
    input:
        psi_fit='splicing/{ref_type}/psi_fit.csv',
        genome_fasta_fn=ancient(config['genome_fasta_fn'])
    output:
        splice_plots=directory('figures/splicing/{ref_type}/{event_type}_logos_and_heatmaps/'),
    conda:
        'env_yamls/nb_seqlogos.yaml'
    log:
        notebook='notebook_processed/{event_type}_sequence_logos_{ref_type}.py.ipynb'
    notebook:
        'notebook_templates/{wildcards.event_type}_sequence_logos.py.ipynb'


def get_gene_track_input(wc):
    sample_names = glob_wildcards('raw_data/{sample_name}.1.fastq.gz').sample_name
    conds = set([sn.rsplit('_', 1)[0] for sn in sample_names])
    bws = expand('coverage_tracks/{cond}.{strand}.bw', cond=conds, strand=['fwd', 'rev'])
    return {
        'psi': f'splicing/{wc.ref_type}/psi.csv',
        'psi_fit': f'splicing/{wc.ref_type}/psi_fit.csv',
        'bws': bws,
        'fasta': ancient(config['genome_fasta_fn'])
    }


rule generate_gene_tracks:
    input:
        unpack(get_gene_track_input)
    output:
        gene_tracks=directory('figures/splicing/gene_tracks/{gene_id}_gene_tracks')
    conda:
        'env_yamls/nb_seqlogos.yaml'
    log:
        notebook='notebook_processed/{gene_id}_splicing_gene_tracks.py.ipynb'
    notebook:
        'notebook_templates/splicing_gene_tracks.py.ipynb'