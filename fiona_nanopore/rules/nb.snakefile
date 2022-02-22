rule notebook_to_markdown:
    input:
        'notebook_processed/{notebook_name}.py.ipynb'
    output:
        'markdown_reports/{notebook_name}.py.md'
    params:
        output_dir='markdown_reports'
    conda:
        'env_yamls/jupyter.yaml'
    shell:
        '''
        jupyter-nbconvert --no-input --output-dir {params.output_dir} --to markdown {input}
        '''