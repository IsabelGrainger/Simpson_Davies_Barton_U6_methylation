import os


def find_fast5s(sample_name):
    fast5_fns = {}
    sample_basedir = f'raw_data/{sample_name}'
    for subdir, _, fns in os.walk(sample_basedir, followlinks=True):
        for fn in fns:
            name, ext = os.path.splitext(fn)
            if ext == '.fast5':
                fast5_fns[name] = os.path.join(subdir, fn)
    return fast5_fns


def fast5_finder():
    sample_names = glob_wildcards('raw_data/{sample_name}').sample_name
    sample_fast5s = {}
    for sn in sample_names:
        sample_fast5s[sn] = find_fast5s(sn)
    def _fast5_finder(wc):
        return sample_fast5s[wc.sample_name][wc.fast5_name]
    return _fast5_finder


rule guppy_basecall:
    input:
        fast5_fn=fast5_finder()
    params:
        flowcell=config['flowcell'],
        kit=config['kit'],
    output:
        "basecalling/{sample_name}/{fast5_name}/{fast5_name}.complete"
    threads: 4
    resources:
        job_class='short'
    shell:
        '''
        GUPPY=$(readlink -f ../scripts/ont-guppy-cpu/bin/guppy_basecaller)
        OUTPUT_DIR=$(readlink -f basecalling/{wildcards.sample_name}/{wildcards.fast5_name})
        OUTPUT=$(readlink -f {output})
        cp {input.fast5_fn} $TMPDIR
        cd $TMPDIR
        $GUPPY \
          --flowcell {params.flowcell} \
          --kit {params.kit} \
          --num_callers {threads} \
          --cpu_threads_per_caller 4 \
          --records_per_fastq 0 \
          --reverse_sequence yes \
          --input_path $TMPDIR \
          --save_path $OUTPUT_DIR
        touch $OUTPUT
        '''


def get_all_basecalling_output(wc):
    fast5s = find_fast5s(wc.sample_name)
    dummy_outputs = []
    for fast5_name, fn in fast5s.items():
        dummy_outputs.append(f'basecalling/{wc.sample_name}/{fast5_name}/{fast5_name}.complete')
    return dummy_outputs


rule concatenate_fastqs:
    input:
        get_all_basecalling_output
    output:
        fastq="basecalled_data/{sample_name}.rna.fastq",
        seq_summary="sequencing_summaries/{sample_name}_sequencing_summary.txt"
    threads: 1
    shell:
        '''
        # add header to sequencing summary
        DUMMY="{input[0]}"
        HEADER_FILE="${{DUMMY%/*.complete}}/sequencing_summary.txt"
        head -n1 $HEADER_FILE > {output.seq_summary}

        for DUMMY in {input};
          do
          BASECALL_DIR="${{DUMMY%/*.complete}}"
          for FASTQ in "${{BASECALL_DIR}}/*.fastq";
            do
            cat $FASTQ >> {output.fastq}
          done
          SEQ_SUMMARY="${{BASECALL_DIR}}/sequencing_summary.txt"
          tail -n+2 $SEQ_SUMMARY >> {output.seq_summary}
        done
        '''


rule rna_to_dna:
    input:
        "basecalled_data/{sample_name}.rna.fastq"
    output:
        "basecalled_data/{sample_name}.dna.fastq"
    threads: 1
    conda:
        'env_yamls/seqkit.yaml'
    shell:
        "cat {input} | seqkit seq --rna2dna > {output}"


rule archive_fast5s:
    input:
        data='raw_data/{sample_name}'
    output:
        archive='fast5_archives/{sample_name}.tar.gz'
    shell:
        '''
        OUTPUT=$(readlink -f {output.archive})
        cp -L --parents `find raw_data/{wildcards.sample_name}/ -name '*.fast5'` $TMPDIR
        cd $TMPDIR
        tar -cvzf {wildcards.sample_name}.tar.gz {input.data}
        mv {wildcards.sample_name}.tar.gz $OUTPUT
        '''


rule md5_fast5:
    input:
        data='fast5_archives/{sample_name}.tar.gz'
    output:
        md5='fast5_archives/{sample_name}.tar.gz.md5'
    shell:
        '''
        md5sum {input.data} > {output.md5}
        '''