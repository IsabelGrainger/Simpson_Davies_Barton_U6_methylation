import re
import sys


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError()
    return attr


def gtf_parser(gtf_fn):
    gtf_records = {}
    gtf_ref_genes = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            if record.startswith('#'):
                continue
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'transcript':
                transcript_id = get_gtf_attribute(record, 'transcript_id')
                locus_id = get_gtf_attribute(record, 'gene_id')
                try:
                    ref_gene_id = get_gtf_attribute(record, 'ref_gene_id')
                except ValueError:
                    ref_gene_id = None
                gtf_ref_genes[transcript_id] = ref_gene_id
                idx = (chrom, locus_id, strand)
                if idx not in gtf_records:
                    gtf_records[idx] = {}
            if feat_type == 'exon':
                transcript_id = get_gtf_attribute(record, 'transcript_id')
                locus_id = get_gtf_attribute(record, 'gene_id')
                idx = (chrom, locus_id, strand)
                if transcript_id not in gtf_records[idx]:
                    gtf_records[idx][transcript_id] = []
                gtf_records[idx][transcript_id].append([start, end])
    for (chrom, locus_id, strand), transcripts in gtf_records.items():
        ref_gene_ids = {t_id: gtf_ref_genes[t_id] for t_id in transcripts}
        yield chrom, locus_id, strand, transcripts, ref_gene_ids


def get_intron_set(exons):
    exons.sort()
    introns = {(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)}
    return introns


def reassign_gene_id_from_ref(gtf_fn):
    for chrom, locus_id, strand, transcripts, ref_gene_ids in gtf_parser(gtf_fn):
        unique_ref_gene_ids = set(ref_gene_ids.values())
        try:
            unique_ref_gene_ids.remove(None)
        except KeyError:
            pass
        if not len(unique_ref_gene_ids):
            # gene is novel, use locus_id
            reass_ref_gene_ids = {t_id: locus_id for t_id in ref_gene_ids}
        elif len(unique_ref_gene_ids) == 1:
            # locus is single gene
            new_ref_gene_id = unique_ref_gene_ids.pop()
            reass_ref_gene_ids = {t_id: new_ref_gene_id for t_id in ref_gene_ids}
        elif len(unique_ref_gene_ids) > 1:
            # locus is actually multiple genes
            ref_gene_introns = {}
            for t_id, exons in transcripts.items():
                t_ref_id = ref_gene_ids[t_id]
                if t_ref_id is not None:
                    if t_ref_id not in ref_gene_introns:
                        ref_gene_introns[t_ref_id] = set()
                    ref_gene_introns[t_ref_id].update(get_intron_set(exons))
            # now assign by intron intersection
            reass_ref_gene_ids = {}
            for t_id, exons in transcripts.items():
                if ref_gene_ids[t_id] is not None:
                    reass_ref_gene_ids[t_id] = ref_gene_ids[t_id]
                    continue
                introns = get_intron_set(exons)
                possible_ref_gene_ids = set()
                for r_g_id, r_introns in ref_gene_introns.items():
                    if not r_introns.isdisjoint(introns):
                        possible_ref_gene_ids.add(r_g_id)
                if not len(possible_ref_gene_ids):
                    reass_ref_gene_ids[t_id] = locus_id
                elif len(possible_ref_gene_ids) == 1:
                    reass_ref_gene_ids[t_id] = possible_ref_gene_ids.pop()
                else:
                    # assign by most upstream/downstream
                    get_most_upstream = min if strand == '+' else max
                    most_upstream = {r_g_id: get_most_upstream(ref_gene_introns[r_g_id])[0] for r_g_id in possible_ref_gene_ids}
                    reass_ref_gene_ids[t_id] = get_most_upstream(most_upstream, key=most_upstream.__getitem__)
        i = 1
        for transcript_id, exons in transcripts.items():
            gene_id = reass_ref_gene_ids[transcript_id]
            if not transcript_id.rsplit('.', 1)[0] == gene_id:
                transcript_id = f'{gene_id}.novel{i}'
                i += 1
            exons.sort()
            yield chrom, exons[0][0], exons[-1][1], strand, exons, gene_id, transcript_id


def write_gtf(record_iterator, output_gtf_fn):
    with open(output_gtf_fn, 'w') as gtf:
        for chrom, start, end, strand, exons, gene_id, transcript_id in record_iterator:
            attr = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
            gtf.write(f'{chrom}\ttransdecoder\ttranscript\t{start + 1}\t{end}\t.\t{strand}\t.\t{attr}\n')
            for exstart, exend in exons:
                gtf.write(f'{chrom}\ttransdecoder\texon\t{exstart + 1}\t{exend}\t.\t{strand}\t.\t{attr}\n')


if __name__ == '__main__':
    write_gtf(
        reassign_gene_id_from_ref(sys.argv[1]),
        sys.argv[2]
    )