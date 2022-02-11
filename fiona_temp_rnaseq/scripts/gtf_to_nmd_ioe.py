import re
import sys
from collections import defaultdict
import numpy as np
import pysam
from orfipy_core import orfs as find_orfs


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


def gtf_parser(gtf_fn):
    gtf_records = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            if record.startswith('#'):
                continue
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'exon' or feat_type == 'CDS':
                transcript_id = get_gtf_attribute(record, 'transcript_id')
                gene_id = get_gtf_attribute(record, 'gene_id')
                idx = (chrom, gene_id, transcript_id, strand)
                if idx not in gtf_records:
                    gtf_records[idx] = {}
                if feat_type not in gtf_records[idx]:
                    gtf_records[idx][feat_type] = []
                gtf_records[idx][feat_type].append([start, end])
    for (chrom, gene_id, transcript_id, strand), features in  gtf_records.items():
        features['exon'].sort()
        exons = np.array(features['exon'])
        try:
            features['CDS'].sort()
            cds = np.array(features['CDS'])
        except KeyError:
            cds=None
        yield chrom, gene_id, transcript_id, strand, exons, cds


RC = str.maketrans('ACGT', 'TGCA')


def rev_comp(seq):
    return seq.translate(RC)[::-1]


def get_record_sequence(chrom, strand, exons, fasta):
    seq = []
    for start, end, in exons:
        seq.append(fasta.fetch(chrom, start, end))
    seq = ''.join(seq)
    if strand == '-':
        seq = rev_comp(seq)
    return seq


def genomic_cds_coords_to_transcriptomic(exons, cds, strand):
    ex_len = [e - s for s, e, in exons]
    transcript_len = sum(ex_len)
    ex_cumstart = np.cumsum([0] + ex_len[:-1])
    ex_starts = exons[:, 0]
    cds_first_ex = np.searchsorted(ex_starts, cds[0, 0], side='right') - 1
    cds_start_transcriptomic = ex_cumstart[cds_first_ex] + (cds[0, 0] - ex_starts[cds_first_ex])
    cds_last_ex = np.searchsorted(ex_starts, cds[-1, 0], side='right') - 1
    cds_end_transcriptomic = ex_cumstart[cds_last_ex] + (cds[-1, 1] - ex_starts[cds_last_ex])
    if strand == '-':
        cds_start_transcriptomic, cds_end_transcriptomic = (
            transcript_len - cds_end_transcriptomic,
            transcript_len - cds_start_transcriptomic
        )
    return cds_start_transcriptomic, cds_end_transcriptomic


def has_uorf(seq, dorf_start, min_size=50):
    orfs = find_orfs(
        seq, minlen=min_size * 3,
        starts=['ATG'],
        stops=['TAA', 'TAG', 'TGA'],
        partial3=False,
        partial5=False,
        include_stop=True,
        strand='f'
    )   
    for start, end, *_ in orfs:
        if start < dorf_start:
            return True
    else:
        return False


def stop_codon_to_ejc_distance(exons, cds, strand):
    if strand == '+':
        stop = cds[-1, 1]
        ejc_idx = np.searchsorted(exons[:, 1], stop, side='left')
        if ejc_idx == len(exons) - 1:
            return np.inf
        else:
            return exons[ejc_idx, 1] - stop
    elif strand == '-':
        stop = cds[0, 0]
        ejc_idx = np.searchsorted(exons[:, 0], stop, side='right')
        if ejc_idx == 1:
            return np.inf
        else:
            return stop - exons[ejc_idx - 1, 0]


def ptc_near_ejc(exons, cds, strand, min_dist=50):
    if stop_codon_to_ejc_distance(exons, cds, strand) <= min_dist:
        return True
    else:
        return False


def classify_nmd_transcripts(gtf_fn, fasta_fn, min_uorf_size=50, min_dist_to_ejc=50):
    nmd_class = {}
    with pysam.FastaFile(fasta_fn) as fasta:
        for chrom, gene_id, transcript_id, strand, exons, cds in gtf_parser(gtf_fn):
            idx = (chrom, gene_id)
            if idx not in nmd_class:
                nmd_class[idx] = defaultdict(list)
            seq = get_record_sequence(chrom, strand, exons, fasta)
            if cds is None:
                # assume NMD (may be some ncRNAs but we are not so interested in those)
                nmd_class[idx][True].append(transcript_id)
            else:
                cds_start, cds_end = genomic_cds_coords_to_transcriptomic(exons, cds, strand)
                if has_uorf(seq, cds_start, min_uorf_size):
                    nmd_class[idx][True].append(transcript_id)
                elif ptc_near_ejc(exons, cds, strand, min_dist_to_ejc):
                    nmd_class[idx][True].append(transcript_id)
                else:
                    nmd_class[idx][False].append(transcript_id)
    return nmd_class


def to_ioe(nmd_class, ioe_fn):
    with open(ioe_fn, 'w') as ioe:
        ioe.write('seqname\tgene_id\tevent_id\talternative_transcripts\ttotal_transcripts\n')
        for (chrom, gene_id), transcript_classes in nmd_class.items():
            if len(transcript_classes[False]) and len(transcript_classes[True]):
                incl_transcripts = ','.join(transcript_classes[False]) # non NMD transcripts
                total_transcripts = ','.join(transcript_classes[False] + transcript_classes[True])
                ioe.write('{chrom}\t{gene_id}\t{gene_id}\t{incl}\t{total}\n'.format(
                    chrom=chrom, gene_id=gene_id, incl=incl_transcripts, total=total_transcripts
                ))


if __name__ == '__main__':
    nmd_class = classify_nmd_transcripts(
        sys.argv[1], sys.argv[2]
    )
    to_ioe(nmd_class, sys.argv[3])