import sys
import csv


def get_transfeat_classes(transfeat_fn):
    nmd_class = {}
    with open(transfeat_fn) as f:
        reader = csv.DictReader(f)
        for record in reader:
            gene_id = record['Gene_ID']
            transcript_id = record['Transcript_ID']
            classi = record['Coding_potentiality'] != "Coding"
            chrom = record['Transcript_coordinates'].split(':')[0]
            idx = (chrom, gene_id)
            if idx not in nmd_class:
                nmd_class[idx] = {True: [], False: []}
            nmd_class[idx][classi].append(transcript_id)
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
    nmd_class = get_transfeat_classes(
        sys.argv[1]
    )
    to_ioe(nmd_class, sys.argv[2])