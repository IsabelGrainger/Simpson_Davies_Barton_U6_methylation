import re
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-g", "--gtf-file", help="Input gtf file",
                    required=True)
parser.add_argument("-o", "--output-prefix", required=True,
                    help="Prefix of the ouput file")

def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search('{} "(.+?)";'.format(attribute), gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError('gtf does not contain attribute {}'.format(attribute))
    return attr


def incrementing_index(index, identifier, fstring):
    try:
        idx = index[identifier]
    except KeyError:
        idx = fstring.format(idx=len(index))
        index[identifier] = idx
    return idx


def main():
    args = parser.parse_args()

    new_event_mapping = {
        'A5': {}, 'A3': {}, 'SE': {}, 'RI': {}, 'MX': {}, 'AF': {}, 'AL': {}
    }
    output_gtf = '{}.gtf'.format(args.output_prefix)
    output_mapping = '{}.tsv'.format(args.output_prefix)
    with open(args.gtf_file) as gtf, open(output_gtf, 'w') as o1:
        for line in gtf:
            gtf_record = line.split('\t')
            if len(gtf_record) != 9 and line.startswith("track name"):
                o1.write(line)
                continue
            event_id = get_gtf_attribute(gtf_record, 'gene_id')
            event_type = event_id.split(':', 1)[0]
            friendly_id = incrementing_index(new_event_mapping[event_type], event_id, '{}.{{idx}}'.format(event_type))
            line = re.sub(
                'transcript_id ".+?:alternative([12])"',
                'transcript_id "{}.alt\\1"'.format(friendly_id),
                line
            )
            o1.write(line)
    with open(output_mapping, 'w') as o2:
        for _, mapping in new_event_mapping.items():
            for event_id, friendly_id in mapping.items():
                o2.write('{}\t{}\n'.format(event_id, friendly_id))


if __name__ == '__main__':
    main()