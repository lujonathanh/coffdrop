__author__ = 'jlu96'

import csv
import sys


def get_genes(segment):
    genes = segment[:-4].split('_')
    return set(genes)

def get_parser():

    import argparse

    description = 'Find and replace every instance of one file with another'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-c', '--COSMICFILE', default = "COSMICCensus.txt",
                        help='File name with conversions.')

    parser.add_argument('-ic', '--input_ID', default='Name')
    parser.add_argument('-oc', '--output_column', default='COSMIC_GENES')
    parser.add_argument('-i', '--annotate_file', required=True)
    parser.add_argument('-o', '--output_file', default=None)

    return parser


def run(args):

    # Load the parameters

    # Gene list
    COSMICFILE = args.COSMICFILE

    with open(COSMICFILE, 'rU') as f:
        cosmic_gene_list = set(l.rstrip().split()[0] for l in f if not l.startswith("#") and l)



    annotate_file = args.annotate_file
    genecol = args.input_ID
    output_file = args.output_file
    if not output_file:
        output_file = annotate_file + '_COSMIC.txt'

    new_col = args.output_column

    records = []

    with open(annotate_file, 'rU') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        header = reader.fieldnames
        for row in reader:
            segment = row[genecol]
            genes = get_genes(segment)
            # print genes
            cosmic_genes = genes.intersection(cosmic_gene_list)
            if cosmic_genes:
                row[new_col] = '_'.join(list(cosmic_genes)) + segment[-4:]
            else:
                row[new_col] = 'NA' + segment[-4:]
            records.append(row)


    new_header = header + [new_col]
    with open(output_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=new_header)
        writer.writeheader()
        for row in records:
            writer.writerow(row)

def main():
    run(get_parser().parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()