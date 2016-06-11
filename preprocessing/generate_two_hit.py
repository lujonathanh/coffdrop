__author__ = 'jlu96'



__author__ = 'jlu96'

"""Given a list of genes,
    generate a list of genes and their copy-numbers."""

import sys
import csv

# read in the file

def get_parser():
    # Parse arguments
    import argparse
    description = 'Given empirical gene mutation frequencies,    generate a random mutation matrix out of the existing one.'
    parser = argparse.ArgumentParser(description = description)


    parser.add_argument('-i', '--input_file', required=True)

    parser.add_argument('-o', '--output_file', required=True)

    return parser


def run(args):
    input = args.input_file
    output = args.output_file


    genes = set()

    with open(input, 'rU') as genefile:
        for row in genefile:
            genes.add(row.rstrip())
    print genes

    gene_list = list(genes)
    loss_gene_list = list([gene + 'loss' for gene in genes])

    new_gene_list = gene_list + loss_gene_list

    with open(output, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for gene in new_gene_list:
            writer.writerow([gene])

    print "Two-hit gene list written to ", output


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
