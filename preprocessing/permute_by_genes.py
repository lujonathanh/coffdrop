__author__ = 'jlu96'

"""Given empirical gene mutation frequencies,
    generate a random mutation matrix out of the existing one."""

import matrixreader as mr
import sys

# read in the file

def get_parser():
    # Parse arguments
    import argparse
    description = 'Given empirical gene mutation frequencies,    generate a random mutation matrix out of the existing one.'
    parser = argparse.ArgumentParser(description = description)


    parser.add_argument('-i', '--input_matrix', required=True)

    parser.add_argument('-r', '--random_seed', required=True, type=int)

    parser.add_argument('-o', '--output_name', required=True)

    return parser


def run(args):
    input = args.input_matrix
    output = args.output_name
    seed = args.random_seed

    matrixreader = mr.MatrixReader(input)

    newmr = matrixreader.permute_by_gene(seed)
    newmr.writemutationmatrix(output)

    print "Gene permuted matrix written to ", output


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
