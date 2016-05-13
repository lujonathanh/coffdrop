__author__ = 'jlu96'

import matrixreader as mr
import sys

# read in the file

def get_parser():
    # Parse arguments
    import argparse
    description = 'Combine several matrix files to one matrix file.'
    parser = argparse.ArgumentParser(description = description)


    parser.add_argument('-i', '--input_matrix', required=True)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-e', '--ending', default="gain", help="Genes that end in this will be removed")

    return parser


def run(args):
    input = args.input_matrix
    output = args.output_name
    ending = args.ending

    matrixreader = mr.MatrixReader(input)

    newmr = matrixreader.removeEndingIn(ending=ending, length=len(ending))

    newmr.writemutationmatrix(output)

    print "Filtered matrix written to ", output


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()


# load the mutation matrix

# parse the args. if any end in loss..

# write to new mutation matrix