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

    return parser


def run(args):
    input = args.input_matrix
    output = args.output_name


    matrixreader = mr.MatrixReader(input)

    matrixreader.write_weSME_mut_list(filename=output)

    print "WeSME mut_list written to ", output


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
