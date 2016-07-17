__author__ = 'jlu96'
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

    parser.add_argument('-g', '--geneFile', default=None)

    parser.add_argument('-p', '--patientFile', default=None)

    parser.add_argument('-o', '--output_name', required=True)

    return parser


def run(args):
    input = args.input_matrix
    geneFile = args.geneFile
    patientFile = args.patientFile
    output = args.output_name


    matrixreader = mr.MatrixReader(input)

    if geneFile != None:
        matrixreader = matrixreader.filter(geneFile)
    if patientFile != None:
        matrixreader = matrixreader.filter_samples(patientFile)

    matrixreader.writemutationmatrix(output)

    print "Filtered matrix written to ", output


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
