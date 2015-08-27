__author__ = 'jlu96'
import maf
import sys
import time
import os

def get_parser():
    # Parse arguments
    import argparse
    description = 'Convert one maf file to one matrix file.'
    parser = argparse.ArgumentParser(description = description)

    # General parameters
    parser.add_argument('-i', '--input_maf', required=True, help="Input maf file. ")

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-l', '--level', default='patient', help='Last index of tumor barcode'
                                                                 'to limit to')

    #parser.add_argument('-v', '--variant_type', default='patient')

    parser.add_argument('-al', '--append_loss', default=0, help='If 1, append "loss" to the '
                                                              'end of each gene name.')

    return parser

def run(args):
    input = args.input_maf
    output = args.output_name
    level = args.level
    appendloss = args.append_loss

    t1 = time.time()
    mafreader = maf.MafReader()
    print "maf file size: ", os.path.getsize(input)

    mafreader.__init__(filename=input)
    t2 = time.time()
    print "number of records: ", mafreader.size
    print "time used to read in: ", t2 - t1

    try:
        eval(level)
        level = eval(level)
    except NameError:
        pass

    mafreader.writemutationmatrix(filename=output, level=level, gainloss=appendloss)
    t3 = time.time()

    print "Mutation matrix written to ", output

def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()