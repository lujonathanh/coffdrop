__author__ = 'jlu96'

import edgereader as edg
import sys


def get_parser():
    # Parse arguments
    import argparse
    description = 'Convert old edges to new edges.'
    parser = argparse.ArgumentParser(description = description)

    # General parameters
    parser.add_argument('-i', '--input', required=True, help="File to limit")

    parser.add_argument('-c', '--column', default='Probability')

    parser.add_argument('-t', '--threshold', type=float, default=0.05)

    parser.add_argument('-lt', '--less_than', type=int, default=1)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-gcn', '--gene_column_names', nargs='+', default=frozenset(['Gene0', 'Gene1']))


    return parser

def run(args):
    edgefile = args.input
    gene_column_names = args.gene_column_names
    column = args.column
    output = args.output_name
    threshold = args.threshold
    less_than = args.less_than

    Reader = edg.EdgeReader(edgefile, gene_column_names=gene_column_names)

    # Limit the edges
    if less_than:
        limited_edges = [edge for edge in Reader.edge_records if Reader.edge_records[edge][column] <= threshold]
    else:
        limited_edges = [edge for edge in Reader.edge_records if Reader.edge_records[edge][column] >= threshold]

    # Write to output
    Reader.writetooutput(output, only_edges=limited_edges)


def main():
    run(get_parser().parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()