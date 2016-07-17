__author__ = 'jlu96'

import pandas as pd
import sys

def get_parser():
    # Parse arguments
    import argparse
    description = 'Convert several edge file to one edge file.'
    parser = argparse.ArgumentParser(description = description)

    # General parameters
    parser.add_argument('-i', '--input_edges', nargs='+', required=True, help="List of files to integrate."
                                                                      "won't be overwritten.")

    parser.add_argument('-o', '--output_name', default=None)

    parser.add_argument('-gcn', '--gene_column_names', nargs='+', default=frozenset(['Gene0', 'Gene1']))

    return parser

def run(args):
    inputs = args.input_edges
    output = args.output_name
    gene_column_names = args.gene_column_names


# read in files
    dfs = []

    for input_file in inputs:
        dfs.append(pd.read_csv(input_file, sep='\t'))


    # make column combining gene names by alphabetic order
    for df in dfs:
        gene_columns = [df[gene_column_name] for gene_column_name in gene_column_names]
        zipped = zip(*gene_columns)

        sorted_zipped = [sorted(zip_pair) for zip_pair in zipped]
        combined_zipped = ["-".join(s) for s in sorted_zipped]
        df['Genes'] = combined_zipped


    genes_set = set(dfs[0]['Genes'].values)
    for df in dfs:
        genes_set = genes_set.intersection(set(df['Genes'].values))

    print "Number of overlapps among files", len(genes_set)

    new_dfs = []
    for df in dfs:

        new_df = df[df['Genes'].isin(genes_set)]

        new_df = new_df.sort_values(by='Genes')
        new_df.reset_index(inplace=True)
        new_dfs.append(new_df)

    out_df = pd.concat(new_dfs, axis=1)
    out_df = out_df.drop('index', axis=1)

    out_df.to_csv(output, sep='\t')

    print "Combined data frames file written to ", output

def main():
    run(get_parser().parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()