__author__ = 'jlu96'

import edgereader as edg
import sys


def get_parser():
    # Parse arguments
    import argparse
    description = 'Annotate several edge files and get unique files.'
    parser = argparse.ArgumentParser(description = description)

    # General parameters
    parser.add_argument('-i', '--edgefiles', nargs='+', required=True, help="List of files to integrate."
                                                                              "Add the columns/fields of all the"
                                                                              "files to the first file."
                                                                              "If there are repeated columns, they"
                                                                              "won't be overwritten.")
    parser.add_argument('-n', '--names', nargs='+',default=None, help="Shorter names than each filename")

    parser.add_argument('-pn', '--pvalnames', nargs='+', default=None, help="Pvalue column names to calculate "
                                                                            "minimum probabilities over.")

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-gcn', '--gene_column_names', nargs='+', default=frozenset(['Gene0', 'Gene1']))




    return parser

def run(args):
    # Input: List of edge files and pvalue columns, and names to assign (default is same as list of pair files)

    edgefiles = args.edgefiles
    names = args.names
    pvalnames = args.pvalnames
    if not names:
        names = edgefiles
    if not pvalnames:
        pvalnames = [name + 'Probability' for name in names]
    gene_column_names = args.gene_column_names

    output = args.output_name



    # Integrate the edge files, and make each probability field name specific to that file.
    completeEdgeReader = edg.EdgeReader(edgefiles[0], gene_column_names=gene_column_names, name=names[0])
    if pvalnames[0] not in completeEdgeReader.fieldnames:
        completeEdgeReader.rename_fieldname('Probability', pvalnames[0])

    for i in range(1,len(edgefiles)):
        newEdgeReader = edg.EdgeReader(edgefiles[i], gene_column_names=gene_column_names,  name=names[i])

        if pvalnames[i] not in newEdgeReader.fieldnames:
            newEdgeReader.rename_fieldname('Probability', pvalnames[i])

        completeEdgeReader.integrate(newEdgeReader, new_fieldnames=pvalnames[i:i+1], full_integration=True)




    # Add new frequency column and best probability (using pvalnames), and probability diff b/n two columns
    completeEdgeReader.calc_edge_to_frequencies()
    completeEdgeReader.calc_min_values(fieldnames=pvalnames, new_fieldname="MinProbability")

    print "Integrated files and unique edges to each file will be written with prefix ", output
    completeEdgeReader.writetooutput(output)


     # Write out the unique edges to each file
    NameToUniqueEdges = completeEdgeReader.get_name_to_unique_edges()
    print "Shared edges: ", len(completeEdgeReader.edge_records) - sum([len(NameToUniqueEdges[name]) for name in NameToUniqueEdges])
    for i in range(len(names)):
        name = names[i]
        print "Unique edges to ", name, ": ", len(NameToUniqueEdges[name])
        completeEdgeReader.writetooutput(output+ '_' + name + '_uniqueedges.tsv', only_edges=NameToUniqueEdges[name],
                                         only_fieldnames=[pvalnames[i]])


def main():
    run(get_parser().parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()




