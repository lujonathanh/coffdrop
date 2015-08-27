__author__ = 'jlu96'

import csv
import sys
import itertools


# Fieldnames are
# Gene pair is frozenset. Rest is header

# a write method

def get_gene_dict(pairdict):
    gene_dict = {}
    for pair in pairdict:
        gene0, gene1 = tuple(pair)
        if gene0 not in gene_dict:
            gene_dict[gene0] = set()
        gene_dict[gene0].add(gene1)
        if gene1 not in gene_dict:
            gene_dict[gene1] = set()
        gene_dict[gene1].add(gene0)
    return gene_dict



class EdgeReader:
    def __init__(self, filename, delimiter='\t', gene_column_names=frozenset(['Gene0', 'Gene1']), already_column = None,
                 include_gene_cols=False, name=None): #, already_function=None):
        self.edge_records = {}
        self.gene_column_names = gene_column_names
        self.nodes = set()
        self.repeats = set()
        self.gene_dict = {}
        if not name:
            name = filename
        self.name = name
        self.name_to_edges = {}
        self.name_to_edges[self.name] = set()

        is_triplet = len(gene_column_names) == 3
        if is_triplet:
            self.pair_dict = {}
            self.gene_to_pair_dict = {}

        print "Reading in ", filename, " ...."

        with open(filename, 'rU') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=delimiter)
            if include_gene_cols:
                self.fieldnames = [fieldname for fieldname in reader.fieldnames]
            else:
                self.fieldnames = [fieldname for fieldname in reader.fieldnames if fieldname not in gene_column_names]

            for row in reader:
                for gene_column_name in gene_column_names:
                    self.nodes.add(row[gene_column_name])

                edge = frozenset([row[gene_column_name] for gene_column_name in gene_column_names])
                if edge not in self.edge_records:
                    self.edge_records[edge] = {}
                    for fieldname in self.fieldnames:
                        try:
                            self.edge_records[edge][fieldname] = eval(row[fieldname])
                        except (NameError, SyntaxError):
                            self.edge_records[edge][fieldname] = row[fieldname]


                    genes = tuple(edge)

                    for gene in genes:
                        if gene not in self.gene_dict:
                            self.gene_dict[gene] = set()
                        self.gene_dict[gene] = self.gene_dict[gene].union(set(genes).difference(set([gene])))
                        if is_triplet:
                            # Make pair dictionary.
                            new_edge = frozenset(set(edge).difference(set([gene])))
                            if new_edge not in self.pair_dict:
                                self.pair_dict[new_edge] = set()
                            self.pair_dict[new_edge].add(gene)

                            # Make dictionary of all pairs that a gene is in triplet with.
                            if gene not in self.gene_to_pair_dict:
                                self.gene_to_pair_dict[gene] = set()
                            self.gene_to_pair_dict[gene].add(new_edge)



                    # gene1, gene2 = tuple(edge)
                    # if gene1 not in self.gene_dict:
                    #     self.gene_dict[gene1] = {gene2}
                    # else:
                    #     self.gene_dict[gene1].add(gene2)
                    #
                    # if gene2 not in self.gene_dict:
                    #     self.gene_dict[gene2] = {gene1}
                    # else:
                    #     self.gene_dict[gene2].add(gene1)

                else:
                    print "Edge already in records!"
                    self.repeats.add(edge)
                    if already_column:
                        if abs(eval(row[already_column])) > abs(self.edge_records[edge][already_column]):
                            print "replacing by larger value of ", already_column
                            for fieldname in self.fieldnames:
                                self.edge_records[edge][fieldname] = row[fieldname]
                        elif abs(eval(row[already_column])) == abs(self.edge_records[edge][already_column]):
                            print "Equal values of ", already_column, ": ", eval(row[already_column]), self.edge_records[edge][already_column]

        self.name_to_edges[self.name] = set(self.edge_records.keys())

        print 'number of repeats ', len(self.repeats)
        if already_column:
            print "Repeated edges chosen by: ", already_column
        else:
            print "Repeats were not included."


    def rename_fieldname(self, old_fieldname, new_fieldname):
        self.fieldnames.remove(old_fieldname)
        self.fieldnames.append(new_fieldname)

        for edge in self.edge_records:
            self.edge_records[edge][new_fieldname] = self.edge_records[edge][old_fieldname]
            self.edge_records[edge].pop(old_fieldname)



    def integrate(self, otherEdgeReader, new_fieldnames=[], full_integration=False):
        """
        :param otherEdgeReader: Other list of edges
        :param new_fieldnames: Column names to add to the list of edges
        :param full_integration: Append the other list of edges to this list.
        """

        if not new_fieldnames:
            new_fieldnames = [other_fieldname for other_fieldname in otherEdgeReader.fieldnames if other_fieldname not in self.fieldnames]

        self.name_to_edges[otherEdgeReader.name] = otherEdgeReader.name_to_edges[otherEdgeReader.name]



        for edge in self.edge_records:
            if edge in otherEdgeReader.edge_records:

                for fieldname in new_fieldnames:
                    self.edge_records[edge][fieldname] = otherEdgeReader.edge_records[edge][fieldname]
            else:

                for fieldname in new_fieldnames:
                    self.edge_records[edge][fieldname] = None




        for edge in otherEdgeReader.edge_records:
            if edge not in self.edge_records:

                if full_integration:
                    self.edge_records[edge] = {}
                    # Take the other edge fieldnames
                    for fieldname in new_fieldnames:
                        self.edge_records[edge][fieldname] = otherEdgeReader.edge_records[edge][fieldname]

                    # Set the current edge fieldnames to None.
                    for fieldname in self.fieldnames:
                        self.edge_records[edge][fieldname] = None


        self.fieldnames.extend(new_fieldnames)



    def get_name_to_unique_edges(self):
        """
        :return: NameToUniqueEdges: all the unique edges to each file
        """

        NameToUniqueEdges = {}

        for name in self.name_to_edges:
            name_edges = self.name_to_edges[name]

            other_names = [other_name for other_name in self.name_to_edges.keys() if other_name != name]
            other_edges = set.union(*[self.name_to_edges[other_name] for other_name in other_names])

            unique_edges = name_edges.difference(other_edges)

            NameToUniqueEdges[name] = unique_edges

        return NameToUniqueEdges


    def calc_edge_to_frequencies(self):


        self.fieldnames.append('Frequency')

        for edge in self.edge_records:
            frequency = len([name for name in self.name_to_edges if edge in self.name_to_edges[name]])
            self.edge_records[edge]['Frequency'] = frequency

    def calc_min_values(self, fieldnames, new_fieldname):
        self.fieldnames.append(new_fieldname)

        for edge in self.edge_records:
            min_fieldname = min([self.edge_records[edge][fieldname] for fieldname in fieldnames if self.edge_records[edge][fieldname] != None])
            self.edge_records[edge][new_fieldname] = min_fieldname



    def writetooutput(self, filename, delimiter='\t', NoneStandin='NA', only_edges=None, only_fieldnames=None):
        if not only_edges:
            only_edges = set(self.edge_records.keys())
        if not only_fieldnames:
            only_fieldnames = self.fieldnames


        with open(filename, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)

            writer.writerow([gene_column_name for gene_column_name in self.gene_column_names] + only_fieldnames)

            for edge in only_edges:
                writer.writerow([gene for gene in edge] +
                                [(self.edge_records[edge][fieldname] if self.edge_records[edge][fieldname]
                                    else NoneStandin) for fieldname in only_fieldnames])

    def get_pair_dict(self):
        return self.edge_records

    def get_gene_dict(self):
        return self.gene_dict

    def write_gene_dict(self, filename):
        with open(filename, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(['Gene', 'Number Of Sets'])
            for gene in self.gene_dict:
                writer.writerow([gene, len(self.gene_dict[gene])])

    def get_top_gene_lists(self, genelist1len, genelist2len):

        genelist1 = sorted(self.gene_dict.keys(), key=lambda gene: len(self.gene_dict[gene]), reverse=True)[0:genelist1len]

        genes2 = set()

        for gene1 in genelist1:
            genes2 = genes2.union(self.gene_dict[gene1])

        genelist2 = sorted(genes2, key=lambda gene2: len(self.gene_dict[gene2].intersection(set(genelist1))), reverse=True)[0:genelist2len]

        true_tested_pairs = []
        for gene1 in genelist1:
            for gene2 in genelist2:
                pair = frozenset([gene1, gene2])
                if pair in self.edge_records:
                    true_tested_pairs.append(pair)

        print "Fraction of true pairs from testing gene lists of length", len(genelist1), "and", len (genelist2), ": ", len(true_tested_pairs) * 1.0 / (len(genelist1) * len(genelist2))

        return genelist1, genelist2, true_tested_pairs


    def get_top_gene_lists_triplets(self, genelist1len, genelist2len, genelist3len):
        genelist1 = sorted(self.gene_dict.keys(), key=lambda gene: len(self.gene_dict[gene]), reverse=True)[0:genelist1len]

        genes2 = set()

        for gene1 in genelist1:
            genes2 = genes2.union(self.gene_dict[gene1])

        genelist2 = sorted(genes2, key=lambda gene2: len(self.gene_dict[gene2].intersection(set(genelist1))), reverse=True)[0:genelist2len]

        pairs = set()
        genes3 = set()
        for gene1 in genelist1:
            for gene2 in genelist2:
                pair = frozenset([gene1, gene2])
                if pair in self.pair_dict:
                    pairs.add(pair)
                    genes3 = genes3.union(self.pair_dict[pair])

        genelist3 = sorted(genes3, key=lambda gene3: len(self.gene_to_pair_dict[gene3].intersection(pairs)), reverse=True)[0:genelist3len]


        return genelist1, genelist2, genelist3



    # Write the top gene lists of triplets.
    def write_top_gene_lists_triplets(self, truetripletsfile, genelist1len, genelist2len, genelist3len,
                                      genelistfile=None):

        genelist1, genelist2, genelist3 = self.get_top_gene_lists_triplets(genelist1len, genelist2len, genelist3len)

        num_true_triplets = 0
        with open(truetripletsfile, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(['Gene1', 'Gene2', 'Gene3', 'TrueTriplet'])
            for gene1 in genelist1:
                for gene2 in genelist2:
                    for gene3 in genelist3:
                        is_true_triplet = 1 if frozenset([gene1, gene2, gene3]) in self.edge_records else 0
                        writer.writerow([gene1, gene2, gene3, is_true_triplet])
                        num_true_triplets += is_true_triplet

        print "Number of True Triplets ", num_true_triplets
        print "Percentage of True Triplets ", num_true_triplets * 1.0 / (len(genelist1) * len(genelist2) * len(genelist3))



        if not genelistfile:
            genelistfile = truetripletsfile + '_genelists'

        with open(genelistfile, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(['Gene1', 'Gene2', 'Gene3'])
            for i in range(max([len(genelist1), len(genelist2), len(genelist3)])):
                gene1 = genelist1[i] if i < len(genelist1) else ""
                gene2 = genelist2[i] if i < len(genelist2) else ""
                gene3 = genelist3[i] if i < len(genelist3) else ""

                writer.writerow([gene1, gene2, gene3])





    # Write the top gene lists.
    def write_top_gene_lists(self, truepairsfile, genelist1len, genelist2len, genelistfile=None):
        genelist1, genelist2, true_tested_pairs = self.get_top_gene_lists(genelist1len, genelist2len)


        with open(truepairsfile, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(['Gene1', 'Gene2', 'TruePair'])
            for gene1 in genelist1:
                for gene2 in genelist2:
                    is_true_pair = 1 if frozenset([gene1, gene2]) in true_tested_pairs else 0
                    writer.writerow([gene1, gene2, is_true_pair])

        if not genelistfile:
            genelistfile = truepairsfile + '_genelists'

        with open(genelistfile, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(['Gene1', 'Gene2'])
            for i in range(max([len(genelist1), len(genelist2)])):
                gene1 = genelist1[i] if i < len(genelist1) else ""
                gene2 = genelist2[i] if i < len(genelist2) else ""

                writer.writerow([gene1, gene2])




def get_parser():
    # Parse arguments
    import argparse
    description = 'Convert several edge file to one edge file.'
    parser = argparse.ArgumentParser(description = description)

    # General parameters
    parser.add_argument('-i', '--input_edges', nargs='+', required=True, help="List of files to integrate."
                                                                              "Add the columns/fields of all the"
                                                                              "files to the first file."
                                                                              "If there are repeated columns, they"
                                                                              "won't be overwritten.")

    parser.add_argument('-o', '--output_name', default=None)

    parser.add_argument('-nf', '--new_fieldnames', nargs='+', default=[], help="List of columns to add from one to the other."
                                                                               "If not specified, add all columns.")

    parser.add_argument('-gcn', '--gene_column_names', nargs='+', default=frozenset(['Gene0', 'Gene1']))

    parser.add_argument('-gdf', '--gene_dict_file', default=None, help='Output name to write each gene and number of sets'
                                                                       'it is in')

    parser.add_argument('-bgl', '--best_gene_list_file', default=None, help='Name of output file. Edgereader will generate'
                                                                            'a list of template and query genes for testing'
                                                                            'pairs/triplets. It will write all the sets that will'
                                                                            'be tested from those lists and whether they were'
                                                                            'true sets or not.')

    parser.add_argument('-gll', '--gene_list_lens', nargs='+', type=int, default=None, help="Lengths of template and query"
                                                                                            "gene lists. In descending order"
                                                                                            "of lengths.")

    parser.add_argument('-fi', '--full_integration', default=False, help="When integrating the pairs,"
                                                                         "include the pairs that are not"
                                                                         "in the original file as well.")

    return parser

def run(args):
    inputs = args.input_edges
    output = args.output_name
    new_fieldnames = args.new_fieldnames
    gene_column_names = args.gene_column_names
    gene_dict_file = args.gene_dict_file
    best_gene_lists_file = args.best_gene_list_file
    gene_list_lens = args.gene_list_lens
    full_integration = args.full_integration

    #already_column = args.already_column


    completeEdgeReader = EdgeReader(inputs[0], gene_column_names=gene_column_names)
    for i in range(1,len(inputs)):
        newEdgeReader = EdgeReader(inputs[i], gene_column_names=gene_column_names)
        completeEdgeReader.integrate(newEdgeReader, new_fieldnames=new_fieldnames, full_integration=full_integration)


    if output:
        print "Writing to output", output
        completeEdgeReader.writetooutput(output)

    if gene_dict_file:
        completeEdgeReader.write_gene_dict(gene_dict_file)

    if best_gene_lists_file:
        if len(gene_column_names) == 2:
            completeEdgeReader.write_top_gene_lists(best_gene_lists_file, *gene_list_lens)
        elif len(gene_column_names) == 3:
            completeEdgeReader.write_top_gene_lists_triplets(best_gene_lists_file, *gene_list_lens)


def main():
    run(get_parser().parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()