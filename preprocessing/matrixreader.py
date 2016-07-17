__author__ = 'jlu96'

import csv
import time
import os
import sys
import numpy as np


class MatrixReader:
    def __init__(self, *filenames):

        self.sample_records = {}
        self.mutation_list = set()

        for filename in filenames:
            with open(filename) as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')

                for row in reader:
                    sample = row[0]
                    if sample not in self.sample_records:
                        self.sample_records[sample] = set(row[1:])
                        self.mutation_list = self.mutation_list.union(set(row[1:]))
                    else:
                        self.sample_records[sample] = self.sample_records[sample].union(set(row[1:]))
                        self.mutation_list = self.mutation_list.union(set(row[1:]))


    def getmutationlist(self):
        return self.mutation_list

    def getnumberalterations(self):
        num = 0
        for sample in self.sample_records:
            num += len(self.sample_records[sample])
        return num

    def writemutationmatrix(self, filename):
        with open(filename, 'w') as opener:
            for sample in self.sample_records:
                opener.write('\t'.join([sample] + [record for record in self.sample_records[sample]]) + '\n')

    def filter(self, geneFile=None):


        if geneFile:
            with open(geneFile) as f:
                genes = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
        else:
            genes = set()


        newMatrix = MatrixReader()
        newMatrix.sample_records = self.sample_records.copy()
        newMatrix.mutation_list = genes.intersection(self.mutation_list)


        for sample in self.sample_records:
            newMatrix.sample_records[sample] = genes.intersection(self.sample_records[sample])
            if not newMatrix.sample_records[sample]:
                newMatrix.sample_records.pop(sample)


        return newMatrix

    def filter_samples(self, sampleFile=None):


        if sampleFile:
            with open(sampleFile) as f:
                samples = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
        else:
            samples = set()


        newMatrix = MatrixReader()
        newMatrix.sample_records = {}
        for sample in samples:
            newMatrix.sample_records[sample] = self.sample_records[sample]
        return newMatrix


    def permute_by_gene(self, random_seed=123):
        gene_records = self.get_gene_records()
        geneToFreq = {}

        for gene in gene_records:
            geneToFreq[gene] = len(gene_records[gene])

        samples = self.sample_records.keys()

        new_gene_records = {}

        for gene in gene_records:
            new_gene_records[gene] = set(np.random.choice(samples, geneToFreq[gene], replace=False))

        new_sample_records = convert_sample_gene(new_gene_records)



        newMatrix = MatrixReader()
        newMatrix.sample_records = new_sample_records
        newMatrix.mutation_list = self.mutation_list

        return newMatrix


    def removeEndingIn(self, ending, length):
        genes = set([g for g in self.mutation_list if g[-1 * length:] != ending])

        newMatrix = MatrixReader()
        newMatrix.sample_records = self.sample_records.copy()
        newMatrix.mutation_list = genes.intersection(self.mutation_list)

        for sample in self.sample_records:
            newMatrix.sample_records[sample] = genes.intersection(self.sample_records[sample])
            if not newMatrix.sample_records[sample]:
                newMatrix.sample_records.pop(sample)


        return newMatrix


    def get_gene_records(self):
        gene_records = {}
        for sample in self.sample_records:
            for gene in self.sample_records[sample]:
                if gene not in gene_records:
                    gene_records[gene] = set()
                gene_records[gene].add(sample)

        return gene_records

    def write_weSME_mut_list(self, filename):
        sampleToNumber = {}
        samples = self.sample_records.keys()
        for i in range(len(self.sample_records)):
            sample = samples[i]
            sampleToNumber[sample] = i
        gene_records = self.get_gene_records()
        geneToSampleNumbers = {}
        for gene in gene_records:
            geneToSampleNumbers[gene] = set([sampleToNumber[g] for g in gene_records[gene]])

        with open(filename, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(["samples", ",".join(samples)])
            for gene in geneToSampleNumbers:
                writer.writerow([gene, ",".join([str(num) for num in sorted(geneToSampleNumbers[gene])])])


def convert_sample_gene(sample_records):
    gene_records = {}
    for sample in sample_records:
        for gene in sample_records[sample]:
            if gene not in gene_records:
                gene_records[gene] = set()
            gene_records[gene].add(sample)

    return gene_records




def get_parser():
    # Parse arguments
    import argparse
    description = 'Combine several matrix files to one matrix file.'
    parser = argparse.ArgumentParser(description = description)


    parser.add_argument('-i', '--input_matrices', nargs='+', required=True)

    parser.add_argument('-o', '--output_name', required=True)

    return parser

def run(args):
    inputs = args.input_matrices
    output = args.output_name
    t1 = time.time()
    print "file sizes: ", [os.path.getsize(input) for input in inputs]

    matrixreader = MatrixReader(*inputs)
    t2 = time.time()
    print "time used to read in: ", t2 - t1

    matrixreader.writemutationmatrix(filename=output)
    t3 = time.time()

    print "time used to write matrix: ", t3 - t2
    print "Combined mutation matrix written to ", output

def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()