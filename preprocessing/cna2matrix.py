__author__ = 'jlu96'

import csv
import sys
import time
import os
import collections
import numpy as np


class Segment:
    # Each geneDict is dict of 0s to patients, 1s, 2s, -1s, -2s, as well as position and chromosome
    def __init__(self, gene, altDict):
        self.alterations = [-2, -1, 0, 1, 2]

        self.dict = altDict.copy()

        # print "in seg dict is ", self.dict
        # print "in seg gene is ", gene

        self.genes = {}
        self.genes[gene] = altDict.copy()
        self.name = gene

        self.patientToGenes = {}
        # An alteration is a 1 or 0 or -1, etc.
        for alteration in self.alterations:
            if alteration in altDict:
                for patient in altDict[alteration]:
                    if patient not in self.patientToGenes:
                        self.patientToGenes[patient] = {}
                    if alteration not in self.patientToGenes[patient]:
                        self.patientToGenes[patient][alteration] = set()
                    self.patientToGenes[patient][alteration].add(gene)

        # print "initial segment: ", self.name
        # print "intial gene: ", gene
        # print "initial genes: ", self.genes
        # print "initial patientToGenes ", self.patientToGenes

    def integrate(self, gene, altDict):

        #print "new gene is "
        for alteration in self.alterations:
            if alteration in altDict:
                for patient in altDict[alteration]:
                    if patient not in self.patientToGenes:
                        self.patientToGenes[patient] = {}
                    if alteration not in self.patientToGenes[patient]:
                        self.patientToGenes[patient][alteration] = set()
                    self.patientToGenes[patient][alteration].add(gene)

        # If you do integrate, make sure you only test the concordance between
        # the actual affected genes.
        # for patient in self.patientToGenes:

        self.dict.update(altDict)
        self.genes[gene] = altDict.copy()
        self.name = self.name + '_' + gene

        # update the list based off the patientToGenes

    def calc_representation(self):
        # For each patient
        # Find ratio of what we say about that patient per alteration, over total what
        # each gene says about that patient per alteration
        # Average these, print standard deviation.
        # Find the most representative.

        # print "For segment ", self.name, " genes are ", [gene for gene in self.genes]

        for alteration in self.dict:
            self.dict[alteration] = set()

        self.num_zeros = 0
        self.ratios = []

        for patient in self.patientToGenes:
            # print "Patient to genes is ", self.patientToGenes[patient]
            # print "Keys are ", self.patientToGenes[patient].keys()
            alterations = sorted(self.patientToGenes[patient].keys(), key=lambda alteration: len(self.patientToGenes[patient][alteration]),
                                 reverse=True)
            self.dict[alterations[0]].add(patient)

            num_each_alteration = [len(self.patientToGenes[patient][alteration]) for alteration in alterations]
            ratio = num_each_alteration[0] * 1.0 / sum(num_each_alteration)
            if alterations[0] == 0:
                self.num_zeros += 1
            self.ratios.append(ratio)

        self.rep = sum(self.ratios)/len(self.ratios)


        for gene in self.genes:
            self.genes[gene]['Gene'] = gene
            self.genes[gene]['Segment'] = self.name
            self.genes[gene]['Concordance'] = concordance_factor(self.genes[gene], self.dict)
            self.genes[gene]['NoZeroConcordance'] = concordance_factor(self.genes[gene], self.dict, dictvalues=[1,2,-1,-2])

    def get_writable_dict_list(self):
        # self.calc_representation()

        return self.genes.values()


def concordance_factor(a_dict, b_dict, dictvalues=[0, 1, 2, -1, -2]):
    intersect_set = set()
    union_set = set()
    for dictvalue in dictvalues:
        try:
            intersect_set = intersect_set.union(set.intersection(a_dict[dictvalue], b_dict[dictvalue]))
            union_set = union_set.union(a_dict[dictvalue], b_dict[dictvalue])
        except KeyError:
            pass

    return len(intersect_set) * 1.0 /len(union_set)



def chrom_compare(a, b, chromosomes=['X', 'Y'] + [str(i) for i in range(1, 23)]):
    if a['Chromosome'] not in chromosomes:
        return 1
    if b['Chromosome'] not in chromosomes:
        return -1

    if a['Chromosome'] != b['Chromosome']:
        return chromosomes.index(a['Chromosome']) - chromosomes.index(b['Chromosome'])
    else:
        return a['Start'] - b['Start']


class CNA:
    def __init__(self, filename, genecol = 'Gene Symbol', notsamplecols = set(['Cytoband', 'Locus ID', 'Locus.ID', 'gene_name', 'GENE_ID', 'COMMON', ''])):
        self.records = {}

        with open(filename) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            header = reader.fieldnames
            for sample in header:
                if sample not in notsamplecols:
                # if sample != 'gene_name' and sample != 'Locus.ID' and sample != 'Cytoband':
                    self.records[sample] = {}

            for row in reader:
                gene = row[genecol]
                #print row
                for sample in self.records:
                    try:
                        self.records[sample][gene] = eval(row[sample])
                    except (NameError, SyntaxError):
                        self.records[sample][gene] = row[sample]
                    # if row[sample] == 'NaN':
                    #     self.records[sample][gene] = 0
                    # else:
                    #     self.records[sample][gene] = eval(row[sample])

    def writemutationmatrix(self, filename, gain_suffix = 'gain', loss_suffix = 'loss', level = None):

        num_alterations = 0
        with open(filename, 'w') as opener:
            for sample in self.records:
                amplified_genes = [gene + gain_suffix for gene in self.records[sample] if self.records[sample][gene] > 0]
                lost_genes = [gene + loss_suffix for gene in self.records[sample] if self.records[sample][gene] < 0]

                num_alterations += len(amplified_genes) + len(lost_genes)
                if level:
                    opener.write('\t'.join([sample[:level]] + amplified_genes + lost_genes) + '\n')
                else:
                    opener.write('\t'.join([sample] + amplified_genes + lost_genes) + '\n')

        print "Total number of gene alterations: ", num_alterations

    def convert_to_genes(self):
        self.geneDict = {}

        # Each gene has: chromosome, position, 2 position, 1 patients, etc.

        for sample in self.records:
            for gene in self.records[sample]:
                if gene not in self.geneDict:
                    self.geneDict[gene] = {}
                if self.records[sample][gene] not in self.geneDict[gene]:

                    self.geneDict[gene][self.records[sample][gene]] = set()
                self.geneDict[gene][self.records[sample][gene]].add(sample)

    def load_gene_positions(self, filename='gene_positions.txt', genecol='Associated Gene Name', chromcol='Chromosome Name',
                       startcol='Gene Start (bp)', endcol='Gene End (bp)', delimiter='\t',
                       chromosomes = set(['X', 'Y'] + [str(i) for i in range(1, 23)])):
        t = time.time()
        geneToStat = {}
        with open(filename,'rU') as geneToStatFile:
            statreader = csv.DictReader(geneToStatFile, delimiter=delimiter)
            for row in statreader:
                if row[chromcol] in chromosomes:
                    gene = row[genecol]
                    position = {}
                    position['Chromosome'] = row[chromcol]
                    position['Start'] = eval(row[startcol])
                    position['End'] = eval(row[endcol])

                    if gene not in geneToStat:
                        geneToStat[gene] = position

        self.loaded_genes = set()
        self.not_loaded_genes = set()

        for gene in self.geneDict:
            if gene in geneToStat:
                self.geneDict[gene].update(geneToStat[gene])
                self.loaded_genes.add(gene)
            else:
                #self.geneDict.pop(gene)
                self.geneDict[gene]['Chromosome'] = gene
                self.geneDict[gene]['Start'] = 0
                self.geneDict[gene]['End'] = 0
                self.not_loaded_genes.add(gene)


        print len(self.loaded_genes), "genes were loaded"
        print len(self.not_loaded_genes), "genes were not loaded"


                    # if gene in self.geneDict:
                    #     self.geneDict[gene]['Chromosome'] = row[chromcol]
                    #     self.geneDict[gene]['Start'] = eval(row[startcol])
                    #     self.geneDict[gene]['End'] = eval(row[endcol])
                        # print "Position for gene", gene, "found as ", position, " and ", geneToStat[gene]



        # print "Time used to load gene positions ", time.time() - t
    def bin_genes(self, concordance_thresh=0.96, distance_thresh=10000000):

        print "Binning genes."
        print "Adjacent concordance threshold: ", concordance_thresh
        print "Distance threshold: ", distance_thresh


        # Sort the loaded_genes by chromosome, start, and end.
        sorted_genes = sorted(self.loaded_genes, key=lambda entry: self.geneDict[entry], cmp=chrom_compare)

        # print [(self.geneDict[g]['Chromosome'], self.geneDict[g]['Start']) for g in sorted_genes]

        self.segments = [gene for gene in sorted_genes]
        self.geneToSegments = {}


        num_segments = 0
        self.segments[num_segments] = Segment(sorted_genes[0], self.geneDict[sorted_genes[0]])
        prevChrom = self.segments[num_segments].dict['Chromosome']
        sclass = self.segments[num_segments].__class__

        iter_genes = iter(sorted_genes)
        iter_genes.next()


        for gene in iter_genes:
            #print "Gene was ", gene
            chrom = self.geneDict[gene]["Chromosome"]
            if chrom != prevChrom or \
                abs(self.segments[num_segments].dict['Start'] - self.geneDict[gene]['Start']) > distance_thresh \
                or concordance_factor(self.segments[num_segments].dict, self.geneDict[gene]) < concordance_thresh:
                #print "Gene is now ", gene
                self.segments[num_segments].calc_representation()
                #print "Now gene is ", gene

                for past_gene in self.segments[num_segments].genes:
                    self.geneToSegments[past_gene] = self.segments[num_segments].name

                num_segments += 1
                # print prevChrom
                # print "Before: ", self.segments[num_segments]
                # print "gene is ", gene
                # print "self.geneDict[gene] is ", self.geneDict[gene]
                self.segments[num_segments] = sclass(gene, self.geneDict[gene])
                # print "After: ", self.segments[num_segments]
                # print "After, gene is ", gene

                # print "Now gene is ", gene
                prevChrom = chrom


            else:
                self.segments[num_segments].integrate(gene, self.geneDict[gene])
                # print "Gene ", gene, " integrated"

        # For the last segment.
        self.segments[num_segments].calc_representation()
        for past_gene in self.segments[num_segments].genes:
            self.geneToSegments[past_gene] = self.segments[num_segments].name
        num_segments += 1

        print "Number of position segments ", num_segments

        # Now load the unloaded genes


        # self.not_added_genes = []


        self.segments = self.segments[:num_segments]


        print "Number of segments ", len(self.segments)
        print "Number of genes ", len(self.geneDict)




    # def write_segments_GISTIC(self, filename):

    def write_segments_matrix(self, filename, gain_suffix = 'gain', loss_suffix = 'loss', level = None):

        majority_ratios = []

        num_equal = 0
        num_alterations = 0

        with open(filename, 'w') as opener:
            for sample in self.records:
                # Number of each of the genes...
                gain_segs = [self.geneToSegments[gene] for gene in self.records[sample] if self.records[sample][gene] > 0
                             and gene in self.loaded_genes]
                loss_segs = [self.geneToSegments[gene] for gene in self.records[sample] if self.records[sample][gene] < 0
                             and gene in self.loaded_genes]

                gain_seg_set = set(gain_segs)
                loss_seg_set = set(loss_segs)

                if gain_seg_set.intersection(loss_seg_set):
                    for seg in gain_seg_set.intersection(loss_seg_set):
                        # print "Segment ", seg, "is annotated as gain and loss for sample ", sample
                        if gain_segs.count(seg) > loss_segs.count(seg):
                            majority_ratios.append(gain_segs.count(seg) * 1.0 / (gain_segs.count(seg) + loss_segs.count(seg)))
                            loss_seg_set.remove(seg)
                        elif gain_segs.count(seg) < loss_segs.count(seg):
                            majority_ratios.append(loss_segs.count(seg) * 1.0 / (gain_segs.count(seg) + loss_segs.count(seg)))
                            gain_seg_set.remove(seg)
                        else:
                            # print "Equal # gains and losses for seg: ", seg, " and sample ", sample
                            # print "Writing as gains. "
                            loss_seg_set.remove(seg)
                            num_equal += 1
                            # majority_ratios.append(0.5)


                gain_seg_set = set([seg + gain_suffix for seg in gain_seg_set])
                loss_seg_set = set([seg + loss_suffix for seg in loss_seg_set])

                num_alterations += (len(gain_seg_set) + len(loss_seg_set))

                if level:
                    opener.write('\t'.join([sample[:level]] + list(gain_seg_set) + list(loss_seg_set)) + '\n')
                else:
                    opener.write('\t'.join([sample] + list(gain_seg_set) + list(loss_seg_set)) + '\n')

        print "****************"
        print "WRITING SEGMENT MATRIX:"
        print "-----------------------"
        print "Total number of segment alterations: ", num_alterations
        print num_equal, "segments had equal numbers of gain and losses for  a sample were found, written as gain alterations"
        print len(majority_ratios), " had some gains and some losses, but were overridden by majority"
        print "Average majority ratio: ", sum(majority_ratios)/len(majority_ratios)
        print "Standard Dev of majority ratios: ", np.std(majority_ratios)
        print "*****************"


    # Finish writing this dictionary. Not yet.
    def write_gene_to_segments(self, filename):


        concordances = []
        no_zero_concordances = []

        fieldnames = ['Gene', 'Segment', 'Concordance', 'NoZeroConcordance']
        with open(filename, 'w') as opener:
            writer = csv.DictWriter(opener, delimiter='\t', extrasaction='ignore', fieldnames=fieldnames)
            writer.writeheader()
            for segment in self.segments:
                dictList = segment.get_writable_dict_list()
                for dict in dictList:
                    writer.writerow(dict)
                    concordances.append(dict['Concordance'])
                    no_zero_concordances.append(dict['NoZeroConcordance'])

        print "********************"
        print "WRITING SEG2GENES:"
        print "--------------------"
        print "Number of genes ", len(self.geneToSegments)
        print "Number of segments ", len(self.segments)
        print "Average concordance: ", sum(concordances)/len(concordances)
        print "Standard deviation of concordance: ", np.std(concordances)
        print "Average nonzero concordances: ", sum(no_zero_concordances)/ len(no_zero_concordances)
        print "Standard deviation of nonzero concordances: ", np.std(no_zero_concordances)




def get_parser():
    # Parse arguments
    import argparse
    description = 'Convert one cna file to one matrix file.'
    parser = argparse.ArgumentParser(description = description)

    # General parameters
    parser.add_argument('-i', '--input_cna', required=True, help='Input GISTIC gene-thresholded file.')

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-s', '--segment_output_name', default=None, help='Name of output segment matrix.')

    parser.add_argument('-gl', '--gain_loss', type=int, default=1, help='Set to 1 to append "gain" to genes that are'
                                                                        'gained and "loss" to genes that are lost.')

    parser.add_argument('-l', '--level', type=int, default=12, help='Last index of tumor barcode'
                                                                 'to limit to')

    parser.add_argument('-ct', '--concordance_thresh', type=float, default=0.95, help='The concordance threshold '
                                                                                     'to use for each gene as you '
                                                                                     'go along through the segment'
                                                                                     '(not the actual concordance threshold'
                                                                                     'for each gene')
    parser.add_argument('-dt', '--distance_thresh', type=float, default=10000000, help='Distance threshold between genes')

    return parser

def run(args):
    input = args.input_cna
    output = args.output_name
    gain_loss = args.gain_loss
    if args.segment_output_name:
        segmatrix = args.segment_output_name
        gene2seg = args.segment_output_name + '.gene2seg'
    concordance_thresh = args.concordance_thresh
    distance_thresh = args.distance_thresh

    level = args.level
    t1 = time.time()
    print "cna file size: ", os.path.getsize(input)
    print "Reading in CNA file. This may take a few minutes..."

    cnareader = CNA(input)
    t2 = time.time()
    print "time used to read in: ", t2 - t1




    if not args.segment_output_name:
        if gain_loss == 1:
            gain_suffix = 'gain'
            loss_suffix = 'loss'
        else:
            gain_suffix = ''
            loss_suffix = ''

        print "gain suffix: ", gain_suffix
        print "loss suffix: ", loss_suffix

        if level == 0:
            level = None

        cnareader.writemutationmatrix(filename=output, gain_suffix=gain_suffix, loss_suffix=loss_suffix, level=level)
        t3 = time.time()

        print "time used to write matrix: ", t3 - t2
        print "Mutation matrix written to ", output


        if gain_loss == 2:
            gain_suffix = 'gain'
            loss_suffix = 'loss'
            print "gain suffix: ", gain_suffix
            print "loss suffix: ", loss_suffix

            output = args.output_name + '-gl' + '.m2'
            cnareader.writemutationmatrix(filename=output, gain_suffix=gain_suffix, loss_suffix=loss_suffix, level=level)
            t4 = time.time()
            print "time used to write gain-loss matrix: ", t4 - t3

    if args.segment_output_name:
        t2 = time.time()
        cnareader.convert_to_genes()
        t3 = time.time()
        print "Time used to convert genes: ", t3 - t2
        cnareader.load_gene_positions()
        t4 = time.time()
        print "Time used to get gene positions: ", t4 - t3
        cnareader.bin_genes(concordance_thresh=concordance_thresh, distance_thresh=distance_thresh)
        t5 = time.time()
        print "Time used to bin genes: ", t5 - t4
        cnareader.write_segments_matrix(filename=segmatrix, level=level)
        cnareader.write_gene_to_segments(filename=gene2seg)
        print "Segments matrix written to: ", segmatrix
        print "Gene to segments written to: ", gene2seg

def main():
    run(get_parser().parse_args(sys.argv[1:]))



if __name__ == '__main__':
    main()