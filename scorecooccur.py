__author__ = 'jlu96'

import csv
import operator
from scipy.misc import comb
from scipy.misc import factorial
import math

def load_gene_lengths(filename):
    geneToLength = {}
    with open(filename, 'rU') as geneToLengthFile:
        lengthreader = csv.reader(geneToLengthFile, delimiter='\t')
        for row in lengthreader:
            geneToLength[row[0]] = eval(row[1])
    return geneToLength

# def score_cooccur_set(genes, geneToCases, patientToGenes):
#     # Find all the patients that have all the genes cooccurring
#     patients = set.intersection(*[set(geneToCases[gene]) for gene in genes])
#
#     # print patients
#     # print [len(patientToGenes[patient]) for patient in patients]
#
#     patient_scores = set()
#
#     # For each patient,
#     for patient in patients:
#         # Get the number of mutations.
#         patient_mutation_num = len(patientToGenes[patient])
#
#         # Calculate a proxy for the probability of cooccurrence: Product of  the k gene lengths * nPk
#         patient_score = patient_cooccur_score(patient_mutation_num, genes)
#         if not patient_score:
#             print patient_mutation_num
#             print genes
#             print [gene_to_length(gene) for gene in genes]
#         patient_scores.add(patient_score)
#
#     # The total score is the logarithm of the reciprocal of the sum of the reciprocals of the individual scores.
#     # Think the logarithm of parallel circuits.
#
#     #print sum([1.0/score for score in patient_scores])
#     # [1.0/score for score in patient_scores]
#     total_score = math.log(1.0 / sum([1.0/score for score in patient_scores]))
#
#     return total_score



def score_cooccur_set(genes, geneToCases, patientToGenes):
    if len(genes) == 3: return score_cooccur_set_three(genes, geneToCases, patientToGenes)
    # Find all the patients that have all the genes cooccurring
    patients = set.intersection(*[set(geneToCases[gene]) for gene in genes])

    gene0, gene1 = tuple(genes)
    l0, l1 = gene_to_length(gene0), gene_to_length(gene1)
    # print patients
    # print [len(patientToGenes[patient]) for patient in patients]

    patient_scores = set()

    # For each patient,
    for patient in patients:
        # Get the number of mutations.
        patient_mutation_num = len(patientToGenes[patient])

        # Calculate a proxy for the probability of cooccurrence: Product of  the k gene lengths * nPk

        patient_score = patient_mutation_num * (patient_mutation_num - 1) * l0 * l1
        patient_scores.add(patient_score)


    # The total score is the logarithm of the reciprocal of the sum of the reciprocals of the individual scores.
    # Think the logarithm of parallel circuits.

    total_score = math.log(1.0 / sum([1.0/score for score in patient_scores]))
    return total_score

def score_cooccur_set_three(genes, geneToCases, patientToGenes):
    # Find all the patients that have all the genes cooccurring
    patients = set.intersection(*[set(geneToCases[gene]) for gene in genes])

    gene0, gene1, gene2 = tuple(genes)
    l0, l1, l2 = gene_to_length(gene0), gene_to_length(gene1), gene_to_length(gene2)
    # print patients
    # print [len(patientToGenes[patient]) for patient in patients]

    patient_scores = set()

    # For each patient,
    for patient in patients:
        # Get the number of mutations.
        patient_mutation_num = len(patientToGenes[patient])

        # Calculate a proxy for the probability of cooccurrence: Product of  the k gene lengths * nPk

        patient_score = patient_mutation_num * (patient_mutation_num - 1) * (patient_mutation_num - 2) * l0 * l1 * l2
        patient_scores.add(patient_score)


    # The total score is the logarithm of the reciprocal of the sum of the reciprocals of the individual scores.
    # Think the logarithm of parallel circuits.

    total_score = math.log(1.0 / sum([1.0/score for score in patient_scores]))
    return total_score




def patient_cooccur_score(patient_mutation_num, genes):
    # try:
    #     geneToLength = patient_cooccur_score.geneToLength
    # except AttributeError:
    #     patient_cooccur_score.geneToLength = load_gene_lengths(filename)
    #     geneToLength = patient_cooccur_score.geneToLength

    lengths = [gene_to_length(gene) for gene in genes]
    k = len(genes)

    # Calculate a proxy for the probability of cooccurrence: Product of  the k gene lengths * nPk
    return prod(lengths) * comb(patient_mutation_num, k) * factorial(k)


# Average gene length is taken from: wikipedia of coding sequence length
# https://en.wikipedia.org/wiki/Human_genome#Coding_sequences_.28protein-coding_genes.29
# 20200 proteins,
# "Over the whole genome, the median size of an exon is 122 bp (mean = 145 bp),
# the median number of exons is 7 (mean = 8.8), and the median coding sequence encodes 367 amino acids (mean = 447 amino acids)
# Use median = 7 * 122 = 854

def gene_in_file(gene, filename='geneToLength_all_firstseen.txt'):
    try:
        geneToLength = gene_to_length.geneToLength
    except AttributeError:
        gene_to_length.geneToLength = load_gene_lengths(filename)
        geneToLength = gene_to_length.geneToLength

    return (gene in geneToLength)


def gene_to_length(gene, filename='geneToLength_all_firstseen.txt',
                          median_gene_length=854):
    try:
        geneToLength = gene_to_length.geneToLength
    except AttributeError:
        gene_to_length.geneToLength = load_gene_lengths(filename)
        geneToLength = gene_to_length.geneToLength

    if gene in geneToLength:
        return geneToLength[gene]
    else:
        #print "Length of ", gene, "not found in file. Using median length instead:", median_gene_length
        return median_gene_length

# def gene_to_length_gainloss(gene, filename='/Users/jlu96/conte/jlu/geneToLength_gainloss_firstseen.txt',
#                           median_gene_length=854):
#     try:
#         geneToLength = gene_to_length.geneToLength
#     except AttributeError:
#         gene_to_length.geneToLength = load_gene_lengths(filename)
#         geneToLength = gene_to_length.geneToLength
#
#     if gene in geneToLength:
#         return geneToLength[gene]
#     else:
#         #print "Length of ", gene, "not found in file. Using median length instead:", median_gene_length
#         return median_gene_length


def prod(iterable):
    return reduce(operator.mul, iterable, 1)