__author__ = 'jlu96'

import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
import operator
import itertools
import time
import mutexprob as mep
import scorecooccur as sco
import numpy as np
import collections
import functools
from scipy.misc import comb
import csv
import bingenesbypairs as bgbp


def aremutex(genes, geneToCases, maxOverlap=0):
    """ONLY WORKS FOR SIZE TWO RIGHT NOW."""

    # set method
    setlist = [set(geneToCases[gene]) for gene in genes]
    overlap = len(set.intersection(*setlist))
    return (overlap <= maxOverlap)


def cooccur(genes, geneToCases, minCooccur=1, pairwise=False):
    if (pairwise):
        for gene1 in genes:
            for gene2 in genes[genes.index(gene1) + 1:]:
                setlist = [set(geneToCases[gene1]), set(geneToCases[gene2])]
                if (len(set.intersection(*setlist)) < minCooccur):
                    return False
        return True
    else:
        setlist = [set(geneToCases[gene]) for gene in genes]
        return len(set.intersection(*setlist)) >= minCooccur

        # def cooccurrence_ratio(overlap, *mutation_frequencies):
        #     return overlap * 1.0 / (sum(mutation_frequencies) - overlap)
        # cooccurrence_ratio = overlap * 1.0 / prod(mutationfrequencies)


def cooccurrence_ratio(genes, geneToCases):
    #print "Coverage is", numcoverage(genes, geneToCases), "for genes", genes, "of lengths", [len(geneToCases[gene]) for gene in genes]
    return numoverlaps(genes, geneToCases) * 1.0 / numcoverage(genes, geneToCases)


def cooccur_array(minCooccur=1, *caselists):
    # Assumes that the genelists are provided as numpy arrays.
    # THIS IS PRETTY INEFFICIENT. DO NOT USE!
    # print caselists
    overlap_limit = len(caselists) - 1
    sumarray = sum(caselists)
    return (sum([i > overlap_limit for i in sumarray]) >= minCooccur)
    # numOverlaps = 0
    # for entry in sumarray:
    #     if (entry > overlap_limit): #check if there are at least as many overlaps
    #         numOverlaps += 1
    #         if (numOverlaps >= minCooccur):
    #             break
    # return (numOverlaps >= minCooccur)


# Left off here. Now, you want to use just a binary array to make the number of cooccurring sets faster.-jlu
# Make the actual mutation matrix with gene mapping to dictionary. Have you done this somewhere already?

def allcooccursets(genes, patients, geneToCases, set_size, minCooccur=1):
    numsets = comb(len(genes), set_size)
    print "Total number of sets to test", numsets
    t = time.time()


    # Sort by the most recurrent genes
    genes.sort(key=lambda gene: len(geneToCases[gene]), reverse=True)

    gene_set_list = itertools.combinations(genes, set_size)
    list = []
    num_tested = 0
    num_found = 0

    for geneset in gene_set_list:
        if cooccur(geneset, geneToCases, minCooccur=minCooccur, pairwise=False):
            list.append(geneset)
            num_found += 1
            print "Working list:", geneset
            print "Number tested: ", num_tested, "\tTime: ", time.time() - t
            print "Number found: ", num_found
        num_tested += 1

    return list


    #
    # caseMatrix = {}
    # for gene in genes:
    #     caseMatrix[gene] = np.array([int(case in geneToCases[gene]) for case in patients])
    #
    # gene_set_list = itertools.combinations(genes, set_size)
    # list = []
    # num_tested = 0
    # num_found = 0
    #
    #
    # for geneset in gene_set_list:
    #     caselists = [caseMatrix[gene] for gene in geneset]
    #     if cooccur_array(minCooccur, *caselists):
    #         list.append(geneset)
    #         num_found += 1
    #         print "Working list:", geneset
    #         print "Number tested: ", num_tested, "\tTime: ", time.time() - t
    #         print "Number found: ", num_found

    # if cooccur(geneset, geneToCases, minCooccur = 20, pairwise = False):
    #     list.append(geneset)
    #     num_found += 1
    #     print "Number tested: ", num_tested, "\tTime: ", time.time() - t
    #     print "Number found: ", num_found
    # num_tested +=1

    return list


def makematrix(genes, patients, geneToCases):
    matrix = {}
    numCases = len(patients)
    for gene in genes:
        row = np.zeros(numCases)
        for i in range(numCases):
            if (patients[i] in geneToCases[gene]):
                row[i] = 1
        matrix[gene] = row
    return matrix


def numcoverage(genes, geneToCases):
    setlist = [set(geneToCases[gene]) for gene in genes]
    return len(set.union(*setlist))


def numoverlaps(genes, geneToCases):
    setlist = [set(geneToCases[gene]) for gene in genes]
    return len(set.intersection(*setlist))


def overlaps(genes, geneToCases):
    setlist = [set(geneToCases[gene]) for gene in genes]
    return set.intersection(*setlist)


def exclusive_dict(genes, geneToCases, patientToGenes):
    exclusive_dict = {}
    for gene in genes:
        exclusive_dict[gene] = set()

    for gene in genes:
        for case in geneToCases[gene]:
            if len(set.intersection(set(genes), set(patientToGenes[case]))) == 1:
                exclusive_dict[gene].add(case)
            # If patient A is in exclusive_dict of gene B, then patient A has only gene B, and no other genes in the set
            # of genes.

    return exclusive_dict


def numexclusive(genes, geneToCases, patientToGenes):
    exclusiveDict = exclusive_dict(genes, geneToCases, patientToGenes)
    return sum([len(exclusiveCases) for exclusiveCases in exclusiveDict.values()])




    # # forloop method DO NOT USE, OVERLAP IS WRONG
    # for gene1 in genes:
    #     for gene2 in genes[genes.index(gene1) + 1 :]:
    #         for case in geneToCases[gene1]:
    #             if case in geneToCases[gene2]:
    #                 return False
    # return True



    # WORKING CODE
    #  freq_dict = {}
    # for gene in genes:
    #     for case in geneToCases[gene]:
    #         if case in freq_dict:
    #             freq_dict[case] += 1
    #         else:
    #             freq_dict[case] = 1
    #
    # for frequency in freq_dict.itervalues():
    #     if (frequency > 1):
    #         maxOverlap -= 1
    #         if (maxOverlap < 0):
    #             return False
    # return True


def prod(iterable):
    return reduce(operator.mul, iterable, 1)


def getmutexsets(numCases, genes, geneToCases, patientToGenes, set_size, p, maxOverlap=0, genes1=None, genes2=None):
    # Sort by the most recurrent genes
    genes.sort(key=lambda gene: len(geneToCases[gene]), reverse=True)

    # A dictionary that gives the number of significant mutually exclusive sets the gene is in, and the gene
    # frequency, and the total coverage of each set that it's found to be in.
    # gene_occurrences = {}
    # for gene in genes:
    #     gene_occurrences[gene] = [0, len(geneToCases[gene]), 0, sco.gene_to_length(gene)]

    gene_occurrences = {}
    for gene in genes:
        gene_occurrences[gene] = collections.OrderedDict()
        gene_occurrences[gene]['Gene'] = gene
        gene_occurrences[gene]['NumberSets'] = 0
        gene_occurrences[gene]['MutationFrequency'] = len(geneToCases[gene])
        gene_occurrences[gene]['TotalCoverage'] = 0
        gene_occurrences[gene]['GeneLength'] = sco.gene_to_length(gene)


    tested_sets = set()

    # Each entry in sigsets is of the form [geneset, mep.mutexprob(numCases, mutationfrequencies, overlap),
    #                                     mutationfrequencies, overlap]
    sigsets = []

    if (set_size == 2):
        for gene1 in (genes1 if genes1 else genes):
            # Compare across lists this way.
            for gene2 in (genes2 if genes2 else genes[genes.index(gene1) + 1:]):
                geneset = [gene1, gene2]
                if gene1 != gene2 and frozenset(geneset) not in tested_sets:

                    tested_sets.add(frozenset(geneset))

                    mutationfrequencies, overlap, coverage, c_ratio, mutexprob = analyze_mutex_set(numCases, geneToCases,
                                                                                                   patientToGenes, geneset)
                    total_gene_length = sum([sco.gene_to_length(gene) for gene in geneset])
                    if mutexprob < p and overlap <= maxOverlap:
                        entry = collections.OrderedDict()
                        entry['GeneSet'] = geneset
                        entry['Probability'] = mutexprob
                        entry['MutationFrequencies'] = mutationfrequencies
                        entry['Overlap'] = overlap
                        entry['Coverage'] = coverage
                        entry['TotalGeneLength'] = total_gene_length

                        sigsets.append(entry)

                        for gene in geneset:
                            gene_occurrences[gene]['NumberSets'] +=1
                            gene_occurrences[gene]['TotalCoverage'] += coverage - len(geneToCases[gene])

    # Sort the sets with lowest pvalue first.
    sigsets.sort(key=lambda entry: entry['Probability'])

    # Remove all those genes that don't occur in any mutually exclusive sets.
    for gene in genes:
        if (gene_occurrences[gene]['NumberSets'] <= 0):
            gene_occurrences.pop(gene, None)

    # Sort the genes by those that appear in the most sets and those that have high frequency.
    #sorted_genes = sorted(gene_occurrences.items(), key=lambda entry: entry[1]['NumberSets'] * entry[1]['MutationFrequency'], reverse=True)

    return sigsets, gene_occurrences




# Add mutex_triangles here 7/4/15 -jlu




def getcooccursets_new(numCases, genes, geneToCases, patientToGenes, set_size, p, minCooccur=1,
                       min_cooccurrence_ratio=0.0,
                       trials=10000, just_sets=False,
                       genes1=None,
                       genes2=None):
    # Sort by the most recurrent genes
    genes.sort(key=lambda gene: len(geneToCases[gene]), reverse=True)

    # A dictionary that gives the number of significant cooccuring sets the gene is in, and the gene
    # frequency
    # gene_occurrences = {}
    # for gene in genes:
    #     gene_occurrences[gene] = [0, len(geneToCases[gene]), 0, sco.gene_to_length(gene)]

    gene_occurrences = {}
    for gene in genes:
        gene_occurrences[gene] = collections.OrderedDict()
        gene_occurrences[gene]['Gene'] = gene
        gene_occurrences[gene]['NumberSets'] = 0
        gene_occurrences[gene]['MutationFrequency'] = len(geneToCases[gene])
        gene_occurrences[gene]['TotalCoverage'] = 0
        gene_occurrences[gene]['GeneLength'] = sco.gene_to_length(gene)


    #Keep track of the ones that have been tested.
    tested_sets = set()

    # An entry in the isgset is: [geneset, mep.cooccurprob(numCases, mutationfrequencies, overlap),
    #                                      mutationfrequencies, overlap, cooccurrence ratio]
    cooccur_gene_dict = {}
    for gene in genes:
        cooccur_gene_dict[gene] = set()

    sigsets = []

    # Find the significant pairs first.
    for gene1 in (genes1 if genes1 else genes):
        for gene2 in (genes2 if genes2 else genes[genes.index(gene1) + 1:]):
            geneset = [gene1, gene2]
            if gene1 != gene2 and frozenset(geneset) not in tested_sets:

                tested_sets.add(frozenset(geneset))

                if numoverlaps(geneset, geneToCases) >= minCooccur:
                    mutationfrequencies, overlap, coverage, c_ratio, cooccurprob, \
                    score, overlap_pmn, average_overlap_pmn, combined_score = analyze_cooccur_set(numCases, geneToCases, patientToGenes,
                                                                                  geneset)
                    total_gene_length = sum([sco.gene_to_length(gene) for gene in geneset])


                    if cooccurprob < p and overlap >= minCooccur and c_ratio >= min_cooccurrence_ratio:
                        if (set_size == 2):
                            entry = collections.OrderedDict()
                            entry['GeneSet'] = geneset
                            entry['Probability'] = cooccurprob
                            entry['MutationFrequencies'] = mutationfrequencies
                            entry['Overlap'] = overlap
                            entry['CooccurrenceRatio'] = c_ratio
                            entry['Coverage'] = coverage
                            entry['SetScore'] = score
                            entry['AverageOverlapPMN'] = average_overlap_pmn
                            entry['TotalGeneLength'] = total_gene_length
                            entry['CombinedScore'] = combined_score

                            sigsets.append(entry)
                            for gene in geneset:
                                gene_occurrences[gene]['NumberSets'] +=1
                                gene_occurrences[gene]['TotalCoverage'] += coverage - len(geneToCases[gene])

                        geneset_set = set(geneset)
                        for gene in geneset:
                            # geneset_set.remove(gene)
                            # Add the rest of the gene elements to the gene's dictionary of cooccurring genes.
                            cooccur_gene_dict[gene] = cooccur_gene_dict[gene].union(geneset_set.difference(set([gene])))

    if set_size == 3:
        possible_sets = []

        for gene_to_genes in cooccur_gene_dict.items():
            cooccur_pair_list = itertools.combinations(gene_to_genes[1], set_size - 1)
            for cooccur_pair in cooccur_pair_list:
                cooccur_pair_set = set(cooccur_pair)
                elem = cooccur_pair_set.pop()
                if cooccur_pair_set.issubset(cooccur_gene_dict[elem]):
                    geneset = cooccur_pair_set
                    geneset.add(elem)
                    geneset.add(gene_to_genes[0])

                    if frozenset(geneset) not in tested_sets:
                        tested_sets.add(frozenset(geneset))

                        if numoverlaps(geneset, geneToCases) >= minCooccur:
                            mutationfrequencies, overlap, coverage, c_ratio, cooccurprob, \
                            score, overlap_pmn, average_overlap_pmn, combined_score = analyze_cooccur_set(numCases, geneToCases,
                                                                                          patientToGenes,
                                                                                          geneset,
                                                                                          compute_prob=(not just_sets))
                            total_gene_length = sum([sco.gene_to_length(gene) for gene in geneset])
                            #
                            # entry = [geneset, cooccurprob,
                            #          mutationfrequencies, overlap, c_ratio, coverage, score,
                            #          average_overlap_pmn, total_gene_length, combined_score]

                            entry = collections.OrderedDict()
                            entry['GeneSet'] = geneset
                            entry['Probability'] = cooccurprob
                            entry['MutationFrequencies'] = mutationfrequencies
                            entry['Overlap'] = overlap
                            entry['CooccurrenceRatio'] = c_ratio
                            entry['Coverage'] = coverage
                            entry['SetScore'] = score
                            entry['AverageOverlapPMN'] = average_overlap_pmn
                            entry['TotalGeneLength'] = total_gene_length
                            entry['CombinedScore'] = combined_score

                            if overlap >= minCooccur and c_ratio >= min_cooccurrence_ratio:
                                possible_sets.append(geneset)
                                if not just_sets:
                                    if cooccurprob < p:
                                        sigsets.append(entry)
                                else:
                                    sigsets.append(entry)
                                for gene in geneset:
                                    gene_occurrences[gene]['NumberSets'] +=1
                                    gene_occurrences[gene]['TotalCoverage'] += coverage - len(geneToCases[gene])

        print "Number of possible cooccuring sets", len(possible_sets)

        print "Number of significant sets of size 3", len(sigsets)


        # Sort the sets with highest product first.
    #    sigsets.sort(key=lambda entry: prod(entry[2]), reverse=True)

    # Sort the sets with lowest pvalue first.
    sigsets.sort(key=lambda entry: entry['Probability'])

    # Remove all those genes that don't occur in any cooccurring sets.
    for gene in genes:
        if (gene_occurrences[gene]['NumberSets'] <= 0):
            gene_occurrences.pop(gene, None)

    # Sort the genes by those that appear in the most sets and those that have high frequency.
    #sorted_genes = sorted(gene_occurrences.items(), key=lambda entry: entry[1]['NumberSets'] * entry[1]['MutationFrequency'], reverse=True)

    return sigsets, gene_occurrences


def mutexpvalues(numCases, genes, geneToCases, set_size, maxOverlap=0):
    p_values = []

    frequency_dict = {}
    for gene in genes:
        frequency_dict[gene] = len(geneToCases[gene])

    genes.sort(key=lambda gene: frequency_dict[gene], reverse=True)

    gene_set_list = itertools.combinations(genes, set_size)

    num_tested = 0
    # ismutextime = 0.0
    # probtime = 0.0
    t = time.time()
    for geneset in gene_set_list:
        mutationfrequencies = [frequency_dict[gene] for gene in geneset]
        overlap = numoverlaps(geneset, geneToCases)
        if (sum(mutationfrequencies) < numCases):

            # t1 = time.time()
            # ismutex = aremutex(geneset, geneToCases, maxOverlap = maxOverlap)
            # t2 = time.time()

            if aremutex(geneset, geneToCases, maxOverlap=maxOverlap):
                p_values.append(mep.mutexprob(numCases, mutationfrequencies, overlap))
                num_tested += 1
                if (num_tested % 1000000 == 0):
                    print "Num tested: ", num_tested, " time: ", time.time() - t
                    t = time.time()
                    # t3 = time.time()

                    # ismutextime += t2 - t1
                    # mep.probtime += t3 - t2
                    #
                    # print "ismutextime: ", ismutextime
                    # print "mep.probtime: ", mep.probtime

    return p_values


def cooccurpvalues(numCases, genes, geneToCases, set_size, minCooccur=1):
    p_values = []

    frequency_dict = {}
    for gene in genes:
        frequency_dict[gene] = len(geneToCases[gene])

    genes.sort(key=lambda gene: frequency_dict[gene], reverse=True)

    gene_set_list = itertools.combinations(genes, set_size)

    num_tested = 0
    # ismutextime = 0.0
    # probtime = 0.0
    t = time.time()
    for geneset in gene_set_list:
        mutationfrequencies = [frequency_dict[gene] for gene in geneset]
        overlap = numoverlaps(geneset, geneToCases)

        # t1 = time.time()
        # ismutex = aremutex(geneset, geneToCases, maxOverlap = maxOverlap)
        # t2 = time.time()

        if cooccur(geneset, geneToCases, minCooccur=minCooccur):
            p_values.append(mep.cooccurprob(numCases, mutationfrequencies, overlap))
            num_tested += 1
            if (num_tested % 1000000 == 0):
                print "Num tested: ", num_tested, " time: ", time.time() - t
                t = time.time()
                # t3 = time.time()

                # ismutextime += t2 - t1
                # mep.probtime += t3 - t2
                #
                # print "ismutextime: ", ismutextime
                # print "mep.probtime: ", mep.probtime

    return p_values


def plotmutexdistribution(numCases, genes, geneToCases, set_size, maxOverlap=0, bins=100):
    p_values = mutexpvalues(numCases, genes, geneToCases, set_size, maxOverlap=maxOverlap)
    plt.figure()
    plt.hist(p_values, bins)
    plt.xlabel('pvalue')
    plt.ylabel('Frequency')
    plt.show()


def plotcooccurdistribution(numCases, genes, geneToCases, set_size, minCooccur=1, bins=100):
    p_values = cooccurpvalues(numCases, genes, geneToCases, set_size, minCooccur=minCooccur)
    plt.figure()
    plt.hist(p_values, bins)
    plt.xlabel('pvalue')
    plt.ylabel('Frequency')
    plt.show()


def runmutex(mutationmatrix, geneFile=None, patientFile=None,
             minFreq=0, set_size=2, maxOverlap=0, p_value=0.05):
    """
      * **m** (*int*) - number of patients.
      * **n** (*int*) - number of genes.
      * **genes** (*list*) - genes in the mutation data.
      * **patients** (*list*) - patients in the mutation data.
      * **geneToCases** (*dictionary*) - mapping of genes to the patients in which they are mutated.
      * **patientToGenes** (*dictionary*) - mapping of patients to the genes they have mutated.
      """

    # Load the mutation data
    mutations = load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)
    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mutations

    sigsets, gene_occurrences = getmutexsets(numCases, genes, geneToCases, patientToGenes, set_size, p_value,
                                             maxOverlap)

    return sigsets, gene_occurrences


# def runcooccur(mutationmatrix, geneFile=None, patientFile=None,
#                minFreq=0, set_size=2, minCooccur=1, p_value=0.05):
#     """
#       * **m** (*int*) - number of patients.
#       * **n** (*int*) - number of genes.
#       * **genes** (*list*) - genes in the mutation data.
#       * **patients** (*list*) - patients in the mutation data.
#       * **geneToCases** (*dictionary*) - mapping of genes to the patients in which they are mutated.
#       * **patientToGenes** (*dictionary*) - mapping of patients to the genes they have mutated.
#       """
#
#     # Load the mutation data
#     mutations = load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)
#     numGenes, numCases, genes, patients, geneToCases, patientToGenes = mutations
#
#     sigsets, gene_occurrences = getcooccursets(numCases, genes, geneToCases, set_size, p_value, minCooccur)
#     return sigsets, gene_occurrences


def load_mutation_data(filename, patientFile=None, geneFile=None, minFreq=0):
    """Loads the mutation data in the given file.

    :type filename: string
    :param filename: Location of mutation data file.
    :type patient_file: string
    :param patient_file: Location of patient (whitelist) file.
    :type gene_file: string
    :param gene_file: Location of gene (whitelist) file.
    :rtype: Tuple

    **Returns:**
      * **m** (*int*) - number of genes.
      * **n** (*int*) - number of patients.
      * **genes** (*list*) - genes in the mutation data.
      * **patients** (*list*) - patients in the mutation data.
      * **geneToCases** (*dictionary*) - mapping of genes to the patients in which they are mutated.
      * **patientToGenes** (*dictionary*) - mapping of patients to the genes they have mutated.
    """
    # Load the whitelists (if applicable)
    if patientFile:
        with open(patientFile) as f:
            patients = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
    else:
        patients = None

    if geneFile:
        with open(geneFile) as f:
            genes = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
    else:
        genes = set()

    # Parse the mutation matrix
    from collections import defaultdict

    geneToCases, patientToGenes = defaultdict(set), defaultdict(set)
    with open(filename) as f:
        arrs = [l.rstrip().split("\t") for l in f if not l.startswith("#")]
        for arr in arrs:
            patient, mutations = arr[0], set(arr[1:])

            if not patients or patient in patients:
                if genes: mutations &= genes
                # else: genes |= mutations

                patientToGenes[patient] = mutations
                for gene in mutations:
                    geneToCases[gene].add(patient)

    genes = geneToCases.keys()
    # Remove genes with fewer than min_freq mutations
    toRemove = [g for g in genes if len(geneToCases[g]) < minFreq]
    for g in toRemove:
        for p in geneToCases[g]:
            patientToGenes[p].remove(g)
        del geneToCases[g]
        genes.remove(g)

    # Format and return output
    genes, patients = geneToCases.keys(), patientToGenes.keys()
    m, n = len(genes), len(patients)
    return m, n, genes, patients, geneToCases, patientToGenes


def remove_blacklists(gene_blacklist, patient_blacklist, numGenes, numCases, genes, patients, geneToCases,
                      patientToGenes):
    if gene_blacklist == None and patient_blacklist == None:
        return numGenes, numCases, genes, patients, geneToCases, patientToGenes
    else:
        if patient_blacklist:
            with open(patient_blacklist) as f:
                blacklist_patients = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
        else:
            blacklist_patients = None

        if gene_blacklist:
            with open(gene_blacklist) as f:
                blacklist_genes = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
        else:
            blacklist_genes = set()

        newgeneToCases = geneToCases
        newpatientToGenes = patientToGenes
        if blacklist_genes:
            for item in geneToCases.items():
                if item[0] in blacklist_genes:
                    for patient in item[1]:
                        newpatientToGenes[patient].remove(item[0])
                        if not newpatientToGenes[patient]:
                            newpatientToGenes.pop(patient)
                    newgeneToCases.pop(item[0])
                    print "gene ", item[0], ' popped'
        if blacklist_patients:
            for item in patientToGenes.items():
                if item[0] in blacklist_patients:
                    for gene in item[1]:
                        newgeneToCases[gene].remove(item[0])
                        if not newgeneToCases[gene]:
                            newgeneToCases.pop(gene)
                    newpatientToGenes.pop(item[0])
                    print "patient ", item[0], ' popped'

        newgenes = newgeneToCases.keys()
        newpatients = newpatientToGenes.keys()
        newnumGenes = len(newgenes)
        newnumCases = len(newpatients)
        return newnumGenes, newnumCases, newgenes, newpatients, newgeneToCases, newpatientToGenes

def remove_extra_genes(extra_genes, numGenes, numCases, genes, patients, geneToCases, patientToGenes):
    newgeneToCases = geneToCases
    newpatientToGenes = patientToGenes

    if extra_genes:
        for item in geneToCases.items():
            if item[0] in extra_genes:
                for patient in item[1]:
                    newpatientToGenes[patient].remove(item[0])
                    if not newpatientToGenes[patient]:
                        newpatientToGenes.pop(patient)
                newgeneToCases.pop(item[0])
    newgenes = newgeneToCases.keys()
    newpatients = newpatientToGenes.keys()
    newnumGenes = len(newgenes)
    newnumCases = len(newpatients)
    return newnumGenes, newnumCases, newgenes, newpatients, newgeneToCases, newpatientToGenes



def graph_mutation_distribution(numCases, genes, geneToCases, filename, top_percentile=10, bins=100):
    mutationfrequencies = [len(geneToCases[gene]) for gene in genes]
    plt.figure()
    plt.hist(mutationfrequencies, bins)
    percentile_score = stats.scoreatpercentile(mutationfrequencies, 100 - top_percentile)
    print "Score is ", percentile_score
    plt.axvline(percentile_score, color='k', linestyle='dashed', linewidth=2,
                label='Top ' + str(top_percentile) + ' Percentile at ' + str(percentile_score))
    plt.title(filename.split('.m2')[0] + "'s Distribution of Mutation Frequencies, n = " + str(numCases), fontsize=20)
    plt.xlabel('Mutation Frequency (# Samples with Mutations)', fontsize=20)
    plt.ylabel('Amount of Mutations', fontsize=20)
    plt.legend()
    # plt.show()
    plt.savefig(filename + '_' + str(bins) + '.pdf', bbox_inches='tight')


def graph_patient_distribution(numGenes, cases, patientToGenes, filename, top_percentile=10, bins=100):
    patientfrequencies = [len(patientToGenes[case]) for case in cases]
    plt.figure()
    plt.hist(patientfrequencies, bins)
    percentile_score = stats.scoreatpercentile(patientfrequencies, 100 - top_percentile)
    # plt.axvline(percentile_score, color='k', linestyle='dashed', linewidth=2,
    #             label='Top ' + str(top_percentile) + ' Percentile at ' + str(percentile_score))
    # plt.title(
    #     filename.split('.m2')[0] + "'s Distribution of Mutation Number per Patient, Total Genes  = " + str(numGenes),
    #     fontsize=20)
    plt.xlabel('# Somatic Mutations In Tumor', fontsize=20)
    plt.ylabel('Number of Samples', fontsize=20)
    plt.legend()
    plt.show()


def graph_cooccurrenceratio_distribution(csigsets, filename, top_percentile=10):
    c_ratios = [entry['CooccurrenceRatio'] for entry in csigsets]
    plt.figure()
    plt.hist(c_ratios, 100)
    percentile_score = stats.scoreatpercentile(c_ratios, 100 - top_percentile)
    # print "Score is ", percentile_score
    plt.axvline(percentile_score, color='k', linestyle='dashed', linewidth=2,
                label='Top ' + str(top_percentile) + ' Percentile at ' + str(percentile_score))
    plt.title(filename.split('.m2')[0] + "'s Distribution of Cooccurrence Ratios, " + str(len(c_ratios)) +
              " sets, Size = " + str(len(csigsets[0]['GeneSet'])), fontsize=20)
    plt.xlabel('Cooccurrence Ratio', fontsize=20)
    plt.ylabel('Frequency of sets', fontsize=20)
    plt.legend()
    plt.show()


def writemutexsets(file, sigsets, gene_occurrences, p_value):
    import decimal

    file.write(str(len(sigsets)) + ' Significant Mutually Exclusive Gene Sets at p = ' + str(p_value) + '\n')

    if sigsets:
        keys = sigsets[0].keys()

        writer = csv.DictWriter(file, delimiter='\t', fieldnames=keys)
        writer.writeheader()
        for entry in sigsets:
            writer.writerow(entry)
        #
        # file.write('Gene Set\tProbability\tMutation Frequencies\tOverlap\tCoverage\tTotalGeneLength\n')
        # for entry in sigsets:plo
        #     file.write(
        #         str(entry[0]) + '\t' + str(round(entry[1], abs(decimal.Decimal(str(p_value)).as_tuple().exponent) + 2))
        #         + '\t' + str(entry[2]) + '\t' + str(entry[3]) + '\t' + str(entry[4]) + '\t' + str(entry[5]) + '\n')

        # Sorted by highest coverage.
        sigsets.sort(key=lambda entry: entry['Coverage'], reverse=True)
        file.write('\nSets sorted by highest coverage.\n')
        writer.writeheader()
        for entry in sigsets:
            writer.writerow(entry)


        #
        # file.write('Gene Set\tProbability\tMutation Frequencies\tOverlap\tCoverage\tTotalGeneLength\n')
        # for entry in sigsets:
        #     file.write(
        #         str(entry[0]) + '\t' + str(round(entry[1], abs(decimal.Decimal(str(p_value)).as_tuple().exponent) + 2))
        #         + '\t' + str(entry[2]) + '\t' + str(entry[3]) + '\t' + str(entry[4]) + '\t' + str(entry[5]) + '\n')



        # DENISOVICH
        # Sorted by maximum number of mutually exclusive sets
        file.write('\nList of ' + str(len(gene_occurrences)) + ' Genes, sorted by descending product of number' + \
                   'of mutually exclusive sets and frequency\n')

        keys = gene_occurrences[gene_occurrences.keys()[0]].keys()
        writer = csv.DictWriter(file, delimiter='\t', fieldnames=keys)
        writer.writeheader()
        sorted_keys = sorted(gene_occurrences.keys(), key=lambda gene: gene_occurrences[gene]['NumberSets'] * gene_occurrences[gene]['MutationFrequency'], reverse=True)
        for entry in sorted_keys:
            writer.writerow(gene_occurrences[entry])


def writemutex_cytoscape(file_prefix, sigsets, gene_occurrences, p_value):
    filename = file_prefix + '_cytoscape.tsv'
    opener = open(filename, 'w')
    set_size = len(sigsets[0]['GeneSet'])
    # opener.write(str(len(sigsets)) +' Significant Mutually Exclusive Gene Sets at p = ' + str(p_value) + '\n')

    newgenekeys = ['Gene' + str(i) for i in range(set_size)]
    fieldnames = sigsets[0].keys()
    fieldnames.remove('GeneSet')
    fieldnames.remove('MutationFrequencies')
    fieldnames = newgenekeys + fieldnames


    writer = csv.DictWriter(opener, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
    writer.writeheader()

    for entry in sigsets:
        for i in range(set_size):
            entry[newgenekeys[i]] = entry['GeneSet'][i]
        writer.writerow(entry)

    opener.close()


    mutationfrequencyfile = open(file_prefix + '_cytoscape_freq.noa', 'w')
    keys = gene_occurrences[gene_occurrences.keys()[0]].keys()
    writer = csv.DictWriter(mutationfrequencyfile, delimiter='\t', fieldnames=keys)
    writer.writeheader()
    for entry in gene_occurrences:
        writer.writerow(gene_occurrences[entry])

    mutationfrequencyfile.close()


def writecooccur_cytoscape(file_prefix, sigsets, gene_occurrences, p_value):
    filename = file_prefix + '_cytoscape.tsv'
    opener = open(filename, 'w')
    set_size = len(sigsets[0]['GeneSet'])

    newgenekeys = ['Gene' + str(i) for i in range(set_size)]
    fieldnames = sigsets[0].keys()
    fieldnames.remove('GeneSet')
    fieldnames.remove('MutationFrequencies')
    fieldnames = newgenekeys + fieldnames

    writer = csv.DictWriter(opener, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
    writer.writeheader()

    for entry in sigsets:
        for i in range(set_size):
            entry[newgenekeys[i]] = entry['GeneSet'][i]
        writer.writerow(entry)

    opener.close()

    #
    # mutationfrequencyfile = open(file_prefix + '_cytoscape_freq.noa', 'w')
    # keys = gene_occurrences[0][1].keys()
    # writer = csv.DictWriter(mutationfrequencyfile, delimiter='\t', fieldnames=keys)
    # writer.writeheader()
    # # file.write('Gene\tMutuallyExclusiveSets\tMutationFrequency\tTotalCoverage\tGeneLength\n')
    # for entry in gene_occurrences:
    #     writer.writerow(entry[1])
    # mutationfrequencyfile.close()


    mutationfrequencyfile = open(file_prefix + '_cytoscape_freq.noa', 'w')
    keys = gene_occurrences[gene_occurrences.keys()[0]].keys()
    writer = csv.DictWriter(mutationfrequencyfile, delimiter='\t', fieldnames=keys)
    writer.writeheader()
    for entry in gene_occurrences:
        writer.writerow(gene_occurrences[entry])

    mutationfrequencyfile.close()


def writecooccursets(file, sigsets, gene_occurrences, p_value):
    import decimal

    file.write(str(len(sigsets)) + ' Significant Cooccuring Gene Sets at p = ' + str(p_value) + '\n')

    keys = sigsets[0].keys()

    writer = csv.DictWriter(file, delimiter='\t', fieldnames=keys)
    writer.writeheader()
    for entry in sigsets:
        writer.writerow(entry)

    file.write('\nThe same sets, sorted by Cooccurrence ratio\n')

    writer.writeheader()
    sigsets.sort(key=lambda entry: entry['CooccurrenceRatio'], reverse=True)
    for entry in sigsets:
        writer.writerow(entry)

    # Sorted by maximum number of cooccur sets and product
    file.write('\nList of ' + str(len(gene_occurrences)) + ' Genes, sorted by descending product of number' + \
               'of mutually exclusive sets and frequency\n')
    # keys = gene_occurrences[0][1].keys()
    # writer = csv.DictWriter(file, delimiter='\t', fieldnames=keys)
    # writer.writeheader()
    # for entry in gene_occurrences:
    #     writer.writerow(entry[1])
    #

    keys = gene_occurrences[gene_occurrences.keys()[0]].keys()
    writer = csv.DictWriter(file, delimiter='\t', fieldnames=keys)
    writer.writeheader()
    sorted_keys = sorted(gene_occurrences.keys(), key=lambda gene: gene_occurrences[gene]['NumberSets'] * gene_occurrences[gene]['MutationFrequency'], reverse=True)
    for entry in sorted_keys:
        writer.writerow(gene_occurrences[entry])

def writecompleteinfo(filename, sigsets, gene_occurrences):

    opener = open(filename, 'w')
    set_size = len(sigsets[0]['GeneSet'])

    newmutationkeys = ['MutationFrequency' + str(i) for i in range(set_size)]
    newsetskeys = ['NumberSets' + str(i) for i in range(set_size)]
    newgenelengthkeys= ['GeneLength' + str(i) for i in range(set_size)]
    otherkeys = ['MutationFrequencyDifference', 'GeneLengthDifference']

    fieldnames = sigsets[0].keys()
    fieldnames.remove('GeneSet')
    fieldnames.remove('MutationFrequencies')
    fieldnames = fieldnames + newmutationkeys + newsetskeys + newgenelengthkeys + otherkeys

    writer = csv.DictWriter(opener, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
    writer.writeheader()

    for entry in sigsets:
        for i in range(set_size):
            gene = entry['GeneSet'][i]
            entry[newmutationkeys[i]] = entry['MutationFrequencies'][i]
            entry[newsetskeys[i]] = gene_occurrences[gene]['NumberSets']
            entry[newgenelengthkeys[i]] = gene_occurrences[gene]['GeneLength']
        entry['MutationFrequencyDifference'] = abs(entry['MutationFrequencies'][0] - entry['MutationFrequencies'][1])
        entry['GeneLengthDifference'] = abs(entry[newgenelengthkeys[0]] - entry[newgenelengthkeys[1]])


        writer.writerow(entry)

    opener.close()



    #
    #
    # file.write('\nList of  ' + str(len(gene_occurrences)) + ' Genes, sorted by descending product of number of ' + \
    #            'cooccurring sets and frequency\n')
    # file.write('Gene\tCooccurringSets\tMutationFrequency\tTotalCoverage\tGeneLength\n')
    # for entry in gene_occurrences:
    #     file.write(entry[0] + '\t' + str(entry[1][0]) + '\t' + str(entry[1][1]) + '\t' + str(entry[1][2]) +
    #                 '\t' + str(entry[1][3]) + '\n')

def filtermatrix(top_percentile, numGenes, numCases, genes, patients, geneToCases, patientToGenes):
    mutationfrequencies = [len(geneToCases[gene]) for gene in genes]
    cutoff = stats.scoreatpercentile(mutationfrequencies, 100 - top_percentile)
    newgeneToCases = geneToCases
    newpatientToGenes = patientToGenes
    for item in geneToCases.items():
        if len(item[1]) < cutoff:
            # for patient in item[1]:
            #     newpatientToGenes[patient].remove(item[0])
            #     if not newpatientToGenes[patient]:
            #         newpatientToGenes.pop(patient)
            newgeneToCases.pop(item[0])

    newgenes = newgeneToCases.keys()
    newpatients = newpatientToGenes.keys()
    newnumGenes = len(newgenes)
    newnumCases = len(newpatients)
    return newnumGenes, newnumCases, newgenes, newpatients, newgeneToCases, newpatientToGenes


def filtertopnumbermatrix(top_number, numGenes, numCases, genes, patients, geneToCases, patientToGenes):
    mutationfrequencies = [len(geneToCases[gene]) for gene in genes]
    cutoff = stats.scoreatpercentile(mutationfrequencies, 100 - top_number * 100. / numGenes)
    newgeneToCases = geneToCases
    newpatientToGenes = patientToGenes
    for item in geneToCases.items():
        if len(item[1]) < cutoff:
            # for patient in item[1]:
            #     newpatientToGenes[patient].remove(item[0])
            #     # if not newpatientToGenes[patient]:
            #     #     newpatientToGenes.pop(patient)
            newgeneToCases.pop(item[0])

    newgenes = newgeneToCases.keys()
    newpatients = newpatientToGenes.keys()
    newnumGenes = len(newgenes)
    newnumCases = len(newpatients)
    return newnumGenes, newnumCases, newgenes, newpatients, newgeneToCases, newpatientToGenes


def getmutexstats(filename, *geneset, **kwargs):
    minFreq = kwargs['minFreq'] if 'minFreq' in kwargs else 0
    trials = kwargs['trials'] if 'trials' in kwargs else 10000

    mutations = load_mutation_data(filename, minFreq=minFreq)
    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mutations

    mutationfrequencies = [len(geneToCases[gene]) for gene in geneset]
    overlap = numoverlaps(geneset, geneToCases)
    coverage = numcoverage(geneset, geneToCases)
    c_ratio = overlap * 1.0 / coverage
    exclusive = numexclusive(geneset, geneToCases, patientToGenes)

    mutexprob = mep.mutexprob_approximate(numCases, exclusive, trials, *mutationfrequencies)
    cooccurprob = mep.cooccurprob_approximate(numCases, overlap, trials, *mutationfrequencies)

    print "For the genes: ", geneset
    print "Mutation frequencies are: ", mutationfrequencies
    print "Overlap is ", overlap
    print "Coverage is ", coverage
    print "Cooccurrence ratio is ", c_ratio
    print "Number of exclusive alterations is ", exclusive
    print "Mutual exclusivity p-value is ", mutexprob
    print "Cooccurrence p-value is ", cooccurprob

    return mutationfrequencies, overlap, coverage, c_ratio, exclusive, mutexprob, cooccurprob




def analyze_mutex_set(numCases, geneToCases, patientToGenes, geneset, compute_prob=True, trials=10000):
    # Mutation Frequencies
    mutationfrequencies = [len(geneToCases[gene]) for gene in geneset]
    overlap = numoverlaps(geneset, geneToCases)
    coverage = numcoverage(geneset, geneToCases)
    c_ratio = cooccurrence_ratio(geneset, geneToCases)

    if compute_prob:
        if len(geneset) == 2:
            mutexprob = mep.mutexprob(numCases, mutationfrequencies, overlap)
        else:
            exclusive = numexclusive(geneset, geneToCases, patientToGenes)
            mutexprob = mep.mutexprob_approximate(numCases, exclusive, trials, *mutationfrequencies)
    else:
        mutexprob = 0.0

    return mutationfrequencies, overlap, coverage, c_ratio, mutexprob

def analyze_mutex_set_new(numCases, geneToCases, patientToGenes, geneset, compute_prob=True, trials=10000):
    # Mutation Frequencies
    mutationfrequencies = [len(geneToCases[gene]) for gene in geneset]
    overlap = numoverlaps(geneset, geneToCases)
    coverage = numcoverage(geneset, geneToCases)
    c_ratio = overlap * 1.0 / coverage

    if compute_prob:
        if len(geneset) == 2:
            mutexprob = mep.mutexprob(numCases, mutationfrequencies, overlap)
        else:
            exclusive = numexclusive(geneset, geneToCases, patientToGenes)
            mutexprob = mep.mutexprob_approximate(numCases, exclusive, trials, *mutationfrequencies)
    else:
        mutexprob = 0.0

    # distance = bgbp.get_gene_distance(*geneset)


    mstats = {}
    mstats['GeneSet'] = geneset
    mstats['Gene0'], mstats['Gene1'] = tuple(geneset)
    mstats['Type'] = 'MutuallyExclusive'
    mstats['Probability'] = mutexprob
    mstats['MutationFrequencies'] = mutationfrequencies
    mstats['MutationFrequency0'], mstats['MutationFrequency1'] = tuple(mutationfrequencies)
    mstats['MutationFrequencyDifference'] = abs(mutationfrequencies[0] - mutationfrequencies[1])
    mstats['MutationFrequencyDifferenceRatio'] = mstats['MutationFrequencyDifference'] * 1.0 / coverage
    mstats['Overlap'] = overlap
    mstats['Coverage'] = coverage
    mstats['CooccurrenceRatio'] = c_ratio
    mstats['Concordance'] = (mstats['Gene0'][-4:] == mstats['Gene1'][-4:])

    mstats['Somatic'] = 0
    if (mstats['Gene0'][-4:] not in {'loss', 'gain'}):
        mstats['Somatic'] += 1
    if (mstats['Gene1'][-4:] not in {'loss', 'gain'}):
        mstats['Somatic'] += 1

    mstats['RoundedLogPCov'] = round(-np.log10(mutexprob/coverage), 1)
    # mstats['Weight'] = 20.0 -
    # mstats['Distance'] = distance

    return mstats

#
# def analyze_mutex_set_new_network(numCases, geneToCases, patientToGenes, geneset, compute_prob=True, trials=10000):
#     try:
#         analyze_mutex_set_new_network.pm
#     except AttributeError:
#         analyze_mutex_set_new_network.pm = mun.PermutationMatrices(geneToCases, patientToGenes, num_permutations=200)
#
#
#     # Mutation Frequencies
#     mutationfrequencies = [len(geneToCases[gene]) for gene in geneset]
#     overlap = numoverlaps(geneset, geneToCases)
#     coverage = numcoverage(geneset, geneToCases)
#     c_ratio = overlap * 1.0 / coverage
#
#     if compute_prob:
#         if len(geneset) == 2:
#             mutexprob = mep.mutexprob(numCases, mutationfrequencies, overlap)
#
#             # Calculate probability
#             condition = {}
#             condition['Genes'] = geneset
#             condition['Overlap'] = overlap
#             condition['Mutex'] = True
#
#             condition_function = mun.Condition([condition])
#             networkprob = analyze_mutex_set_new_network.pm.pvalue(condition_function)
#             if networkprob < mutexprob:
#                 print "Old mutexprob is ", mutexprob, "and newmutexprob is ", networkprob
#         #
#         # else:
#         #     exclusive = numexclusive(geneset, geneToCases, patientToGenes)
#         #     mutexprob = mep.mutexprob_approximate(numCases, exclusive, trials, *mutationfrequencies)
#     else:
#         mutexprob = 0.0
#
#     # distance = bgbp.get_gene_distance(*geneset)
#
#
#     mstats = {}
#     mstats['GeneSet'] = geneset
#     mstats['Gene0'], mstats['Gene1'] = tuple(geneset)
#     mstats['Type'] = 'MutuallyExclusive'
#     # mstats['Probability'] = mutexprob
#
#     mstats['Probability'] = networkprob
#
#     mstats['MutationFrequencies'] = mutationfrequencies
#     mstats['MutationFrequency0'], mstats['MutationFrequency1'] = tuple(mutationfrequencies)
#     mstats['MutationFrequencyDifference'] = abs(mutationfrequencies[0] - mutationfrequencies[1])
#     mstats['MutationFrequencyDifferenceRatio'] = mstats['MutationFrequencyDifference'] * 1.0 / coverage
#     mstats['Overlap'] = overlap
#     mstats['Coverage'] = coverage
#     mstats['CooccurrenceRatio'] = c_ratio
#     mstats['Concordance'] = (mstats['Gene0'][-4:] == mstats['Gene1'][-4:])
#
#     mstats['Somatic'] = 0
#     if (mstats['Gene0'][-4:] not in {'loss', 'gain'}):
#         mstats['Somatic'] += 1
#     if (mstats['Gene1'][-4:] not in {'loss', 'gain'}):
#         mstats['Somatic'] += 1
#
#     mstats['RoundedLogPCov'] = round(-np.log10(mutexprob/coverage), 1)
#     # mstats['Weight'] = 20.0 -
#     # mstats['Distance'] = distance
#
#     return mstats



def analyze_cooccur_set(numCases, geneToCases, patientToGenes, geneset, compute_prob=True, trials=10000):
    mutationfrequencies = [len(geneToCases[gene]) for gene in geneset]
    overlap = numoverlaps(geneset, geneToCases)

    if compute_prob:
        if len(geneset) == 2:
            cooccurprob = mep.cooccurprob(numCases, mutationfrequencies, overlap)
        else:
            cooccurprob = mep.cooccurprob_approximate(numCases, overlap, trials, *mutationfrequencies)
    else:
        cooccurprob = 0.0

    coverage = numcoverage(geneset, geneToCases)
    c_ratio = cooccurrence_ratio(geneset, geneToCases)

    if overlap:
        score = sco.score_cooccur_set(geneset, geneToCases, patientToGenes)

        overlap_pmn = [len(patientToGenes[patient]) for patient in overlaps(geneset, geneToCases)]
        average_overlap_pmn = sum(overlap_pmn) * 1.0 / len(overlap_pmn)
        combined_score = score / c_ratio
    else:
        score = 0
        overlap_pmn = 0
        average_overlap_pmn = 0
        combined_score = 0


    return mutationfrequencies, overlap, coverage, c_ratio, cooccurprob, score, overlap_pmn, average_overlap_pmn, combined_score


def analyze_cooccur_set_new(numCases, geneToCases, patientToGenes, geneset, compute_prob=True, trials=10000,
                            getdistance=False):
    mutationfrequencies = [len(geneToCases[gene]) for gene in geneset]
    overlap = numoverlaps(geneset, geneToCases)

    if compute_prob:
        if len(geneset) == 2:
            cooccurprob = mep.cooccurprob(numCases, mutationfrequencies, overlap)
        else:
            cooccurprob = mep.cooccurprob_approximate(numCases, overlap, trials, *mutationfrequencies)
    else:
        cooccurprob = 0.0

    coverage = numcoverage(geneset, geneToCases)
    c_ratio = overlap * 1.0 / coverage

    if overlap:
        score = sco.score_cooccur_set(geneset, geneToCases, patientToGenes)

        overlap_pmn = [len(patientToGenes[patient]) for patient in overlaps(geneset, geneToCases)]
        average_overlap_pmn = sum(overlap_pmn) * 1.0 / len(overlap_pmn)
        combined_score = score / c_ratio
    else:
        score = 0
        overlap_pmn = 0
        average_overlap_pmn = 0
        combined_score = 0


    if getdistance:
        distance = bgbp.get_gene_distance(*geneset)

    cstats = {}
    cstats['GeneSet'] = geneset
    if len(geneset) == 3:
        cstats['Gene0'], cstats['Gene1'], cstats['Gene2'] = tuple(geneset)
        cstats['MutationFrequency0'], cstats['MutationFrequency1'], cstats['MutationFrequency2'] = tuple(mutationfrequencies)
    else:
        cstats['Gene0'], cstats['Gene1'] = tuple(geneset)
        cstats['MutationFrequency0'], cstats['MutationFrequency1'] = tuple(mutationfrequencies)
    cstats['Type'] = 'Cooccurring'
    cstats['Probability'] = cooccurprob
    cstats['MutationFrequencies'] = mutationfrequencies
    cstats['MutationFrequencyDifference'] = abs(mutationfrequencies[0] - mutationfrequencies[1])
    cstats['MutationFrequencyDifferenceRatio'] = cstats['MutationFrequencyDifference'] * 1.0 / coverage
    cstats['Overlap'] = overlap
    cstats['CooccurrenceRatio'] = c_ratio
    cstats['Coverage'] = coverage
    cstats['SetScore'] = score
    cstats['AverageOverlapPMN'] = average_overlap_pmn
    #cstats['TotalGeneLength'] = total_gene_length
    cstats['CombinedScore'] = combined_score
    cstats['Concordance'] = cstats['Gene0'][-4:] == cstats['Gene1'][-4:]

    cstats['Somatic'] = 0
    if (cstats['Gene0'][-4:] not in {'loss', 'gain'}):
        cstats['Somatic'] += 1
    if (cstats['Gene1'][-4:] not in {'loss', 'gain'}):
        cstats['Somatic'] += 1


    cstats['RoundedLogPCov'] = round(-np.log10(cooccurprob/coverage), 1)


    if getdistance:
        cstats['Distance'] = distance

    return cstats
    #return mutationfrequencies, overlap, coverage, c_ratio, cooccurprob, score, overlap_pmn, average_overlap_pmn, combined_score


def writeheatmap(numCases, geneToCases, patientToGenes, file_prefix, heats=['CooccurrenceRatio']):
    """Given the geneToCases and patientToGenes dictionaries, write the heatmaps of all the scores in: heats. Write with
    rows and columns."""

    # Start with all possible gene sets.

    # start with np array
    # For each gene in genes, and each gene in genes
    # analyze the cooccur set. Get the entries

    genes = geneToCases.keys()
    numgenes = len(genes)
    heatmap_list = []
    #All defaults are 0s.
    for h in range(len(heats)):
        heatmap_list.append(np.zeros([numgenes, numgenes]))


    for g in range(numgenes):
        for g2 in range(numgenes):
            if g != g2:
                geneset = [genes[g], genes[g2]]
                mutationfrequencies, overlap, coverage, c_ratio, cooccurprob, \
                score, overlap_pmn, average_overlap_pmn, combined_score = analyze_cooccur_set(numCases, geneToCases, patientToGenes,
                                                                                      geneset)
                total_gene_length = sum([sco.gene_to_length(gene) for gene in geneset])
                entry = collections.OrderedDict()
                entry['GeneSet'] = geneset
                entry['Probability'] = cooccurprob
                entry['MutationFrequencies'] = mutationfrequencies
                entry['Overlap'] = overlap
                entry['CooccurrenceRatio'] = c_ratio
                entry['Coverage'] = coverage
                entry['SetScore'] = score
                entry['AverageOverlapPMN'] = average_overlap_pmn
                entry['TotalGeneLength'] = total_gene_length
                entry['CombinedScore'] = combined_score


                for h in range(len(heats)):
                    heatmap_list[h][g][g2] = entry[heats[h]]


    # Make the row and column headers


    # Write all the heats.

    for h in range(len(heats)):
        filename = file_prefix + '_' + heats[h] + '.tsv'
        heatmap = heatmap_list[h]

        with open(filename, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter = '\t')
            writer.writerow(genes)
            for g in range(len(genes)):
                writer.writerow(heatmap[g].tolist())

            # row_format ="{:>15}" * (len(teams_list) + 1)
            # print row_format.format("", *teams_list)
            # for team, row in zip(teams_list, data):
            #     print row_format.format(team, *row)








    # if (sum(mutationfrequencies) < numCases):
    #     if mep.issigmutex(numCases, mutationfrequencies, p, overlap):
    #         if aremutex(geneset, geneToCases, maxOverlap = maxOverlap):


#
#
#
#
#
#
#     score = sco.score_cooccur_set(geneset, geneToCases, patientToGenes)
#     overlap_pmns = [len(patientToGenes[patient]) for patient in overlaps(geneset, geneToCases)]
#     average_overlap_pmn = sum(overlap_pmns) * 1.0 / len(overlap_pmns)
#     # Cooccur Score
#     # Probability
#         # Can approximate, or not
#
#     # If mutex
#
#     # If cooccur







def runwholepipeline(args):
    """
      * **m** (*int*) - number of genes.
      * **n** (*int*) - number of patients.
      * **genes** (*list*) - genes in the mutation data.
      * **patients** (*list*) - patients in the mutation data.
      * **geneToCases** (*dictionary*) - mapping of genes to the patients in which they are mutated.
      * **patientToGenes** (*dictionary*) - mapping of patients to the genes they have mutated.
      """
    # Record Starting Time
    t = time.time()

    # Parse the arguments into shorter variable handles
    mutationmatrix = args.mutation_matrix
    output_prefix = args.output_prefix
    geneFile = args.gene_file
    patientFile = args.patient_file
    gene_blacklist = args.gene_blacklist_file
    patient_blacklist = args.patient_blacklist_file
    gene_file_1 = args.gene_list_1
    gene_file_2 = args.gene_list_2
    minFreq = args.min_freq
    set_size = args.set_size
    type = args.type
    maxOverlap = args.max_overlaps
    minCooccur = args.min_cooccur
    p_value = args.p_value
    cytoscape_output = args.cytoscape_output
    complete_output = args.complete_output
    top_percentile = args.top_percentile
    top_number = args.top_number
    just_sets = args.just_sets
    plot_distribution = args.plot_distribution
    plot_patient_distribution = args.plot_patient_distribution
    min_cooccurrence_ratio = args.min_cooccurrence_ratio




    #Load the matrix, and remove the blacklists
    mutations = remove_blacklists(gene_blacklist, patient_blacklist,
                                  *load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq))
    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mutations

    # Remove verbose output for convenience when getting summary stats
    if (type != 'sss'):
        print 'Pre-filter Mutation data: %s genes x %s patients' % (numGenes, numCases)



    # Graph distributions first.
    if plot_distribution:
        graph_mutation_distribution(numCases, genes, geneToCases, mutationmatrix, bins=plot_distribution)
    if plot_patient_distribution:
        graph_patient_distribution(numGenes, patients, patientToGenes, mutationmatrix, bins=plot_patient_distribution)


    # Filter by the top percentile or top number.
    numGenes, numCases, genes, patients, geneToCases, patientToGenes = filtermatrix(top_percentile, *mutations)
    if top_number:
        numGenes, numCases, genes, patients, geneToCases, patientToGenes = filtertopnumbermatrix(top_number, *mutations)

    if (type != 'sss'):
        print 'Post-filter Mutation data: %s genes x %s patients' % (numGenes, numCases)
    # print "Average number of mutations per patient: ", sum([len(patientToGenes[patient]) for patient in patients]) * 1.0 / numCases
    # print "Standard deviation of mutations per patient: ", np.std([len(patientToGenes[patient]) for patient in patients])
    # print "Median of mutations per patient:", np.median([len(patientToGenes[patient]) for patient in patients])




    if gene_file_1:
        with open(gene_file_1) as f:
            gene_list_1 = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
    else:
        gene_list_1 = None

    if gene_file_2:
        with open(gene_file_2) as f:
            gene_list_2 = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
    else:
        gene_list_2 = None



    if (type == 'm'):

        msigsets, mgene_occurrences = getmutexsets(numCases, genes, geneToCases, patientToGenes, set_size, p_value,
                                                   maxOverlap, genes1=gene_list_1, genes2=gene_list_2)

        if (output_prefix == None):
            file = sys.stdout
            writemutexsets(file, msigsets, mgene_occurrences, p_value)
        else:
            file_prefix = output_prefix + '.g' + str(numGenes) + '.p' + str(numCases) + '.s' + str(set_size) + \
                          '.ns' + str(len(msigsets)) + '.ng' + str(len(mgene_occurrences)) + '.mf' + \
                          str(minFreq) + '.p' + str(p_value) + '.mo' + str(maxOverlap)
            file = open(file_prefix + '.tsv', 'w')
            writemutexsets(file, msigsets, mgene_occurrences, p_value)
            file.close()

            if cytoscape_output: writemutex_cytoscape(file_prefix, msigsets, mgene_occurrences, p_value)

            if complete_output:
                filename = file_prefix + '_complete' + '.tsv'
                writecompleteinfo(filename, msigsets, mgene_occurrences)

        print 'Number of sets: ', len(msigsets)
        print "Number of genes: ", len(mgene_occurrences)
        print 'Time used: ', time.time() - t

    elif (type == 'c'):

        csigsets, cgene_occurrences = getcooccursets_new(numCases, genes, geneToCases, patientToGenes,
                                                         set_size, p_value, minCooccur,
                                                         just_sets=just_sets,
                                                         min_cooccurrence_ratio=min_cooccurrence_ratio,
                                                         genes1=gene_list_1, genes2=gene_list_2)

        if (output_prefix == None):
            file = sys.stdout
            writecooccursets(file, csigsets, cgene_occurrences, p_value)
        else:
            file_prefix = output_prefix + '.g' + str(numGenes) + '.p' + str(numCases) + '.s' + str(set_size) \
                          + '.ns' + str(len(csigsets)) + '.ng' + str(len(cgene_occurrences)) + '.mf' \
                          + str(minFreq) + '.tn' + str(top_number) + '.p' + str(p_value) + '.mc' + str(minCooccur) + \
                          '.mcr' + str(min_cooccurrence_ratio)

            file = open(file_prefix + '.tsv', 'w')
            writecooccursets(file, csigsets, cgene_occurrences, p_value)
            file.close()

            graph_cooccurrenceratio_distribution(csigsets, mutationmatrix)

            if cytoscape_output: writecooccur_cytoscape(file_prefix, csigsets, cgene_occurrences, p_value)

            if complete_output:
                filename = file_prefix + '_complete' + '.tsv'
                writecompleteinfo(filename, csigsets, cgene_occurrences)

        print 'Number of sets: ', len(csigsets)
        print "Number of genes: ", len(cgene_occurrences)
        print 'Time used: ', time.time() - t

    elif (type == 'h'):
        file_prefix = output_prefix + '.g' + str(numGenes) + '.p' + str(numCases)
        writeheatmap(numCases, geneToCases, patientToGenes, file_prefix,
                     heats=['CooccurrenceRatio', 'Probability', 'SetScore'])
        print 'Time used: ', time.time() - t

    elif (type == 'sss'):
        print "Cancer\tPatient Number\tAverageMutations\tStandardDeviationMutations\tMinimum\tFirstQuartile\tMedianMutations\tThirdQuartile\tMaximum"
        print mutationmatrix, "\t", numCases, "\t", sum([len(patientToGenes[patient]) for patient in patients]) * 1.0 / numCases, \
             "\t", np.std([len(patientToGenes[patient]) for patient in patients]), "\t", np.min([len(patientToGenes[patient]) for patient in patients]), \
            "\t", np.percentile([len(patientToGenes[patient]) for patient in patients], 25), \
            "\t", np.median([len(patientToGenes[patient]) for patient in patients]), "\t", np.percentile([len(patientToGenes[patient]) for patient in patients], 75), \
            "\t", np.percentile([len(patientToGenes[patient]) for patient in patients], 100)


def get_parser():
    # Parse arguments
    import argparse

    description = 'Given an input number of samples n, and probability of mutual exclusivity p, ' \
                  'plots the number of mutations that each sample' \
                  'must have in order to reach that probability p.'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-o', '--output_prefix', default=None,
                        help='Output path prefix (TSV format). Otherwise just prints.')
    parser.add_argument('-mf', '--min_freq', type=int, default=0,
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-tp', '--top_percentile', type=float, default=100,
                        help='Limit to this percentage of mutations of greatest abundance.')
    parser.add_argument('-tn', '--top_number', type=int, default=0,
                        help='Limit to this number of mutations of greatest abundance.')

    parser.add_argument('-pf', '--patient_file', default=None,
                        help='File of patients to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None,
                        help='File of genes to be included (optional).')
    parser.add_argument('-pbf', '--patient_blacklist_file', default=None,
                        help='File of patients to be excluded (optional).')
    parser.add_argument('-gbf', '--gene_blacklist_file', default=None,
                        help='File of genes to be excluded (optional).')


    parser.add_argument('-gl1', '--gene_list_1', default=None,
                        help='First sets of genes to draw from')

    parser.add_argument('-gl2', '--gene_list_2', default=None,
                        help='Second set of genes to draw from')

    parser.add_argument('-s', '--set_size', type=int, default=2, help='Number of genes per set')
    parser.add_argument('-t', '--type', default='m',
                        help='Use m for mutual exclusivity, c for cooccuring')
    parser.add_argument('-mo', '--max_overlaps', type=int, default=10,
                        help='Maximum allowed number of overlapping mutations per mutually exclusive set')
    parser.add_argument('-mc', '--min_cooccur', type=int, default=1,
                        help='Minimum number of cooccurrences per cooccuring set')
    parser.add_argument('-p', '--p_value', type=float, default=0.05,
                        help='Significance Threshold.')


    parser.add_argument('-v', '--cytoscape_output', type=int, default=1, help='Set to 1 to produce cytoscape-compatible'
                                                                              'outputs')
    parser.add_argument('-co', '--complete_output', type=int, default=1, help='Set to 1 to produce a file with complete information')

    #parser.add_argument('-cif', '--check_in_file', default=None, help='Name of file to check if the set is in')


    parser.add_argument('-js', '--just_sets', type=int, default=0, help='Just find the possible cooccurring sets,'
                                                                        'without checking for significance')
    parser.add_argument('-pd', '--plot_distribution', type=int, default=0,
                        help='Plot distribution of mutation frequencies '
                             'with this number of bins before running')
    parser.add_argument('-ppd', '--plot_patient_distribution', type=int, default=0,
                        help='Plot distribution of mutation number'
                             'per gene '
                             'with this number of bins before running')

    parser.add_argument('-mcr', '--min_cooccurrence_ratio', type=float, default=0.0,
                        help='The minimum cooccurrence ratio for a set to be deemed significant.')


    return parser


def main():
    runwholepipeline(get_parser().parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()
