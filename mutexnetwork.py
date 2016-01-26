__author__ = 'jlu96'

import mutex

from comet import permute
import time
import mutex as mex
import mutex_triangles as met
import csv
import random
import os
import parallel_compute_working as pac
import numpy as np
import math
import line_profiler

class TwoWayDict(dict):
    def __setitem__(self, key, value):
        # Remove any previous connections with these values
        if key in self:
            del self[key]
        if value in self:
            del self[value]
        dict.__setitem__(self, key, value)
        dict.__setitem__(self, value, key)

    def __delitem__(self, key):
        dict.__delitem__(self, self[key])
        dict.__delitem__(self, key)

    def __len__(self):
        """Returns the number of connections"""
        return dict.__len__(self) // 2

class PermutationMatrix:
    def __init__(self, geneToCases, patientToGenes):
        # Each is (gene, ID, sum)
        self.geneToCases = geneToCases
        self.patientToGenes = patientToGenes
        self.genes = geneToCases.keys()
        self.patients = patientToGenes.keys()
        self.numGenes = len(geneToCases)
        self.numCases = len(patientToGenes)
        self.gene_ID = TwoWayDict()
        self.patient_ID = TwoWayDict()
        self.gene_ID_sums = {}
        self.patient_ID_sums = {}

        # Each entry is (gene_ID, patient_ID)
        self.ones = set()
        self.zeros = set()

        ID = 0
        for gene in geneToCases:
            self.gene_ID[gene] = ID
            sum = len(geneToCases[gene])
            self.gene_ID_sums[ID] = sum
            ID += 1

        ID = 0
        for patient in patientToGenes:
            self.patient_ID[patient] = ID
            sum = len(patientToGenes[patient])
            self.patient_ID_sums[ID] = sum
            ID += 1



    #@profile
    def permute_data(self):
        # Take gene IDs

        binary_numbers = self.random_binary_numbers()

        conformed_binary_numbers = self.conform_binary_numbers(binary_numbers)

        self.update_gene_patient_dicts(conformed_binary_numbers)




    # FUNCTIONS TO GENERATE THE RANDOM MATRIX---------------------------------------------------------------------------

    ##@profile
    def random_binary_matrix_ones_zeros(self):
        new_zeros = set()
        new_ones = set()


        # Generate random binary matrix for each gene ID
        for g in self.gene_ID_sums:
            sum = self.gene_ID_sums[g]
            binary_array = random_binary_array(sum, self.numCases)
            for p in range(len(binary_array)):
                if binary_array[p]:
                    new_ones.add((g,p))
                else:
                    new_zeros.add((g,p))
        return new_ones, new_zeros

    ##@profile
    def random_binary_matrix(self):
        binary_matrix = []
        for g in self.gene_ID_sums:
            sum = self.gene_ID_sums[g]
            binary_array = random_binary_array(sum, self.numCases)
            binary_matrix.append(binary_array)
        return binary_matrix

    #@profile
    def random_binary_numbers(self):
        binary_numbers = []
        for g in self.gene_ID_sums:
            g_sum = self.gene_ID_sums[g]
            binary_array = random_binary_array(g_sum, self.numCases)
            binary_number = 0
            for i in binary_array:
                binary_number += i
                binary_number <<= 1
            binary_number >>= 1
            binary_numbers.append(binary_number)
        return binary_numbers



    # FUNCTIONS TO CONFORM THE RANDOM MATRIX----------------------------------------------------------------------------

    #@profile
    def conform_binary_numbers(self, binary_numbers):

        for p in range(self.numCases):
            conform_sum = self.patient_ID_sums[p]
            p_mask = 1 << (self.numCases - p - 1)
            p_sum = 0
            p_ones = []
            p_zeros = []

            for g in range(self.numGenes):
                # print bin(binary_numbers[g])
                # print bin(p_mask)

                if (binary_numbers[g] & p_mask):
                    p_sum += 1
                    p_ones.append(g)
                else:
                    p_zeros.append(g)

            ones_to_add = conform_sum - len(p_ones)
            if ones_to_add > 0:

                # Shuffle the rows to search
                random.shuffle(p_zeros)

                # Make a new p_mask to look at later columns
                latter_mask = p_mask
                # print
                for i in range(self.numCases - p):
                    latter_mask >>= 1
                    # print "####################################################################################"
                    # print "BEGINNING!!!!-----------------------"
                    for g in p_zeros:
                        # If the later column has a 1, switch that column to 0 and p_column to 1
                        if (binary_numbers[g] & latter_mask):

                            binary_numbers[g] -= latter_mask
                            binary_numbers[g] += p_mask
                            p_zeros.remove(g)

                            # print "SWITCH HAPPENING FOR # ", g
                            # orig = binary_numbers[g] - p_mask + latter_mask
                            # p_array = binary_number_to_array(p_mask, self.numCases)
                            # latter_array = binary_number_to_array(latter_mask, self.numCases)
                            # orig_array = binary_number_to_array(orig, self.numCases)
                            # new_array = binary_number_to_array(binary_numbers[g], self.numCases)
                            #
                            # for index in range(self.numCases):
                            #     if p_array[index]:
                            #         print "p_array ", index
                            #     if latter_array[index]:
                            #         print "latter_array ", index
                            #     if orig_array[index]:
                            #         print "orig_array ", index
                            #     if new_array[index]:
                            #         print "new_array ", index


                            #
                            if sum(binary_number_to_array(binary_numbers[g], self.numCases)) != self.gene_ID_sums[g]:
                                print "in ones to add"
                                print "Sum violated for column ", i
                                # print bin(p_mask), bin(latter_mask), bin(binary_numbers[g]), self.gene_ID_sums[g]
                                # print p_mask
                                # orig = binary_numbers[g] - p_mask + latter_mask
                                # p_array = binary_number_to_array(p_mask, self.numCases)
                                # latter_array = binary_number_to_array(latter_mask, self.numCases)
                                # orig_array = binary_number_to_array(orig, self.numCases)
                                # new_array = binary_number_to_array(binary_numbers[g], self.numCases)
                                #
                                # for i in range(self.numCases):
                                #     if p_array[i]:
                                #         print "p_array ", i
                                #     if latter_array[i]:
                                #         print "latter_array ", i
                                #     if orig_array[i]:
                                #         print "orig_array ", i
                                #     if new_array[i]:
                                #         print "new_array ", i

                                # p_digit = math.log(p_mask, 2)
                                #
                                # latter_dig = math.log(latter_mask, 2)
                                # orig_dig = math.log(orig, 2)
                                # now_dig = math.log(binary_numbers[g], 2)
                                # print "p_digit", p_digit
                                # print "latter_dig", latter_dig
                                # print "orig_dig", orig_dig
                                # print "now_dig", now_dig

                                exit(1)

                            ones_to_add -= 1
                            if ones_to_add == 0:
                                break
                    if ones_to_add == 0:
                        break

                if ones_to_add != 0:
                    print "Error: ", ones_to_add, " unfinished added ones"
                    exit(1)

            elif ones_to_add < 0:
                zeros_to_add = -ones_to_add

                # Shuffle the rows to search
                random.shuffle(p_ones)

                # Make a new mask to look at later columns
                latter_mask = p_mask

                for i in range(self.numCases - p):
                    latter_mask >>= 1
                    for g in p_ones:
                        # If the later column has a 0, switch that column to 1 and p_column to 0
                        if not (binary_numbers[g] & latter_mask):
                            binary_numbers[g] += latter_mask
                            binary_numbers[g] -= p_mask

                            # This one can no longer be switched
                            p_ones.remove(g)

                            if sum(binary_number_to_array(binary_numbers[g], self.numCases)) != self.gene_ID_sums[g]:
                                print "Sum violated for column ", i
                                print bin(p_mask), bin(latter_mask), bin(binary_numbers[g]), self.gene_ID_sums[g]
                                exit(1)
                            zeros_to_add -= 1
                            if zeros_to_add == 0:
                                break

                    if zeros_to_add == 0:
                        break

                if zeros_to_add != 0:
                    print "Error: ", zeros_to_add, " unfinished added zeros"
                    exit(1)


        return binary_numbers

    #@profile
    def update_gene_patient_dicts(self, binary_numbers):
        self.geneToCases = {}
        self.patientToGenes = {}

        for g in range(len(binary_numbers)):
            binary_array = binary_number_to_array(binary_numbers[g], self.numCases)
            for p in range(len(binary_array)):
                if binary_array[p]:
                    gene = self.gene_ID[g]
                    patient = self.patient_ID[p]

                    if gene not in self.geneToCases:
                        self.geneToCases[gene] = set()
                    if patient not in self.patientToGenes:
                        self.patientToGenes[patient] = set()

                    self.geneToCases[gene].add(patient)
                    self.patientToGenes[patient].add(gene)

    def getGenePatientDicts(self):
        return self.geneToCases, self.patientToGenes

    def test_conformity(self):
        """
        Test if the current dictionaries satisfy the gene row sums and patient column sums
        """
        conform_gene_sums = [self.gene_ID_sums[g] for g in range(self.numGenes)]
        conform_patient_sums = [self.patient_ID_sums[p] for p in range(self.numCases)]


        for g in range(self.numGenes):
            gene = self.gene_ID[g]
            g_sum = len(self.geneToCases[gene])
            if g_sum != conform_gene_sums[g]:
                print "Sum for g ", g, " is ", g_sum, ", not ", conform_gene_sums[g]
                return False

        for p in range(self.numCases):
            patient = self.patient_ID[p]
            p_sum = len(self.patientToGenes[patient])
            if p_sum != conform_patient_sums[p]:
                print "Sum for patient ", self.patient_ID[p], " is ", p_sum, ", not ", conform_patient_sums[p]
                return False

        return True





    ##@profile
    def conform_binary_matrix_ones_zeros(self, ones, zeros):
        """
        Make the binary matrix fit the given patient column sums.
        """

        decided_ones, decided_zeros = set(), set()

        for p in range(self.numCases):

            # Get the ones and zeros of these patients
            p_ones = set([entry for entry in ones if entry[1] == p])
            ones = ones.difference(p_ones)
            p_zeros = set([entry for entry in zeros if entry[1] == p])
            zeros = zeros.difference(p_zeros)

            ones_to_add = self.patient_ID_sums[p] - len(p_ones)

            # If ones are missing from the patient column, switch some number of the zeros with a one, and
            # correspondingly switch some number of ones with zeros in the other

            if ones_to_add > 0:
                new_p_ones = random.sample(p_zeros, ones_to_add)
                new_zeros = random.sample(ones, ones_to_add)

                p_ones = p_ones.union(new_p_ones)
                p_zeros = p_zeros.difference(new_p_ones)
                ones = ones.difference(new_zeros)
                zeros = zeros.union(new_zeros)


            elif ones_to_add < 0:
                zeros_to_add = -ones_to_add
                new_p_zeros = random.sample(p_ones, zeros_to_add)
                new_ones = random.sample(zeros, zeros_to_add)

                p_ones = p_ones.difference(new_p_zeros)
                p_zeros = p_zeros.union(new_p_zeros)
                ones = ones.union(new_ones)
                zeros = zeros.difference(new_ones)

            else:
                pass

            decided_ones = decided_ones.union(p_ones)
            decided_zeros = decided_zeros.union(p_zeros)

        return decided_ones, decided_zeros





    def getGenePatientDicts_ones_zeros(self):
        geneToCases = {}
        patientToGenes = {}

        for g,p in self.ones:
            gene = self.gene_ID[g]
            patient = self.patient_ID[p]

            if gene not in geneToCases:
                geneToCases[gene] = set()
            geneToCases[gene].add(patient)

            if patient not in patientToGenes:
                patientToGenes[patient] = set()
            patientToGenes[patient].add(gene)

        return geneToCases, patientToGenes




    ##@profile
    def test_conformity_ones_zeros(self):
        """
        Test if the current matrix satisfy the gene row sums and patient column sums
        """
        conform_gene_sums = [self.gene_ID_sums[g] for g in range(self.numGenes)]
        conform_patient_sums = [self.patient_ID_sums[p] for p in range(self.numCases)]


        for g in range(self.numGenes):
            g_sum = sum([1 for entry in self.ones if entry[0] == g])
            if g_sum != conform_gene_sums[g]:
                print "Sum for gene ", self.gene_ID[g], " is ", g_sum, ", not ", conform_gene_sums[g]
                return False

        for p in range(self.numCases):
            p_sum = sum([1 for entry in self.ones if entry[1] == p])
            if p_sum != conform_patient_sums[p]:
                print "Sum for patient ", self.patient_ID[p], " is ", p_sum, ", not ", conform_patient_sums[p]
                return False

        return True



def random_binary_array(num_ones, length):

    random_row = [1] * num_ones + [0] * (length - num_ones)
    random.shuffle(random_row)

    return random_row


def binary_array_to_number(array):
    string = "".join([str(i) for i in array])
    return int(string, 2)


def binary_number_to_array(number, length):
    bin_string = str(bin(number))[2:]
    leading_zeros = [0] * (length - len(bin_string))

    return leading_zeros + [int(i) for i in bin_string]








# Matrix class that takes input matrix, holds all relevant matrices

class PermutationMatrices:
    # generate matrices
    #@profile
    def __init__(self, geneToCases, patientToGenes, num_permutations, seeds=[], Q=100, matrixdirectory=None, binary_perm_method=False,
                 write_matrices=False, load_matrices=False):
        t_start = time.time()

        if not seeds:
            for i in range(num_permutations):
                seeds.append(random.random())

        self.seeds = seeds


        genes = geneToCases.keys()
        patients = patientToGenes.keys()
        self.num_permutations = num_permutations


        # The original matrix
        self.geneToCases_orig = geneToCases
        self.patientToGenes_orig = patientToGenes
        self.numCases = len(self.geneToCases_orig)
        self.numGenes = len(self.patientToGenes_orig)

        # A dictionary that holds each gene to a list of sets of cases. Each entry in the list is the gene's cases
        # for one of the permutations
        self.geneToCases_perm = {}
        for gene in genes:
            self.geneToCases_perm[gene] = []

        # Same as above, but from patient to genes.
        self.patientToGenes_perm = {}
        for patient in patients:
            self.patientToGenes_perm[patient] = []



        # Generate output directory to write matrix
        if not os.path.exists(os.path.dirname(matrixdirectory)) and write_matrices:
            os.makedirs(os.path.dirname(matrixdirectory))



        if binary_perm_method:
            PM = PermutationMatrix(geneToCases, patientToGenes)
            print "Using matrix permutation method"
        else:
            G = permute.construct_mutation_graph(geneToCases, patientToGenes)
            print '\t- Graph has', len( G.edges() ), 'edges among', len( G.nodes() ), 'nodes.'


        for i in range(num_permutations):

            t = time.time()

            if binary_perm_method:
                PM.permute_data()
                if not PM.test_conformity():
                    print "Permuted matrix failed test"
                    exit(1)
                newGeneToCases, newPatientToGenes = PM.getGenePatientDicts()
                print i + 1, " matrices generated in ", time.time() - t
                print " number of same alterations is ", sum([len(patientToGenes[patient].intersection(newPatientToGenes[patient])) for patient in patients])
            # Matrix method---------------------------------------------------------------------
            else:
                _, _, _, _, newGeneToCases, newPatientToGenes = permute.permute_mutation_data(G, genes, patients, self.seeds[i], Q)
                print i + 1, " matrices generated in ", time.time() - t
                print "For Q=", Q, " number of same alterations is ", sum([len(patientToGenes[patient].intersection(newPatientToGenes[patient])) for patient in patients])

            # Old graph method------------------------------------------------------------------
                       # Old graph method------------------------------------------------------------------

            for gene in newGeneToCases:
                self.geneToCases_perm[gene].append(newGeneToCases[gene])

            for patient in newPatientToGenes:
                self.patientToGenes_perm[patient].append(newPatientToGenes[patient])


            # Write permuted matrices out
            if write_matrices and matrixdirectory:
                adj_list = [ p + "\t" + "\t".join( sorted(newPatientToGenes[p]) ) for p in patients ]

                permutation_file = "{}/permuted-matrix-{}.m2".format(matrixdirectory, i+1)
                with open(permutation_file, 'w') as outfile: outfile.write('\n'.join(adj_list))

            # Clear dictionaries from stack once they're not used
            newGeneToCases.clear()
            newPatientToGenes.clear()

        print "Time to generate mutation matrices ", time.time() - t_start

    # # Method to iterate over its matrices. Return value: number of satisfying matrices
    # def pvalue(self, condition_function):
    #
    #     num_satisfy = 0
    #
    #     for i in range(self.num_permutations):
    #         geneToCases = self.geneToCases_perm.copy()
    #         for gene in geneToCases:
    #             geneToCases[gene] = self.geneToCases_perm[gene][i]
    #
    #
    #         if condition_function.apply(geneToCases):
    #             num_satisfy += 1
    #
    #     return num_satisfy * 1.0 / self.num_permutations

    def set_to_pvalue(self, set_condition_function_list):

        set_to_pvalue = {}
        for set, condition_function in set_condition_function_list:
            set_to_pvalue[set] = 0


        for i in range(self.num_permutations):
            geneToCases = {}
            for gene in self.geneToCases_perm:
                geneToCases[gene] = self.geneToCases_perm[gene][i]

            for set, condition_function in set_condition_function_list:
                if condition_function.apply(geneToCases):
                    set_to_pvalue[set] += 1

        for set in set_to_pvalue:
            set_to_pvalue[set] /= (1.0 * self.num_permutations)

        return set_to_pvalue

    def complete_cooccurpairs(self, genepairs, p=0.05, minCooccur=1, min_cooccurrence_ratio=0.0, parallel_compute_number=0,
                              compute_scores=True):
        """
        :param genepairs:
        :param cprob:
        :param minCooccur:
        :param min_cooccurrence_ratio:
        :param parallel_compute_number:
        :return: cpairsdict, cgenedict
        """


        print "Generating list of", len(genepairs), " co-occurring hypotheses to test on permutation matrices..."
        # Generate condition functions after analyzing each gene pair for Co-occurrence/min
        cooccur_set_condition_function_list = []

        # Generate list of condition functions to test the permutation matrix for
        for genepair in genepairs:
            ConditionFunction = Condition(None)

            condition_dict = {}
            condition_dict['Genes'] = tuple(genepair)
            condition_dict['Overlap'] = len(set.intersection(*[self.geneToCases_orig[gene] for gene in condition_dict['Genes']]))
            condition_dict['Mutex'] = False

            ConditionFunction.set_params([condition_dict])

            cooccur_set_condition_function_list.append((genepair, ConditionFunction))

        print "Done. Now, calulating p-values of hypotheses..."

        # Generate co-occurring pairs
        if parallel_compute_number:
            cooccur_pair_to_pvalue = pac.parallel_compute_new(self.set_to_pvalue, [cooccur_set_condition_function_list],
                                                         cooccur_set_condition_function_list, 0, pac.partition_inputs, {0: pac.combine_dictionaries},
                                                         number=parallel_compute_number,
                                                         procnumber=parallel_compute_number)
        else:
            cooccur_pair_to_pvalue = self.set_to_pvalue(cooccur_set_condition_function_list)


        print "Done. Now, finding co-occurring pairs"
        # Generate dictionary for each pair. Optionally analyze each cooccur set as well.
        cpairsdict = {}
        cgenedict = {}

        for genepair in cooccur_pair_to_pvalue:
            if cooccur_pair_to_pvalue[genepair] < p:

                cstats = mex.analyze_cooccur_set_new(self.numCases, self.geneToCases_orig, self.patientToGenes_orig,
                                                     geneset=tuple(genepair), compute_scores=compute_scores)


                if cstats['Overlap'] >= minCooccur and cstats['CooccurrenceRatio'] >= min_cooccurrence_ratio:


                    cstats['PermutationProbability'] = cooccur_pair_to_pvalue[genepair]
                    cpairsdict[genepair] = cstats
                    gene1, gene2 = tuple(genepair)
                    if gene1 not in cgenedict:
                        cgenedict[gene1] = set()
                        cgenedict[gene1].add(gene2)
                    else:
                        cgenedict[gene1].add(gene2)

                    if gene2 not in cgenedict:
                        cgenedict[gene2] = set()
                        cgenedict[gene2].add(gene1)
                    else:
                        cgenedict[gene2].add(gene1)

        return cpairsdict, cgenedict


    def complete_mutexpairs(self, genepairs, p=0.05, maxOverlap=200, parallel_compute_number=0):


        print "Generating list of", len(genepairs), " mutually exclusive hypotheses to test on permutation matrices..."
        # Generate condition functions after analyzing each gene pair for Co-occurrence/min
        mutex_set_condition_function_list = []

        # Generate list of condition functions to test the permutation matrix for
        for genepair in genepairs:
            ConditionFunction = Condition(None)

            condition_dict = {}
            condition_dict['Genes'] = tuple(genepair)
            condition_dict['Overlap'] = len(set.intersection(*[self.geneToCases_orig[gene] for gene in condition_dict['Genes']]))
            condition_dict['Mutex'] = True

            ConditionFunction.set_params([condition_dict])

            mutex_set_condition_function_list.append((genepair, ConditionFunction))

        print "Done. Now, calulating p-values of hypotheses..."

        # Generate co-occurring pairs
        if parallel_compute_number:
            mutex_pair_to_pvalue = pac.parallel_compute_new(self.set_to_pvalue, [mutex_set_condition_function_list],
                                                         mutex_set_condition_function_list, 0, pac.partition_inputs, {0: pac.combine_dictionaries},
                                                         number=parallel_compute_number,
                                                         procnumber=parallel_compute_number)
        else:
            mutex_pair_to_pvalue = self.set_to_pvalue(mutex_set_condition_function_list)


        print "Done. Now, finding mutually exclusive pairs"
        # Generate dictionary for each pair. Optionally analyze each mutex set as well.
        mpairsdict = {}
        mgenedict = {}

        for genepair in mutex_pair_to_pvalue:
            if mutex_pair_to_pvalue[genepair] < p:

                mstats = mex.analyze_mutex_set_new(self.numCases, self.geneToCases_orig, self.patientToGenes_orig,
                                                     geneset=tuple(genepair))

                if mstats['Overlap'] <= maxOverlap:

                    mstats['PermutationProbability'] = mutex_pair_to_pvalue[genepair]

                    mpairsdict[genepair] = mstats
                    gene1, gene2 = tuple(genepair)
                    if gene1 not in mgenedict:
                        mgenedict[gene1] = set()
                        mgenedict[gene1].add(gene2)
                    else:
                        mgenedict[gene1].add(gene2)

                    if gene2 not in mgenedict:
                        mgenedict[gene2] = set()
                        mgenedict[gene2].add(gene1)
                    else:
                        mgenedict[gene2].add(gene1)


        return mpairsdict, mgenedict

# A class of test functions for each pair, to asses for mutual exclusivity/co-occurrence. Used for the permutation
# matrices.
class Condition:
    def __init__(self, condition_dicts):
        self.set_params(condition_dicts)

    def set_params(self, condition_dicts):
        self.conditions = condition_dicts

    def apply(self, geneToCases):
        """
        :param geneToCases: Input matrix
        :return: Whether the input matrix satisfies all the conditions.
        """

        satisfies = True

        for condition in self.conditions:
            genes = condition['Genes']
            overlap = condition['Overlap']
            mutex = condition['Mutex']

            found_overlap = len(set.intersection(*[geneToCases[gene] for gene in genes]))

            satisfies = satisfies and (found_overlap <= overlap if mutex else found_overlap >= overlap)
            # print "Genes are ", genes, " found overlap is ", found_overlap, " overlap is ", overlap
            # if found_overlap != overlap:
            #     print "Genes are ", genes, " found overlap is ", found_overlap, " overlap is ", overlap

        return satisfies


def main():

    mutationmatrix = '/Users/jlu96/maf/new/PRAD_broad/PRAD_broad-som.m2'
    patientFile = None
    geneFile = None
    #'/Users/jlu96/conte/jlu/REQUIREDFILES_OnlyLoss2/COSMICGenes_OnlyLoss.txt'
    minFreq = 0
    num_permutations = 15
    binary_perm_method = False
    Q = 100
    write_matrices = True
    matrixdirectory = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/PRAD_broad-som-jl-' + ('matrix' if binary_perm_method else 'network')
    outmutexfile = matrixdirectory + '/mutex' + str(num_permutations) + str(time.time()) + '.tsv'
    outcooccurfile = matrixdirectory + '/cooccur' + str(num_permutations)  + str(time.time()) + '.tsv'
    outseedsfile = matrixdirectory + '/seeds' + str(time.time()) + '.tsv'
    if not os.path.exists(os.path.dirname(matrixdirectory)):
        os.makedirs(os.path.dirname(matrixdirectory))


    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)

    print "numGenes ", numGenes, " and numCases ", numCases

    for patient in patients:
        if not patientToGenes[patient]:
            patientToGenes.pop(patient)
            print patient, "popped"

    # Generate Permutation Matrices
    pm = PermutationMatrices(geneToCases, patientToGenes, num_permutations, Q=Q, matrixdirectory=matrixdirectory,
                             binary_perm_method=binary_perm_method, write_matrices=write_matrices)

    # Make list of pairs from highly mutated genes
    test_genes = [gene for gene in genes if len(geneToCases[gene]) > 5]
    # for test_gene in test_genes:
    #     print test_gene
    genepairs = met.getgenepairs(geneToCases, test_genes)
    print "Number of pairs to test ", len(genepairs)





    # CALCULATE MUTEX

    # Create a list of ConditionFunctions that you must later initialize...
    ConditionFunctions = range(len(genepairs))
    mutex_set_condition_function_list = []

    # Generate set_condition_function_list
    for i in range(len(genepairs)):
        genepair = genepairs[i]

        condition_dict = {}
        condition_dict['Genes'] = tuple(genepair)
        condition_dict['Overlap'] = len(set.intersection(*[geneToCases[gene] for gene in condition_dict['Genes']]))
        condition_dict['Mutex'] = True

        ConditionFunctions[i] = Condition([condition_dict])

        if [condition_dict] != ConditionFunctions[i].conditions:
            print condition_dict, ConditionFunctions[i].conditions


        mutex_set_condition_function_list.append((genepair, ConditionFunctions[i]))

    print "Finished mutex condition function list"

    t= time.time()
    # Calculate pvalues for mutual exclusivity
    pair_to_mutex = {}

    pair_to_mutex_network_pvalue = pm.set_to_pvalue(mutex_set_condition_function_list)
    print "mutex pair network pvalues finished in ", time.time() - t

    for genepair in genepairs:
        pair_to_mutex[genepair] = mex.analyze_mutex_set_new(numCases, geneToCases, patientToGenes, genepair)
        pair_to_mutex[genepair]['NetworkProbability'] = pair_to_mutex_network_pvalue[genepair]




    # Write to output
    with open(outmutexfile, 'w') as csvfile:
        fieldnames = pair_to_mutex[genepairs[0]].keys()
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for genepair in pair_to_mutex:
            writer.writerow(pair_to_mutex[genepair])



    # CALCULATE COOCCUR

    cooccur_set_condition_function_list = []

    # Generate set_condition_function_list
    for genepair in genepairs:
        ConditionFunction = Condition(None)

        condition_dict = {}
        condition_dict['Genes'] = tuple(genepair)
        condition_dict['Overlap'] = len(set.intersection(*[geneToCases[gene] for gene in condition_dict['Genes']]))
        condition_dict['Mutex'] = False

        ConditionFunction.set_params([condition_dict])

        cooccur_set_condition_function_list.append((genepair, ConditionFunction))



    t= time.time()
    # Calculate pvalues for mutual exclusivity
    pair_to_cooccur = {}

    pair_to_cooccur_network_pvalue = pm.set_to_pvalue(cooccur_set_condition_function_list)
    print "cooccur pair network pvalues finished in ", time.time() - t

    for genepair in genepairs:
        pair_to_cooccur[genepair] = mex.analyze_cooccur_set_new(numCases, geneToCases, patientToGenes, genepair)
        pair_to_cooccur[genepair]['NetworkProbability'] = pair_to_cooccur_network_pvalue[genepair]




    # Write to output
    with open(outcooccurfile, 'w') as csvfile:
        fieldnames = pair_to_cooccur[genepairs[0]].keys()
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for genepair in pair_to_cooccur:
            writer.writerow(pair_to_cooccur[genepair])


    # Write seeds to output
    with open(outseedsfile, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for seed in pm.seeds:
            writer.writerow([seed])


if __name__ == '__main__':
    main()


# and check for certain gene overlaps (certain conditions. Return list of tuples?)

# Method to return pvalue of mutexprob. Input: genes and overlap

# Ultimately want this to work with sample-specific permutations...?