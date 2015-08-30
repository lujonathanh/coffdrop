__author__ = 'jlu96'

import mutex

from comet import permute
import time
import mutex as mex
import mutex_triangles as met
import csv
import random
# Matrix class that takes input matrix, holds all relevant matrices

class PermutationMatrices:
    # generate matrices
    def __init__(self, geneToCases, patientToGenes, num_permutations, seeds=[], Q=100, matrixdirectory=None):
        t_start = time.time()

        if not seeds:
            for i in range(num_permutations):
                seeds.append(random.random())

        self.seeds = seeds


        genes = geneToCases.keys()
        patients = patientToGenes.keys()
        self.num_permutations = num_permutations

        # A dictionary that holds each gene to a list of sets of cases. Each entry in the list is the gene's cases
        # for one of the permutations
        self.geneToCases_perm = {}
        for gene in genes:
            self.geneToCases_perm[gene] = []

        # Same as above, but from patient to genes.
        self.patientToGenes_perm = {}
        for patient in patients:
            self.patientToGenes_perm[patient] = []



        G = permute.construct_mutation_graph(geneToCases, patientToGenes)
        print '\t- Graph has', len( G.edges() ), 'edges among', len( G.nodes() ), 'nodes.'

        for i in range(num_permutations):

            _, _, _, _, geneToCases, patientToGenes = permute.permute_mutation_data(G, genes, patients, self.seeds[i], Q)
            print i, " matrices generated.", self.seeds[i]

            for gene in geneToCases:
                self.geneToCases_perm[gene].append(geneToCases[gene])

            for patient in patientToGenes:
                self.patientToGenes_perm[patient].append(patientToGenes[patient])


            # Make directory to hold the temporary files.
            if matrixdirectory:
                adj_list = [ p + "\t" + "\t".join( sorted(patientToGenes[p]) ) for p in patients ]

                permutation_file = "{}/permuted-matrix-{}.m2".format(matrixdirectory, i+1)
                with open(permutation_file, 'w') as outfile: outfile.write('\n'.join(adj_list))

            # Clear dictionaries from stack once they're not used
            geneToCases.clear()
            patientToGenes.clear()

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
    minFreq = 2
    num_permutations = 100
    Q = 100
    matrixdirectory = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/matrix'
    outmutexfile = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/confirmPRAD_testoutmutex' + str(num_permutations) + 'Q' + str(Q) + str(time.time()) + '.tsv'
    outcooccurfile = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/confirmPRAD_testoutcooccur' + str(num_permutations) + 'Q' + str(Q) + str(time.time()) + '.tsv'
    outseedsfile = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/confirmPRAD_seeds' + str(time.time()) + '.tsv'

    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)

    for patient in patients:
        if not patientToGenes[patient]:
            patientToGenes.pop(patient)
            print patient, "popped"

    # Generate Permutation Matrices
    pm = PermutationMatrices(geneToCases, patientToGenes, num_permutations, Q=Q, matrixdirectory=matrixdirectory)

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
        fieldnames = pair_to_mutex[genepair].keys()
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
        fieldnames = pair_to_cooccur[genepair].keys()
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