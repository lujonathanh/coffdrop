__author__ = 'jlu96'

import mutex as mex
import mutex_triangles as met
import csv
from scipy import stats
import partition
import time


def get_patient_ChiP(patient, geneToCases, patientToGenes):
    """
    :param patient: 
    :param numGenes: 
    :param numCases: 
    :param geneToCases: 
    :param patientToGenes: 
    :return: The chi-squared value of the sample, given the expected probabilities of the genes
    """
    patient_genes = patientToGenes[patient]
    numCases = len(patientToGenes)

    f_obs = [1. if gene in patient_genes else 0. for gene in geneToCases]


    # The expected value is the marginal probability of the gene's occurrence
    f_exp = [len(geneToCases[gene]) * 1.0 / numCases for gene in geneToCases]



    chisq, p = stats.chisquare(f_obs, f_exp)

    if p < 0.05:
        print patient
        # print "Observed: "
        # print f_obs[0:50]
        # print "Expected: "
        # print f_exp[0:50]

    return p


    

def get_patientToChiP(patients, geneToCases, patientToGenes):
    """
    :param patients:
    :param geneToCases:
    :param patientToGenes:
    :return: patientToChi
    """
    patientToChi = {}

    for patient in patients:
        chi = get_patient_ChiP(patient, geneToCases, patientToGenes)
        patientToChi[patient] = chi

    return patientToChi

    

def get_pair_ChiP(geneset, geneToCases, patientToGenes):
    """
    :param geneset: 
    :param geneToCases: 
    :param patientToGenes:

    The overlap probability for a sample is p_g_i * p_g_j (marginal probabilities of genes i and j involved in overlap)
                                * p_s_k ** 2 (marginal probability of sample k mutation)

    :return: The chi-squared value of the gene set, given the overlap events
    """
    set_size = len(geneset)


    numGenes = len(geneToCases)

    patients = patientToGenes.keys()
    patientlist = [set(geneToCases[gene]) for gene in geneset]

    # Overlapping patients
    overlap_patients = set.intersection(*patientlist)

    # Marginal probability of genes
    p_gs = [len(geneToCases[gene]) * 1.0 / len(patients) for gene in geneset]

    # Probability of overlap of genes
    p_g_product = 1.
    for p_g in p_gs:
        p_g_product *= p_g


    f_obs = [1. if patient in overlap_patients else 0. for patient in patients]

    f_exp = []
    for patient in patients:
        p_s = len(patientToGenes[patient]) * 1.0 / numGenes
        p_exp = p_g_product * (p_s ** set_size)
        f_exp.append(p_exp)

    print f_obs[0:50]
    print f_exp[0:50]

    chisq, p = stats.chisquare(f_obs, f_exp)

    print chisq, p
    return p


def get_triplet_ChiP(gene_triplet, geneToCases, patientToGenes):
    gene0, gene1, gene2 = gene_triplet
    gene_list = list(gene_triplet)

    gene_triplet = set(gene_triplet)
    numGenes = len(geneToCases)
    patients = patientToGenes.keys()
    numCases = len(patients)

    #Marginal probability of genes
    p_gs = [len(geneToCases[gene]) * 1.0 / len(patients) for gene in gene_list]
    print "Marginal probs ", p_gs

    f_obs = []
    f_exp = []

    # All those with all of them
    p_alls = []
    p_all = 1.0
    for p_g in p_gs:
        p_all *= p_g
    p_alls.append(p_all)

    f_obs.append(len([p for p in patientToGenes if gene_triplet.issubset(patientToGenes[p])]))
    f_exp.append(numCases * 1.0 * p_all)

    # All those with two of them
    for geneA in gene_list:
        index = gene_list.index(geneA)
        other_genes = gene_triplet.difference(set([geneA]))
        geneB, geneC = other_genes

        #print "In gene ", geneA, "Second gene ", geneB, "last gene ", geneC

        p_all = 1.0
        for i in range(3):
            if i == index:
                p_all *= 1.0 - p_gs[i]
            else:
                p_all *= p_gs[i]
        f_obs.append(len([p for p in patientToGenes if other_genes.issubset(patientToGenes[p]) and geneA not in patientToGenes[p]]))
        f_exp.append(numCases * 1.0 * p_all)

        p_alls.append(p_all)


    # All those with one of them

    #print "only contain 1"
    for geneA in gene_list:
        index = gene_list.index(geneA)
        other_genes = gene_triplet.difference(set([geneA]))
        geneB, geneC = other_genes


        #print "In gene ", geneA, "Second gene ", geneB, "last gene ", geneC
        p_all = 1.0
        for i in range(3):
            if i != index:
                p_all *= 1.0 - p_gs[i]
            else:
                p_all *= p_gs[i]
        f_obs.append(len([p for p in patientToGenes if not other_genes.intersection(patientToGenes[p]) and geneA in patientToGenes[p]]))
        f_exp.append(numCases * 1.0 * p_all)

        p_alls.append(p_all)

    # All those with none
    p_all = 1.0
    for p_g in p_gs:
        p_all *= 1.0 - p_g
    p_alls.append(p_all)

    f_obs.append(len([p for p in patientToGenes if not gene_triplet.intersection(patientToGenes[p])]))
    f_exp.append(numCases * 1.0 * p_all)

    print "Observed ", f_obs, "expected ", f_exp
    print "Probabilities ", p_alls
    print "total probability ", sum(p_alls)

    chisq, p = stats.chisquare(f_obs, f_exp)
    print "chi-squared test "
    print chisq, p
    return p






def add_ChiP_all_pairs(pairsdict, geneToCases, patientToGenes):

    for pair in pairsdict:
        geneset = set(pair)
        pairsdict[pair]['Chi-SquaredProbability'] = get_pair_ChiP(geneset, geneToCases, patientToGenes)

    return pairsdict





# Binomial Test Functions
#
# def get_pair_BinomP(geneset, geneToCases, patientToGenes):
#     """
#     :param geneset:
#     :param geneToCases:
#     :param patientToGenes:
#
#     Test the hypothesis that the paired probability is equal to the expected probability
#
#     :return:
#     """
#
#     numCases = len(patientToGenes)
#
#     patientlist = [set(geneToCases[gene]) for gene in geneset]
#
#     # Overlapping patients
#     overlap_patients = set.intersection(*patientlist)
#
#     # Marginal probability of genes
#     p_gs = [len(geneToCases[gene]) * 1.0 / numCases for gene in geneset]
#
#     # Probability of overlap of genes
#     p_g_overlap = 1.
#     for p_g in p_gs:
#         p_g_overlap *= p_g
#
#
#
#     # cooccurprob = stats.binom.sf(len(overlap_patients) - 1, numCases,   p_g_overlap)
#     # mutexprob = stats.binom.cdf(len(overlap_patients), numCases,   p_g_overlap)
#
#
#     return p


# def get_triplet_BinomP(geneset, geneToCases, patientToGenes, cpairs, mpairs, triplet_type='CooccurringMutuallyExclusiveMutuallyExclusive'):
#     numCases = len(patientToGenes)
#
#     if triplet_type == 'CooccurringMutuallyExclusiveMutuallyExclusive':
#         cpair = cpairs[0]
#         g_m = geneset.difference(set(cpair)).pop()
#         g_c1, g_c2 = cpair
#
#         p_g_m = len(geneToCases[g_m]) * 1.0 / numCases
#         #p_g_c = len([p for p in patientToGenes if g_c1 in patientToGenes[p] and g_c2 in patientToGenes[p]]) * 1.0 / numCases
#         p_g_c1 = len(geneToCases[g_c1]) * 1.0 / numCases
#         p_g_c2 = len(geneToCases[g_c2]) * 1.0 / numCases
#         # Claculate probaiblity of CMM: p(g_c1)p(g_c2)(1-p(g_m) + (1-p(g_c1))(1-p(g_2))p(g_m)
#
#
#
#         p_cmm = p_g_c1 * p_g_c2 * (1 - p_g_m) + (1 - p_g_c1) * (1 - p_g_c2) * p_g_m
#
#         overlap_patients = [p for p in patientToGenes if
#                             (g_m in patientToGenes[p] and g_c1 not in patientToGenes[p] and g_c2 not in patientToGenes[p])
#                             or (g_m not in patientToGenes[p] and g_c1 in patientToGenes[p] and g_c2 in patientToGenes[p])]
#
#         p = 0.5 * stats.binom_test(len(overlap_patients), n=numCases, p=p_cmm)
#
#         print "Coccur pair : ", len(geneToCases[g_c1]), len(geneToCases[g_c2]), "Mutex: ", len(geneToCases[g_m])
#
#         print "p_g_cs", p_g_c1, p_g_c2, "p_g_m", p_g_m
#         print "p_cmm ", p_cmm, "expected ", p_cmm * numCases
#
#         print "observed prop", len(overlap_patients)*1.0/numCases, "num overlap ", len(overlap_patients), "out of ", numCases, "total "
#
#         return p
#
# # Make sure to only do tail probabilities
# def get_triplet_BinomP_ab(geneset, geneToCases, patientToGenes, cpairs, mpairs, triplet_type='CooccurringMutuallyExclusiveMutuallyExclusive'):
#     numCases = len(patientToGenes)
#
#     if triplet_type == 'CooccurringMutuallyExclusiveMutuallyExclusive':
#         cpair = cpairs[0]
#         g_m = geneset.difference(set(cpair)).pop()
#         g_c1, g_c2 = cpair
#
#         p_g_m = len(geneToCases[g_m]) * 1.0 / numCases
#         p_g_c = len([p for p in patientToGenes if g_c1 in patientToGenes[p] and g_c2 in patientToGenes[p]]) * 1.0 / numCases
#         # p_g_c1 = len(geneToCases[g_c1]) * 1.0 / numCases
#         # p_g_c2 = len(geneToCases[g_c2]) * 1.0 / numCases
#         # Claculate probaiblity of CMM: p(g_c1)p(g_c2)(1-p(g_m) + (1-p(g_c1))(1-p(g_2))p(g_m)
#
#         p_cmm = p_g_m * (1 - p_g_c) + p_g_c * (1 - p_g_m)
#
#
#         # p_cmm = p_g_c1 * p_g_c2 * (1 - p_g_m) + (1 - p_g_c1) * (1 - p_g_c2) * p_g_m
#
#         overlap_patients = [p for p in patientToGenes if
#                             (g_m in patientToGenes[p] and g_c1 not in patientToGenes[p] and g_c2 not in patientToGenes[p])
#                             or (g_m not in patientToGenes[p] and g_c1 in patientToGenes[p] and g_c2 in patientToGenes[p])]
#
#         p = 0.5 * stats.binom_test(len(overlap_patients), n=numCases, p=p_cmm)
#
#         print "Coccur pair : ", len(geneToCases[g_c1]), len(geneToCases[g_c2]), "Mutex: ", len(geneToCases[g_m])
#
#         print "p_g_cs", p_g_c, "p_g_m", p_g_m
#         print "p_cmm ", p_cmm, "expected ", p_cmm * numCases
#
#         print "observed prop", len(overlap_patients)*1.0/numCases, "num overlap ", len(overlap_patients), "out of ", numCases, "total "
#
#         return p






# def add_BinomP_all_pairs(pairsdict, geneToCases, patientToGenes):
#
#     for pair in pairsdict:
#         geneset = set(pair)
#         pairsdict[pair]['BinomProbability'] = get_pair_BinomP(geneset, geneToCases, patientToGenes)
#
#     return pairsdict
#
#
#
#





def get_pair_BinomP_cohort(geneset, geneToCases, patientToGenes, cohort):
    """
    :param geneset:
    :param geneToCases:
    :param patientToGenes:
    :param patient_cohorts:
    :return: Binomial Probability under each cohorts
    """

    numCases = len(cohort)

    patientlist = [set(geneToCases[gene]).intersection(cohort) for gene in geneset]

    # Overlapping patients
    overlap_patients = set.intersection(*patientlist)

    # Marginal probability of genes
    p_gs = [len(geneToCases[gene].intersection(cohort)) * 1.0 / numCases for gene in geneset]

    # Probability of overlap of genes
    p_g_overlap = 1.
    for p_g in p_gs:
        p_g_overlap *= p_g

    # Get p-value
    # rv = stats.binom(numCases, p_g_overlap)
    cooccurprob = stats.binom.sf(len(overlap_patients) - 1, numCases,   p_g_overlap)
    mutexprob = stats.binom.cdf(len(overlap_patients), numCases,   p_g_overlap)
    return numCases, [len(p) for p in patientlist], len(overlap_patients), cooccurprob, mutexprob



def add_BinomP_cohorts_all_pairs(pairsdict, geneToCases, patientToGenes, cohort_dict, pvalue=0.05):

    num_cohorts = len(cohort_dict)

    for pair in pairsdict:
        geneset = set(pair)
        all_c_sig = True
        all_m_sig = True

        # Do over all of them
        cohort_size, freqs, overlap, cprob, mprob = get_pair_BinomP_cohort(geneset, geneToCases, patientToGenes, patientToGenes.keys())
        pairsdict[pair]['AllSize'] = cohort_size
        pairsdict[pair]['AllFreqs'] = freqs
        pairsdict[pair]['AllOverlap'] = overlap
        pairsdict[pair]['AllCBinomProb'] = cprob
        pairsdict[pair]['AllMBinomProb'] = mprob

        # Then over each of the individual cohorts
        for cohort_num in cohort_dict:
            cohort_size, freqs, overlap, cprob, mprob = get_pair_BinomP_cohort(geneset, geneToCases, patientToGenes, cohort_dict[cohort_num])
            pairsdict[pair][str(num_cohorts) + 'Size' + str(cohort_num)] = cohort_size
            pairsdict[pair][str(num_cohorts) + 'Freqs' + str(cohort_num)] = freqs
            pairsdict[pair][str(num_cohorts) + 'Overlap' + str(cohort_num)] = overlap
            pairsdict[pair][str(num_cohorts) + 'CBinomProb' + str(cohort_num)] = cprob
            pairsdict[pair][str(num_cohorts) + 'MBinomProb' + str(cohort_num)] = mprob
            all_c_sig = all_c_sig and (cprob < pvalue)
            all_m_sig = all_m_sig and (mprob < pvalue)
        pairsdict[pair][str(num_cohorts) + 'CAllSig'] = all_c_sig
        pairsdict[pair][str(num_cohorts) + 'MAllSig'] = all_m_sig

    return pairsdict




"""Add all the information within the cohorts."""

def add_cohorts_all_pairs(pairsdict, geneToCases, patientToGenes, cohort_dict):

    num_cohorts = len(cohort_dict)

    for pair in pairsdict:
        geneset = set(pair)
        for cohort_num in cohort_dict:
            pairsdict[pair][str(num_cohorts) + 'CohortObservedPairs' + str(cohort_num)], \
            pairsdict[pair][str(num_cohorts) + 'CohortExpectedPairs' + str(cohort_num)] = analyze_pair_cohort(geneset, geneToCases, patientToGenes, cohort_dict[cohort_num])

            cohort_size = len(cohort_dict[cohort_num])
            pairsdict[pair][str(num_cohorts) + 'CohortObservedRatio' + str(cohort_num)] = pairsdict[pair][str(num_cohorts) + 'CohortObservedPairs' + str(cohort_num)] * 1.0 / cohort_size
            pairsdict[pair][str(num_cohorts) + 'CohortExpectedRatio' + str(cohort_num)] = pairsdict[pair][str(num_cohorts) + 'CohortExpectedPairs' + str(cohort_num)] * 1.0 / cohort_size
    return pairsdict



def add_BinomP_min_cohort_all_pairs(pairsdict, geneToCases, patientToGenes, cohort_dict, min_cohort, pvalue=0.05,
                                    min_p_thresh=0.05):

    num_cohorts = len(cohort_dict)

    # limit to the pairs that are present in min_cohort
    min_cohort_genes = set.union(*(patientToGenes[p] for p in min_cohort))
    min_cohort_pairs = set([pair for pair in pairsdict if set(pair).issubset(min_cohort_genes)])
    print "Original pairs ", len(pairsdict), ". Min cohort pairs: ", len(min_cohort_pairs)
    min_pthresh = pvalue * 1.0 / len(min_cohort_pairs) # the multiple-testing corrected pvalue threshold

    for pair in pairsdict:
        geneset = set(pair)
        all_c_sig = True
        all_m_sig = True

        if pair in min_cohort_pairs:
            # Add their value in the minimum binom prob
            cohort_size, freqs, overlap, cprob, mprob = get_pair_BinomP_cohort(geneset, geneToCases, patientToGenes, min_cohort)
            pairsdict[pair][str(num_cohorts) + 'Size' + 'Min'] = cohort_size
            pairsdict[pair][str(num_cohorts) + 'Freqs' + 'Min'] = freqs
            pairsdict[pair][str(num_cohorts) + 'Overlap' + 'Min'] = overlap
            pairsdict[pair][str(num_cohorts) + 'Min' + 'CBinomProb'] = cprob
            pairsdict[pair][str(num_cohorts) + 'Min' + 'MBinomProb'] = mprob
            pairsdict[pair]['MinSig'] = 1 if ((cprob < min_pthresh) or (mprob < min_pthresh)) else 0
        else:
            pairsdict.pop(pair)

    return pairsdict


def analyze_pair_cohort(geneset, geneToCases, patientToGenes, cohort):
    """ Return the observed and expected number of overlapping pairs within this cohort
    """

    numCases = len(cohort)
    # Actual number of pairs in cohort

    patientlist = [set(geneToCases[gene]).intersection(cohort) for gene in geneset]

    overlap_patients = set.intersection(*patientlist)

    obs = len(overlap_patients)

    # Expected number of pairs in cohort

        # Marginal probability of genes
    p_gs = [len(geneToCases[gene].intersection(cohort)) * 1.0 / numCases for gene in geneset]

        # Probability of overlap of genes
    p_g_overlap = 1.
    for p_g in p_gs:
        p_g_overlap *= p_g

    exp = numCases * p_g_overlap

    return obs, exp




# def get_expected_overlaps_cohort(geneset, geneToCases, patientToGenes, cohort):
#     """
#     :param geneset:
#     :param geneToCases:
#     :param patientToGenes:
#     :param patient_cohorts:
#     :return: obs_overlap, exp_overlap
#     """
#
#
#
#     patientlist = [set(geneToCases[gene]) for gene in geneset]
#
#     # Overlapping patients
#     overlap_patients = set.intersection(cohort, *patientlist)
#
#     obs_overlap = len(overlap_patients)
#
#
#
#
#
#
# def add_expected_pairs_cohorts(geneset, geneToCases, patientToGenes, patient_cohorts):




# Patient Cohort Functions

def generate_patient_cohorts(patientToGenes, num_cohorts):
    """
    :param patientToGenes:
    :param num_cohorts:
    :return: patient_cohort_dict: Dictionary of patient cohorts
    """

    sorted_patients = sorted(patientToGenes.keys(), key= lambda patient: len(patientToGenes[patient]))
    numCases = len(patientToGenes)

    patient_cohorts = [sorted_patients[i * numCases / num_cohorts: (i+1) * numCases/num_cohorts] for i in range(num_cohorts)]

    patient_cohort_dict = {}

    for i in range(len(patient_cohorts)):
        patient_cohort_dict[i] = patient_cohorts[i]

    return patient_cohort_dict



# def add_ChiP_all_pairs_cohorts(pairsdict, geneToCases, patientToGenes, patient_cohorts):






def main():


    mutationmatrix = '/Users/jlu96/maf/new/OV_broad/OV_broad-cna-jl.m2'
    patientFile = '/Users/jlu96/maf/new/OV_broad/shared_patients.plst'
    cpairfile = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/LorenzoModel/Binomial/OV_broad-cna-jl-cpairs-min_cohort.txt'
    partitionfile = '/Users/jlu96/maf/new/OV_broad/OV_broad-cna-jl.ppf'
    load_partitions = True
    do_min_cohort = True

    geneFile = None
    minFreq = 0
    test_minFreq = 100
    compute_mutex = True



    include_cohort_info = False
    num_cohorts_list = [1,3, 5, 7]


    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)

    print "number of genes is ", numGenes


    if do_min_cohort:
        cohort_dict, clusterToProp, min_cohort = partition.load_patient_cohorts(partitionfile, patientToGenes)
        min_cohort_genes = set.union(*(patientToGenes[p] for p in min_cohort))

        print "getting pairs"
        genepairs = met.getgenepairs(geneToCases, min_cohort_genes, test_minFreq=test_minFreq)

        print "Number of pairs ", len(genepairs)


        print "Normal cooccur test"
        t = time.time()
        cpairsdict, cgenedict = met.cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, compute_mutex=compute_mutex)
        print "Normal cooccur done in ", time.time() - t

        print "Beginning cohorts"
        t = time.time()
        cpairsdict = add_BinomP_min_cohort_all_pairs(cpairsdict, geneToCases, patientToGenes, cohort_dict, min_cohort)
        print "Cohorts done in ", time.time() - t

    else:
        genepairs = met.getgenepairs(geneToCases, genes, test_minFreq=test_minFreq)
        print "Number of pairs ", len(genepairs)


        print "Normal cooccur test"
        cpairsdict, cgenedict = met.cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, compute_mutex=compute_mutex)

        # print "Add binomial probability"
        # cpairsdict = add_BinomP_all_pairs(cpairsdict, geneToCases, patientToGenes)

        # undo
        print "Beginning cohorts"





        if load_partitions:
            cohort_dict = partition.load_patient_cohorts(partitionfile)
            cpairsdict = add_BinomP_cohorts_all_pairs(cpairsdict, geneToCases, patientToGenes, cohort_dict)

        else:
            for num_cohorts in num_cohorts_list:
                # get cohorts
                cohort_dict = generate_patient_cohorts(patientToGenes, num_cohorts)

                cpairsdict = add_BinomP_cohorts_all_pairs(cpairsdict, geneToCases, patientToGenes, cohort_dict)

                if include_cohort_info:
                    cpairsdict = add_cohorts_all_pairs(cpairsdict, geneToCases, patientToGenes, cohort_dict)

    print "Writing to file..."
    met.writeanydict(cpairsdict, cpairfile)



if __name__ == '__main__':
    main()