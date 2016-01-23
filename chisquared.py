__author__ = 'jlu96'

import mutex as mex
import mutex_triangles as met
import csv
from scipy import stats


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



def add_ChiP_all_pairs(pairsdict, geneToCases, patientToGenes):

    for pair in pairsdict:
        geneset = set(pair)
        pairsdict[pair]['Chi-SquaredProbability'] = get_pair_ChiP(geneset, geneToCases, patientToGenes)

    return pairsdict





# Binomial Test Functions

def get_pair_BinomP(geneset, geneToCases, patientToGenes):
    """
    :param geneset:
    :param geneToCases:
    :param patientToGenes:

    Test the hypothesis that the paired probability is equal to the expected probability

    :return:
    """

    numCases = len(patientToGenes)

    patientlist = [set(geneToCases[gene]) for gene in geneset]

    # Overlapping patients
    overlap_patients = set.intersection(*patientlist)

    # Marginal probability of genes
    p_gs = [len(geneToCases[gene]) * 1.0 / numCases for gene in geneset]

    # Probability of overlap of genes
    p_g_overlap = 1.
    for p_g in p_gs:
        p_g_overlap *= p_g

    # Get p-value
    p = stats.binom_test(len(overlap_patients), n=numCases, p=p_g_overlap)


    return p




def add_BinomP_all_pairs(pairsdict, geneToCases, patientToGenes):

    for pair in pairsdict:
        geneset = set(pair)
        pairsdict[pair]['BinomProbability'] = get_pair_BinomP(geneset, geneToCases, patientToGenes)

    return pairsdict









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
    p = stats.binom_test(len(overlap_patients), n=numCases, p=p_g_overlap)

    return p



def add_BinomP_cohorts_all_pairs(pairsdict, geneToCases, patientToGenes, cohort_dict, pvalue=0.05):

    num_cohorts = len(cohort_dict)

    for pair in pairsdict:
        geneset = set(pair)
        all_sig = True
        for cohort_num in cohort_dict:
            BinomP = get_pair_BinomP_cohort(geneset, geneToCases, patientToGenes, cohort_dict[cohort_num])
            pairsdict[pair][str(num_cohorts) + 'BinomProbability' + str(cohort_num)] = BinomP
            all_sig = all_sig and (BinomP < pvalue)
        pairsdict[pair][str(num_cohorts) + 'AllSig'] = all_sig


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


    mutationmatrix = '/Users/jlu96/maf/new/PRAD_broad/PRAD_broad-som.m2'
    patientFile = '/Users/jlu96/maf/new/PRAD_broad/shared_patients.plst'
    cpairfile = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/LorenzoModel/Binomial/PRAD_broad-som-Cpairs.txt'

    geneFile = None
    minFreq = 0
    test_minFreq = 5
    num_cohorts_list = [1,2,3]


    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)

    print "number of genes is ", numGenes

    # patientToChi = get_patientToChiP(patients, geneToCases, patientToGenes)
    # patient_chi_file = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/Chi-Squared/PRAD_broad-som-patients.txt'
    # with open(patient_chi_file, 'w') as pcf:
    #     writer = csv.writer(pcf, delimiter='\t')
    #     writer.writerow(['Patient', 'Chi-Squared', 'MutationFrequency'])
    #     for patient, chi in patientToChi.items():
    #         writer.writerow([patient, chi, len(patientToGenes[patient])])


    genepairs = met.getgenepairs(geneToCases, genes, test_minFreq=test_minFreq)
    print "Number of pairs ", len(genepairs)


    print "Normal cooccur test"
    cpairsdict, cgenedict = met.cooccurpairs(numCases, geneToCases, patientToGenes, genepairs)

    # print "Add binomial probability"
    # cpairsdict = add_BinomP_all_pairs(cpairsdict, geneToCases, patientToGenes)

    # undo
    # print "Beginning cohorts"
    #
    #
    # for num_cohorts in num_cohorts_list:
    #     # get cohorts
    #     cohort_dict = generate_patient_cohorts(patientToGenes, num_cohorts)
    #
    #     cpairsdict = add_BinomP_cohorts_all_pairs(cpairsdict, geneToCases, patientToGenes, cohort_dict)
    #
    #     cpairsdict = add_cohorts_all_pairs(cpairsdict, geneToCases, patientToGenes, cohort_dict)
    #
    # print "Writing to file..."
    met.writeanydict(cpairsdict, cpairfile)



if __name__ == '__main__':
    main()