__author__ = 'jlu96'

import random
from scipy import stats

class CooccurPatientDistribution:
    def __init__(self, patients, patientToGenes, dist_num=1000):
        self.patients = patients
        self.patientToGenes = patientToGenes

        self.mut_freqs = [len(patientToGenes[patient]) for patient in patients]

        # Dictionary with low_mutated_score distributions
        self.low_mutated_score_dist_dict = {}
        self.dist_num = dist_num


    def calc_low_mutated_score_dist(self, numPatients):

        low_mutated_score_dist = []

        for i in range(self.dist_num):
            set_patients = random.sample(self.patients, numPatients)
            score = low_mutated_score(set_patients, self.patientToGenes)
            low_mutated_score_dist.append(score)


        self.low_mutated_score_dist_dict[numPatients] = low_mutated_score_dist


    def get_percentile(self, low_mutated_score, numPatients):
        if numPatients not in self.low_mutated_score_dist_dict:
            self.calc_low_mutated_score_dist(numPatients)


        return stats.percentileofscore(self.low_mutated_score_dist_dict[numPatients], low_mutated_score)



def low_mutated_score(patients, patientToGenes):
    mut_freqs = [len(patientToGenes[patient]) for patient in patients]

    return sum([1.0/mut_freq for mut_freq in mut_freqs])



def add_low_mutated_scores(cpairsdict, geneToCases, patientToGenes, dist_num=1000):
    """
    :return: The co-occurring pairs with low mutated scores and percentiles added.
    """

    cooccur_patients = set()

    for pair in cpairsdict:
        gene0, gene1 = cpairsdict[pair]['Gene0'], cpairsdict[pair]['Gene1']
        shared_patients = geneToCases[gene0].intersection(geneToCases[gene1])
        cpairsdict[pair]['LowMutatedScore'] = low_mutated_score(shared_patients, patientToGenes)

        cooccur_patients = cooccur_patients.union(shared_patients)

    CPD = CooccurPatientDistribution(cooccur_patients, patientToGenes, dist_num=dist_num)

    for pair in cpairsdict:
        cpairsdict[pair]['LowMutatedScorePercentile'] = CPD.get_percentile(cpairsdict[pair]['LowMutatedScore'], cpairsdict[pair]['Overlap'])

    return cpairsdict