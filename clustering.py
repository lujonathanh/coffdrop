__author__ = 'jlu96'

import mutex as mex
import scipy.cluster.hierarchy as hac
import numpy as np
import matplotlib.pyplot as plt
import time
import csv
from scipy import stats


def distance_matrix(keyToBinary, keys = None, distance_type='Jaccard'):
    if not keys:
        keys = keyToBinary.keys()

    matrix = [[0 for j in range(len(keys))] for i in range(len(keys))]

    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if distance_type == 'Jaccard':
                distance = mex.cooccurrence_ratio([keys[i], keys[j]], keyToBinary)
                matrix[i][j] = distance
                matrix[j][i] = distance

    return matrix

def patient_distance_matrix(patientToGenes, patients = None, distance_type='Jaccard'):
    if not patients:
        patients = patientToGenes.keys()

    return distance_matrix(patientToGenes, patients, distance_type)

def gene_distance_matrix(geneToCases, genes = None, distance_type='Jaccard'):
    if not genes:
        genes = geneToCases.keys()

    return distance_matrix(geneToCases, genes, distance_type)

def write_matrix(filename, matrix, columns):
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(columns)
        for row in matrix:
            writer.writerow(row)

def gene_freq_association(gene, patients, geneToCases, patientToGenes):
    in_patients = geneToCases[gene]
    out_patients = set(patients).difference(in_patients)

    in_freq = [len(patientToGenes[patient]) for patient in in_patients]
    in_avg = sum(in_freq)/len(in_freq)
    out_freq = [len(patientToGenes[patient]) for patient in out_patients]
    out_avg = sum(out_freq)/len(out_freq)

    u, p = stats.mannwhitneyu(in_freq, out_freq)

    return p, in_avg, out_avg

def all_gene_freq_association(geneToCases, patientToGenes):
    genes = geneToCases.keys()
    patients = patientToGenes.keys()

    genePValue = set()

    # get all gene frequency associations
    for gene in genes:
        p, in_avg, out_avg = gene_freq_association(gene, patients, geneToCases, patientToGenes)
        genePValue.add((gene, p, len(geneToCases[gene]), in_avg, out_avg))

    genePValue = sorted(genePValue, key= lambda a: a[1])

    return genePValue

def main():

    mutationmatrix = '/Users/jlu96/maf/new/BRCA_wustl/BRCA_wustl-som-cna-jl.m2'
    patientFile = '/Users/jlu96/maf/new/BRCA_wustl/shared_patients.plst'
    geneFile = None
    minFreq = 0
    outpatientfile = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/LorenzoModel/Heatmaps/BRCA_patientmatrix_log.csv'

    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)


    # ---- All Gene Frequency Association -------------------------------------------------------------------------
    genePvalue_file = '/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/LorenzoModel/geneMakeBin/BRCA_genemakebin.csv'
    genePValue = all_gene_freq_association(geneToCases, patientToGenes)

    with open(genePvalue_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Gene', 'Pvalue', 'Frequency', 'InAvg', 'OutAvg'])
        for entry in genePValue:
            writer.writerow(entry)




    # ---- Patient Distance Matrix --------------------------------------------------------------------------------
    # t = time.time()
    # pdm = np.array(patient_distance_matrix(patientToGenes))
    # print "Done with patient distance matrix in time ", time.time() - t
    #
    # write_matrix(outpatientfile, -1 * np.log(pdm + 1e-10), patientToGenes.keys())



    #----------------------------------------------------------------------------------------------------------------

    # gdm = np.array(gene_distance_matrix(geneToCases))

    # code for plotting
    # fig, axes23 = plt.subplots(2,3)
    #
    # z_pdm = hac.linkage(pdm)
    #
    # for axes in axes23:
    #     # Plotting
    #     axes[0].plot(range(1, len(z_pdm)+1), z_pdm[::-1, 2])
    #     knee = np.diff(z_pdm[::-1, 2], 2)
    #     axes[0].plot(range(2, len(z_pdm)), knee)
    #
    #     num_clust1 = knee.argmax() + 2
    #     knee[knee.argmax()] = 0
    #     num_clust2 = knee.argmax() + 2
    #
    #     axes[0].text(num_clust1, z_pdm[::-1, 2][num_clust1-1], 'possible\n<- knee point')
    #
    #     part1 = hac.fcluster(z_pdm, num_clust1, 'maxclust')
    #     part2 = hac.fcluster(z_pdm, num_clust2, 'maxclust')
    #
    #     #clr = ['#2200CC' ,'#D9007E' ,'#FF6600' ,'#FFCC00' ,'#ACE600' ,'#0099CC' ,
    #     #'#8900CC' ,'#FF0000' ,'#FF9900' ,'#FFFF00' ,'#00CC01' ,'#0055CC']
    #
    #     for part, ax in zip([part1, part2], axes[1:]):
    #         for cluster in set(part):
    #             ax.scatter(pdm[part == cluster, 0], pdm[part == cluster, 1])
    #
    #     m = '\n(method: {})'.format('normal')
    #     plt.setp(axes[0], title='Screeplot{}'.format(m), xlabel='partition',
    #              ylabel='{}\ncluster distance'.format(m))
    #     plt.setp(axes[1], title='{} Clusters'.format(num_clust1))
    #     plt.setp(axes[2], title='{} Clusters'.format(num_clust2))
    #
    # plt.tight_layout()
    # plt.show()

if __name__ == '__main__':
    main()