__author__ = 'jlu96'
import mutex as mex
import matplotlib.pyplot as plt
import csv
import numpy as np
from sklearn import mixture
from sklearn.cluster import KMeans
from scipy.stats import poisson
import collections
import os

def partition_EM(patientToGenes, k):
    """
    :param geneToCases:
    :param patientToGenes:
    :param k: Number of partitions
    :return: cohort_list
    """

    # partition the patients, and intersect the geneToCases
    return



def partition_gene(patientToGenes, genes):
    """
    :param geneToCases:
    :param patientToGenes:
    :param genes:
    :return: cohorts by each gene. Size 2^(#genes)
    """

    cohorts = [patientToGenes.keys()]
    for gene in genes:
        new_cohorts = []
        for cohort in cohorts:
            new_cohort_1 = [patient for patient in patientToGenes if gene not in patientToGenes[patient]]
            if new_cohort_1:
                new_cohorts.append(new_cohort_1)
            new_cohort_2 = list(set(cohort).difference(set(new_cohort_1)))
            if new_cohort_2:
                new_cohorts.append(new_cohort_2)
        cohorts = new_cohorts
    # print genes
    # print cohorts

    return cohorts

def partition_gene_list(patientToGenes, genes, binary=True):
    """
    :param patientToGenes:
    :param genes:
    :return: The cohorts, ordered from least to greatest in number of those genes they have.
    If binary = True, return just those with, those without.
    """

    gene_set = set(genes)
    cohort_dict = {}

    for patient in patientToGenes:
        num = len(set.intersection(gene_set, patientToGenes[patient]))

        # just 0 and 1
        if binary:
            if num > 0:
                num = 1

        if num not in cohort_dict:
            cohort_dict[num] = []
        cohort_dict[num].append(patient)


    return cohort_dict




def draw_partitions_cohorts(geneToCases, patientToGenes, cohort_pairings, title=None, num_bins=50):
    # LEFT OF HERE, JLU. Finish this, then above. Make plots in parallel, compare.
    # Work with: TP53? Others?

    numGenes = len(geneToCases.keys())
    numCohorts = len(cohort_pairings)

    cohort_frequencies = [[len(patientToGenes[case]) for case in cohort_pair[1]] for cohort_pair in cohort_pairings]
    cohort_names = [cohort_pair[0] for cohort_pair in cohort_pairings]

    draw_partitions(patientToGenes, cohort_names, cohort_frequencies, title=title, num_bins=num_bins)


def draw_partitions(patientToGenes, cohort_names, cohort_frequencies, title=None, num_bins=50):

    numCohorts = len(cohort_frequencies)
    bins = range(0, max([len(p_gene) for p_gene in patientToGenes.values()]), max([len(p_gene) for p_gene in patientToGenes.values()])/num_bins)

    plt.figure()


    for i in range(len(cohort_frequencies)):
        plt.hist(cohort_frequencies[i], bins, alpha=1.0/numCohorts, label=str(cohort_names[i]))


    plt.title(title, fontsize=20)
    plt.xlabel('# Somatic Mutations In Tumor', fontsize=20)
    plt.ylabel('Number of Samples', fontsize=20)
    plt.legend()
    plt.show()

def norm(x, height, center, std):
    return(height*np.exp(-(x - center)**2/(2*std**2)))

def integrate_cohorts(cohort_dict, numCases, num_integrated):
    cohorts_int = {}
    start_index = 0
    num_in_cohort = 0
    new_cohort = []
    for i in cohort_dict.keys():
        num_in_cohort += len(cohort_dict[i])
        new_cohort.extend(cohort_dict[i])
        if (num_in_cohort > numCases/num_integrated):
            cohorts_int[start_index] = new_cohort
            start_index = i+1
            new_cohort = []
            num_in_cohort = 0

    if new_cohort:
        cohorts_int[start_index] = new_cohort

    return cohorts_int

def partition_GMM(patientToGenes, num_components, num_bins, title=None, do_plot=True):
    g = mixture.GMM(n_components=num_components)
    mut_num_list = [len(patientToGenes[p]) for p in patientToGenes]
    obs = np.array([[entry] for entry in mut_num_list])
    g.fit(obs)

    print "***********************************"
    print "COMPONENTS: ", num_components
    print "Weights: " + str(np.round(g.weights_,2))
    print "Means: " + str(np.round(g.means_,2))
    print "Covariates: " + str(np.round(g.covars_,2))

    print "Total log probability: " + str(sum(g.score(obs)))
    print "AIC: " + str(g.aic(obs))
    print "BIC: ", g.bic(obs)

    score, respon = g.score_samples(obs)

    for i in range(num_components):
        print "Model ", np.round(g.means_, 2)[i], " explains ", np.round(len([in_w for in_w in respon if in_w[i] == max(in_w)])) * 1.0 /len(respon)


    # Simulate gaussians
    # sim_samples = g.sample(len(patientToGenes))
    bins = range(0, max([len(p_gene) for p_gene in patientToGenes.values()]), max([len(p_gene) for p_gene in patientToGenes.values()])/num_bins)
    histogram = np.histogram([len(patientToGenes[p]) for p in patientToGenes], bins=bins)

    # get the scale of the gaussians from the biggest one
    # max_comp = g.weights_.index(max(g.weights_))
    # max_mean = g.means_[max_comp]

    which_bins = [[bin for bin in bins if bin > mean][0] for mean in g.means_]
    print which_bins
    print bins
    print histogram
    print bins.index(which_bins[0]) - 1
    bin_heights = [histogram[0][bins.index(which_bin) - 1] for which_bin in which_bins]
    # max_height = max(histogram)

    if do_plot:
        plt.figure()
        plt.hist([len(patientToGenes[p]) for p in patientToGenes], bins=bins)
        for i in range(num_components):
            X = np.arange(0, max(mut_num_list), 1)
            Y = norm(X, bin_heights[i], g.means_[i], np.sqrt(g.covars_[i]))
            plt.plot(X, Y, label=str(np.round(g.weights_[i], 3)), linewidth=5)
        plt.title("GMM size " + str(num_components), fontsize=20)
        plt.xlabel('# Somatic Mutations In Tumor', fontsize=20)
        plt.ylabel('Number of Samples', fontsize=20)
        plt.legend()
        plt.show()
        # draw_partitions(patientToGenes, ['Original', 'Simulated'], [[len(patientToGenes[p]) for p in patientToGenes], sim_samples],
        #                 num_bins=num_bins, title=title)

    data = {}
    data['Components'] = num_components
    data['Weights'] = np.round(g.weights_,2)
    data['Means'] = np.round(g.means_,2)
    # data['Covariates'] = np.round(g.covars_,2)
    # data["Total log probability"] = sum(g.score(obs))
    data["AIC"] = g.aic(obs)
    data["BIC"] = g.bic(obs)
    data['Explained'] = [np.round([len([in_w for in_w in respon if in_w[i] == max(in_w)]) * 1.0 /len(respon) for i in range(num_components)], 2)]

    return data

def partition_gene_kmeans(geneToCases, patientToGenes, gene_list, num_components, num_bins, title=None, do_plot=True):

    # get gene index mapping
    giv = getgiv(geneToCases.keys(), gene_list)

    # convert patients into vectors
    patientToVector = getpatientToVector(patientToGenes, giv)

    vectors = patientToVector.values()

    print vectors[0]
    print "Length of vectors is ", len(vectors[0])

    km = KMeans(num_components)

    km.fit(vectors)

    clusterToPatient = {}

    for patient in patientToVector:
        cluster = km.predict(patientToVector[patient])[0]
        if cluster not in clusterToPatient:
            clusterToPatient[cluster] = set()
        clusterToPatient[cluster].add(patient)

    # plot patients in each cluster


    if do_plot:
        bins = range(0, max([len(p_gene) for p_gene in patientToGenes.values()]), max([len(p_gene) for p_gene in patientToGenes.values()])/num_bins)
        plt.figure()
        for cluster in clusterToPatient:
            plt.hist([len(patientToGenes[p]) for p in clusterToPatient[cluster]], bins=bins, label=str(cluster), alpha = 1.0/num_components)
        plt.xlabel('# Somatic Mutations In Tumor', fontsize=20)
        plt.ylabel('Number of Samples', fontsize=20)
        plt.legend()
        plt.title("Kmeans size " + str(num_components), fontsize=20)
        plt.show()



    data = {}
    data['Score'] = km.score(vectors)
    data['Number'] = num_components
    data['% Explained'] = np.round([100 * len(clusterToPatient[cluster]) * 1.0 / len(patientToGenes) for cluster in clusterToPatient], 2)
    data['Vector size'] = len(vectors[0])
    # data['Covariates'] = np.round(g.covars_,2)
    # data["Total log probability"] = sum(g.score(obs))
    # data["AIC"] = g.aic(obs)
    # data["BIC"] = g.bic(obs)
    # data['Explained'] = [np.round([len([in_w for in_w in respon if in_w[i] == max(in_w)]) * 1.0 /len(respon) for i in range(num_components)], 2)]

    return data


def getgiv(all_genes, gene_list):
    """
    :param all_genes:
    :param gene_list:
    :return: A list of the genes in common, the gene_index_vector.
    """
    giv = list(set(all_genes).intersection(set(gene_list)))

    return giv



def getpatientToVector(patientToGenes, gene_index_vector):
    patientToVector = {}
    for patient in patientToGenes:
        patient_genes = patientToGenes[patient]
        patientToVector[patient] = []
        for gene in gene_index_vector:
            patientToVector[patient].append(1 if gene in patient_genes else 0)

    return patientToVector



# 2/18/16 -Jlu

def best_pmm(patientToGenes, num_components, max_iter=30, rand_num=5, far_rand_num=5, min_cluster_size=0):

    data_record = []
    lls_record = []

    # Do normal
    first_data, lls = partition_pmm(patientToGenes, num_components,  max_iter=max_iter, min_cluster_size=min_cluster_size)

    data_record.append(first_data)
    lls_record.append(lls)

    # Do best rand init
    for i in range(rand_num):
        data, lls = partition_pmm(patientToGenes, num_components, rand_init=True, max_iter=max_iter, min_cluster_size=min_cluster_size)
        data_record.append(data)
        lls_record.append(lls)

    for i in range(far_rand_num):
        data, lls = partition_pmm(patientToGenes, num_components, far_rand_init=True, max_iter=max_iter, min_cluster_size=min_cluster_size)
        data_record.append(data)
        lls_record.append(lls)

    combined_record = zip(data_record, lls_record)

    combined_record = sorted(combined_record, key=lambda entry: entry[0]['Likelihood'], reverse=True)

    data_record, lls_record = zip(*combined_record)

    best_data = data_record[0]

    if (best_data['Likelihood'] > first_data['Likelihood'] + 10):
        print "First data not best!"
        best_data['IsFirst'] = False
    else:
        best_data['IsFirst'] = True


    clusterToPatient = pmm_to_cluster(patientToGenes, best_data['Classes'], best_data['Means'], best_data['Probabilities'])

    plot_pmm_clusters(patientToGenes, clusterToPatient, num_components)

    plot_likelihoods(lls_record)

    return best_data, clusterToPatient
    # Return clusters


def pmm_to_cluster(patientToGenes, classes, lam, p_k):
    clusterToPatient = {}

    for k in classes:
        clusterToPatient[k] = set()

    clusterToPatient[-1] = set()


    for patient in patientToGenes:
        d = len(patientToGenes[patient])

        max_class = -1
        max_ll = -np.inf
        for k in classes:
            if (np.log(p_k[k]) + np.log(poisson(lam[k]).pmf(d))) > -np.inf:
                if (np.log(p_k[k]) + np.log(poisson(lam[k]).pmf(d))) > max_ll:
                    max_class = k
                    max_ll = (np.log(poisson(lam[k]).pmf(d)))


        clusterToPatient[max_class].add(patient)

    for cluster in clusterToPatient:
        if not clusterToPatient[cluster]:
            clusterToPatient[cluster].add('EMPTY PATIENTS')

    return clusterToPatient



def pmm_cross_validate(num_components, patientToGenes, test_size):
    """
    :param num_components:
    :param patientToGenes:
    :param test_size:
    :return: The average likelihood of the model when applied to a new test set, and its BIC
    """

def pmm_get_likelihood(patientToGenes, lam, p_k):
    D = [len(patientToGenes[p]) for p in patientToGenes]
    numCases = len(D)
    classes = range(len(lam))

    ll_kd = np.array([ [np.log(p_k[k]) + np.log(poisson(lam[k]).pmf(d)) for d in D] for k in classes])
    likelihood_sums = np.zeros(numCases)

    for i in range(numCases):
        likelihood_sums[i] = sum([(np.exp(ll_kd[k][i]) if ll_kd[k][i] > -np.inf else 0) for k in range(num_components)] )

    # complete log likelihood

    ll_new = sum(np.log(np.array([ls for ls in likelihood_sums if ls > 0])))
    
    return ll_new




def partition_pmm(patientToGenes, num_components, diff_thresh=10, num_bins=50, max_iter=100, by_iter=True,
                  rand_init=False, far_rand_init=False, do_plot=False, get_best=True, min_cluster_size=0):


    # get the whole data distribution


    # D = [1,2,3,4,5, 100, 150, 200, 1000]
    D = [len(patientToGenes[p]) for p in patientToGenes]
    numCases = len(D)
    data = collections.OrderedDict()

    # print "D is ", D

    # get the lambdas at equal-spaced intervals


    lam = [np.percentile(D, (i + 1) * 100.0 / (num_components + 1)) for i in range(num_components)]
    p_k = [1.0 / num_components for i in range(num_components)]
    classes = range(num_components)

    if rand_init:
        old_lam = lam
        old_p_k = p_k
        #random sample  in a range centered at the quartiles
        lam = [np.random.uniform(l - 0.5 * old_lam[0], l + 0.5 * old_lam[0]) for l in old_lam]
        rand_freq = [2**np.random.uniform(-1, 1) * pk for pk in old_p_k]
        p_k = list(np.array(rand_freq)/sum(rand_freq))
        classes = range(num_components)

    if far_rand_init:
        lam = [np.random.uniform(min(D), max(D)) for l in lam]
        rand_freq = [np.random.uniform(0, 1) for l in lam]
        p_k = list(np.array(rand_freq)/sum(rand_freq))


    print "Initial Lambda is ", lam
    print "Initial p_k is", p_k

    data['Initial Means'] = np.round(lam,1)
    data['Initial p_k'] = np.round(p_k, 2)

    ll = -3e100
    num_iter = 0

    # stupid inital values
    p_k_d= np.zeros(num_components)
    lam_prev = np.zeros(num_components)
    p_k_prev = np.zeros(num_components)

    # for the best values
    ll_best = -np.inf
    p_k_best = None
    lam_best = None

    lls = []

    while 1:


        # We have the log-likelihood of data d and class k in matrix
        #            data 1 data 2 data 3
        # clsss 1   ll_11   ll_12
        # class 2
        ll_kd = np.array([ [np.log(p_k[k]) + np.log(poisson(lam[k]).pmf(d)) for d in D] for k in classes])


        # Likelihood_sums: the total likelihood of each data, summed across class k
        likelihood_sums = np.zeros(numCases)

        for i in range(numCases):
            likelihood_sums[i] = sum([(np.exp(ll_kd[k][i]) if ll_kd[k][i] > -np.inf else 0) for k in range(num_components)] )

        # complete log likelihood

        ll_new = sum(np.log(np.array([ls for ls in likelihood_sums if ls > 0])))

        if num_iter == 0:
            data['Initial LL'] = np.round(ll_new)

        print "ll_new is ", ll_new

        # if ll_new > -27000:
        #     print "P-K is ", p_k, "Lam is ", lam, "Likelihood is ", ll_new
        #     clusterToPatient = pmm_to_cluster(patientToGenes, classes, lam, p_k)
        #     plot_pmm_clusters(patientToGenes, clusterToPatient, num_components, num_bins=100)
        #     print "Missing patients: ", len(clusterToPatient[-1])


        if ll_new > ll_best:
            ll_best = ll_new
            p_k_best = p_k
            lam_best = lam

        # When we break out of the loop, take previous value since it might have jumped out
        if (by_iter):
            if num_iter > max_iter:
                break
        else:
            if (ll_new - ll < diff_thresh):

                p_k_d = p_k_d_prev
                lam = lam_prev
                p_k = p_k_prev

            break

        p_k_d_prev = p_k_d
        lam_prev = lam
        p_k_prev = p_k


        # Calculate p_k_d. This is p(data d | class k) * p(class k)/sum(p(data|class i) *p(class i);
        # i.e. prob of this class given this data

        p_k_d = np.zeros(ll_kd.shape)

        for i in range(numCases):
            # Use max class likelihood to divide all the likelihoods by
            max_val = np.amax(ll_kd, axis=0)[i]

            # sum the likekhoods for every class, make this the denominator of probability
            denom = sum([(np.exp(ll_kd[k][i] - max_val) if ll_kd[k][i] > -np.inf else 0) for k in range(num_components)])

            for k in range(num_components):
                p_k_d[k][i] = (np.exp(ll_kd[k][i] - max_val) / denom if ll_kd[k][i] > -np.inf else 0)
                # print "numerator is ", np.exp(ll_kd[k][i] - max), " prob is ", p_k_d[k][i]

        # print "p_k_d is ", p_k_d

        # sum probabilities of each data being each class over all data
        Z_k = p_k_d.sum(axis=1)


        # see derivation

        lam = [sum([p_k_d[k][i] * D[i] for i in range(numCases)]) * 1.0 / Z_k[k] for k in classes]
        p_k = Z_k * 1.0 / numCases

        p_k = p_k/p_k.sum()


        # print "New lambda is ", lam
        # print "New p_k is ", p_k


        ll = ll_new

        lls.append(ll)
        num_iter += 1



    if get_best:
        p_k = p_k_best
        lam = lam_best
        ll = ll_best





    clusterToPatient = pmm_to_cluster(patientToGenes, classes, lam, p_k)


    # plot patients in each cluster

    if do_plot:
        plot_pmm_clusters(patientToGenes, clusterToPatient, num_components, num_bins=100)

    print "Missing patients: ", len(clusterToPatient[-1])



    data['Number'] = num_components
    data['Means'] = np.round(lam, 1)
    data['Probabilities'] = np.round(p_k, 2)
    data['Likelihood'] = np.round(ll)
    data['Classes'] = classes
    data['AIC'] = np.round(2 * (len(p_k) + len(lam)) - 2 * ll)
    data['BIC'] = np.round(-2 * ll + (len(p_k) + len(lam)) * np.log(numCases))
    data['Missing'] = len(clusterToPatient[-1]) if -1 in clusterToPatient else 0
    data['MinClusterSize'] = min([len(clusterToPatient[c]) if c != -1 else np.inf  for c in clusterToPatient])
    data['MoreThanMin'] = 1 if data['MinClusterSize'] > min_cluster_size else 0

    return data, lls

def plot_pmm_clusters(patientToGenes, clusterToPatient, num_components, num_bins=100):

    D = [len(patientToGenes[p]) for p in patientToGenes]

    bins = range(0, max(list(D)), max(list(D))/num_bins)
    plt.figure()
    for cluster in clusterToPatient:
        plt.hist([len(patientToGenes[p]) for p in clusterToPatient[cluster]], bins=bins, label=str(cluster), alpha = 1.0/num_components)
    plt.xlabel('# Somatic Mutations In Tumor', fontsize=20)
    plt.ylabel('Number of Samples', fontsize=20)
    plt.legend()
    plt.title("Cluster size " + str(num_components), fontsize=20)
    plt.show()

def plot_likelihoods(ll_record):
    plt.figure()
    for i in range(len(ll_record)):
        plt.plot(ll_record[i], label=str(i))
    plt.title("Log-likelihood change in EM", fontsize=20)
    plt.legend(loc=4)
    plt.show()






def main():
    # INDEX BY LOSSES
    mutationmatrix = '/Users/jlu96/maf/new/OV_broad/OV_broad-cna-jl.m2'
    patientFile = '/Users/jlu96/maf/new/OV_broad/shared_patients.plst'
    out_file = '/Users/jlu96/conte/jlu/Analyses/CancerMutationDistributions/OV_broad-cna-jl-PMM.csv'
    partition_file = '/Users/jlu96/maf/new/OV_broad/OV_broad-cna-jl.ppf'
    min_cluster_size = 30
    num_init = 9
    minComp = 4
    maxComp = 6
    do_plot = True
    do_gmm = False
    do_dna = True
    num_integrated = 4
    do_kmeans = False
    do_pmm = False
    geneFile = None
    minFreq = 0
    dna_gene_file = '/Users/jlu96/conte/jlu/Analyses/CancerGeneAnalysis/DNADamageRepair_loss.txt'


    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)

    p_gene_list = []

    with open(dna_gene_file, 'rU') as row_file:
        reader = csv.reader(row_file, delimiter='\t')
        for row in reader:
            p_gene_list.append(row[0])


    if do_kmeans:
        datas = []
        for i in np.arange(minComp, maxComp, 1):
            datas.append(partition_gene_kmeans(geneToCases, patientToGenes, p_gene_list, i, num_bins=50, title=None, do_plot=True))

        with open(out_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=datas[0].keys())
            writer.writeheader()
            for row in datas:
                writer.writerow(row)


    if do_dna:
        cohort_dict = partition_gene_list(patientToGenes, p_gene_list, binary=not bool(num_integrated))
        # Make new cohorts over this
        if num_integrated:
            cohort_dict = integrate_cohorts(cohort_dict, numCases, num_integrated)


        cohort_pairings = [(key, cohort_dict[key]) for key in cohort_dict]
        draw_partitions_cohorts(geneToCases, patientToGenes, cohort_pairings, title='DNADamageGenes',
                        num_bins=100 if mutationmatrix[-9:] == 'cna-jl.m2' else 50)


    if do_gmm:
        datas = []
        for i in np.arange(minComp, maxComp, 1):
            datas.append(partition_GMM(patientToGenes, i, num_bins=50, title='GMM size ' + str(i), do_plot=do_plot))

        with open(out_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=datas[0].keys())
            writer.writeheader()
            for row in datas:
                writer.writerow(row)


    if do_pmm:
        datas = []
        clusters = []
        for num_components in np.arange(minComp, maxComp, 1):
            best_data, clusterToPatient = best_pmm(patientToGenes, num_components, rand_num=5, far_rand_num=5,
                                                   min_cluster_size=min_cluster_size)
            datas.append(best_data)
            clusters.append(clusterToPatient)
            # data, lls = partition_pmm(patientToGenes, i, num_bins=50, max_iter=20, rand_init=False, do_plot=True)
            # datas.append(data)
            # all_lls.append(lls)
            # for j in range(num_init):
            #     data, lls = partition_pmm(patientToGenes, i, num_bins=50, max_iter=20)
            #     datas.append(data)
            #     all_lls.append(lls)

        # os.system('say "Jonathan your program has finished"')


        # get the best BIC
        combined = zip(datas, clusters)
        combined = sorted(combined, key=lambda entry: (entry['MoreThanMin'], entry['BIC']))
        datas, clusters = zip(*combined)


        with open(out_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=datas[0].keys())
            print datas
            writer.writeheader()
            for row in datas:
                writer.writerow(row)


        best_data = datas[0]
        clusterToPatient = clusters[0]

        # code to parition by best clusters
        with open(partition_file, 'w') as csvfile:
            writer = csv.writer(csvfile)

            writer.writerow(['Likelihood', best_data['Likelihood']])
            writer.writerow(['BIC', best_data['BIC']])
            writer.writerow(['NumComponents', best_data['Number']])
            writer.writerow(['Cluster', 'Lambda', 'Probability', 'Patients'])
            for k in clusterToPatient:
                if k != -1:
                    lam = best_data['Means'][k]
                    p_k = best_data['Probabilities'][k]
                else:
                    lam = None
                    p_k = None
                writer.writerow([k, lam, p_k] + list(clusterToPatient[k]))

        load_patient_cohorts(partition_file)



# If there are any patients that aren't assigned, i.e. in cluster -1
# Throw them out?
def load_patient_cohorts(partitionfile, patientToGenes, add_to_closest=True):
    clusterToPatient = {}

    with open(partitionfile, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if (row[0] == 'Cluster'): break
        for row in reader:
            clusterToPatient[row[0]] = {}
            clusterToPatient[row[0]]['Mean'] = row[1]
            clusterToPatient[row[0]]['Probability'] = row[2]
            clusterToPatient[row[0]]['Patients'] = set(row[3:])


    if -1 in clusterToPatient:
        if add_to_closest:
            other_cs = clusterToPatient.keys()
            other_cs.remove(-1)
            for patient in clusterToPatient[-1]:
                sims = [(abs(len(patientToGenes[patient]) - clusterToPatient[c]['Mean']), c) for c in other_cs]
                sims = sorted(sims, key = lambda entry: entry[0])
                best_c = sims[0][1]
                clusterToPatient[best_c].add(patient)

        clusterToPatient.pop(-1)

    cohort_dict = {}

    for c in clusterToPatient:
        cohort_dict[c] = clusterToPatient[c]['Patients']

    return cohort_dict


if __name__ == '__main__':
    main()