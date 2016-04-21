__author__ = 'jlu96'
import mutex as mex
import matplotlib.pyplot as plt
import csv
import numpy as np
from sklearn import mixture
from sklearn.cluster import KMeans
from sklearn.cross_validation import KFold
from scipy.stats import poisson
from scipy import stats
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


def get_patients_gene_mut_num(patients, genes, patientToGenes):
    return [set.intersection(patientToGenes[p], genes) for p in patients]

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

def integrate_cohorts_sizes(cohort_dict, sizes):
    cohorts_int = {}
    size_index = 0
    num_in_cohort = 0
    new_cohort = []
    for i in cohort_dict.keys():
        num_in_cohort += len(cohort_dict[i])
        new_cohort.extend(cohort_dict[i])
        if (num_in_cohort > sizes[size_index]):
            cohorts_int[size_index] = new_cohort
            size_index += 1
            new_cohort = []
            num_in_cohort = 0

    if new_cohort:
        cohorts_int[size_index] = new_cohort

    return cohorts_int


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


def get_cluster_gTC_pTG(geneToCases, patientToGenes, patients):
    new_pTG = dict([c for c in patientToGenes.items() if c[0] in patients])
    new_genes = set.union(*new_pTG.values())
    new_gTC = dict([g for g in geneToCases.items() if g[0] in new_genes])
    for g in new_gTC:
        new_gTC[g] = new_gTC[g].intersection(patients)
    
    for g in new_genes:
        if g in new_gTC and not new_gTC[g]:
            new_gTC.pop(g)
    
    new_genes = new_genes.intersection(set(new_gTC.keys()))
    
    return new_genes, new_gTC, new_pTG










# 3/12/16-Jlu


class PMM:

    def __init__(self, filename=None, delimiter='\t', lam=None, p_k=None, classes=None, patientToGenes=None,
                data = None, clusterToPatient = None, do_fit=True):

        if filename:
            with open(filename, 'rU') as csvfile:
                reader = csv.DictReader(csvfile, delimiter=delimiter)
                row = reader.next()
                print row
                self.lam = eval(row['Means'])
                self.p_k = eval(row['Probabilities'])
                self.classes = eval(row['Classes']) if 'Classes' in row else range(len(self.lam))
                self.num_components = len(self.classes)
        else:
            self.lam = lam
            self.p_k = p_k
            self.classes = classes
            if not classes:
                self.classes = range(len(self.lam))
            self.num_components = len(self.classes)


        self.data = data
        self.clusterToPatient = clusterToPatient
        print "Class is ", self.classes, "Keys are ", self.clusterToPatient.keys()

        self.patientToGenes = patientToGenes

        if patientToGenes and do_fit:
            self.fit_to_data(patientToGenes)

    def fit_to_data(self, patientToGenes, min_cluster_size=0):
        self.patientToGenes = patientToGenes
        self.data, self.clusterToPatient = pmm_fit_to_data(patientToGenes, classes=self.classes, lam=self.lam, p_k=self.p_k,
                                                           min_cluster_size=min_cluster_size)
        return self.data, self.clusterToPatient


    def plot_clusters(self, title):
        plot_pmm_clusters(self.patientToGenes, self.clusterToPatient, self.num_components, title=title)


    def write_clusters(self, partition_file):
        with open(partition_file, 'w') as csvfile:
            writer = csv.writer(csvfile)

            writer.writerow(['Likelihood', self.data['Likelihood']])
            writer.writerow(['BIC', self.data['BIC']])
            writer.writerow(['NumComponents', self.data['Number']])
            writer.writerow(['Cluster', 'Lambda', 'Probability', 'Patients'])
            for k in self.clusterToPatient:
                if k != -1:
                    lam = self.data['Means'][k]
                    p_k = self.data['Probabilities'][k]
                else:
                    lam = None
                    p_k = None
                writer.writerow([k, lam, p_k] + list(self.clusterToPatient[k]))

    def compare_dna(self, dna_cohort_dict, do_KS=False):

        partition_stats_list = []

        sizes = [len(self.clusterToPatient[c]) for c in self.clusterToPatient]

        # partition by genes
        dna_cohorts = integrate_cohorts_sizes(dna_cohort_dict, sizes)

        pmm_cluster_list = []
        dna_cluster_list = []
        
        print "In partition stats Class is ", self.classes, "Keys are ", self.clusterToPatient.keys()
        
        for i in range(len(self.classes)):
            partition_stats = collections.OrderedDict()
            partition_stats['Class'] = self.classes[i]
            partition_stats['Mean'] = self.lam[i]
            partition_stats['Probability'] = self.p_k[i]


            partition_stats['PMM_patients'] = self.clusterToPatient[self.classes[i]]
            partition_stats['DNA_patients'] = dna_cohorts[i]

            pmm_cluster_list.append(partition_stats['PMM_patients'])
            dna_cluster_list.append(partition_stats['DNA_patients'])
            
            dna_pmn = [len(self.patientToGenes[p]) for p in partition_stats['DNA_patients']]
            pmm_pmn = [len(self.patientToGenes[p]) for p in partition_stats['PMM_patients']]

            if do_KS:
                poisson_cdf.mu = self.lam[i]
                partition_stats['KS'] = stats.kstest(dna_pmn, poisson_cdf)

            #qq plot of the dna and then the poisson
            poisson_q = get_quantiles(dna_pmn, pmm_pmn)
            dna_q = get_quantiles(dna_pmn, dna_pmn)

            plot_pmm_clusters(self.patientToGenes, {'PMM': partition_stats['PMM_patients'], 'DNA': partition_stats['DNA_patients'] },
                              2, num_bins=100, title='DNA VS PMN')

            plt.figure()
            plt.plot(dna_q, poisson_q, 'bo')
            plt.plot([0, 100], [0,100], 'r-', label = 'y=x')
            plt.title('QQ for ' + str(self.classes[i]), fontsize=20)
            plt.xlabel('DNA_Q', fontsize=20)
            plt.ylabel('PMM_Q', fontsize=20)
            plt.legend()
            plt.show()

            partition_stats_list.append(partition_stats)

        if do_KS:
            self.data['KS_geom_mean'] = mex.prod([partition_stats['KS'][1] for partition_stats in partition_stats_list]) ** (1.0/ len(partition_stats_list))

            print "KS average is ", self.data['KS_geom_mean']
            
        self.data['CohenKappa'] = cohen_kappa(pmm_cluster_list, dna_cluster_list)


        return partition_stats_list



def cohen_kappa(cluster_list_1, cluster_list_2):
    # assume same categories each
    num_agree = 0
    prob_agree = 0
    total = len(set.union(*[set(c) for c in cluster_list_1]))
    
    num_classes = len(cluster_list_1)
    
    cluster_list_1 = [set(c) for c in cluster_list_1]
    cluster_list_2 = [set(c) for c in cluster_list_2]
    
    for k in range(num_classes):
        a = cluster_list_1[k]
        b = cluster_list_2[k]
        num_agree += len(a.intersection(b))
        prob_agree += (len(a) * len(b) * 1.0) / (total ** 2)
    

    obs_agree = num_agree * 1.0 / total
    
    ck = (obs_agree - prob_agree)/(1.0 - prob_agree)
    
    print "Number agreements ", num_agree
    print "Total ", total
    print "Prob agreements ", prob_agree
    print "Cohen kappa ", ck
    
    return ck
        
    
    




def poisson_cdf(x):
    if not hasattr(poisson_cdf, 'mu'):
        poisson_cdf.mu = 0
    print "X is ", x, "and mu is ", poisson_cdf.mu
    return poisson.cdf(x, poisson_cdf.mu)

def get_quantiles(test_dist, base_dist):
    return [stats.percentileofscore(base_dist, t) for t in test_dist]

def assign_missing(clusterToPatient, patientToGenes):
    if -1 not in clusterToPatient:
        print "No missing patients in clusters"
        return clusterToPatient
    missing_patients = clusterToPatient[-1]
    cluster_means = [(sum([len(patientToGenes[p]) for p in clusterToPatient[c]]) * 1.0 /len(clusterToPatient[c]), c) for c in clusterToPatient if c != -1]
    print cluster_means, cluster_means[0][0]
    for patient in missing_patients:
        num = len(patientToGenes[patient])
        correct_cluster = sorted(cluster_means, key=lambda entry: abs(num - entry[0]))[0][1]
        clusterToPatient[correct_cluster].add(patient)
    clusterToPatient.pop(-1)

    return clusterToPatient



def best_pmm(patientToGenes, num_components, max_iter=30, rand_num=5, far_rand_num=5, min_cluster_size=0,
             plot_clusters=True):

    data_record = []
    lls_record = []

    # Do normal
    first_data, lls = partition_pmm(patientToGenes, num_components,  max_iter=max_iter, min_cluster_size=min_cluster_size)

    data_record.append(first_data)
    lls_record.append(lls)

    # Do best rand init
    for i in range(rand_num):
        data, lls = partition_pmm(patientToGenes, num_components, rand_init=True, max_iter=max_iter, min_cluster_size=min_cluster_size,
                                 verbose=False)
        data_record.append(data)
        lls_record.append(lls)

    for i in range(far_rand_num):
        data, lls = partition_pmm(patientToGenes, num_components, far_rand_init=True, max_iter=max_iter, min_cluster_size=min_cluster_size,
                                 verbose=False)
        data_record.append(data)
        lls_record.append(lls)

    combined_record = zip(data_record, lls_record)

    combined_record = sorted(combined_record, key=lambda entry: (-1 * entry[0]['Missing'], entry[0]['Likelihood']), reverse=True)

    data_record, lls_record = zip(*combined_record)

    best_data = data_record[0]

    if (best_data['Likelihood'] > first_data['Likelihood'] + 10):
        print "First data not best!"
        best_data['IsFirst'] = False
    else:
        best_data['IsFirst'] = True


    clusterToPatient = pmm_to_cluster(patientToGenes, best_data['Classes'], best_data['Means'], best_data['Probabilities'])

    if plot_clusters:
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

    missing_clusters = set()
    for cluster in clusterToPatient:
        if not clusterToPatient[cluster]:
            print '**********NO PATIENTS IN CLUSTER ', lam[cluster], p_k[cluster]
            missing_clusters.add(cluster)
            #clusterToPatient[cluster].add('NO PATIENTS IN CLUSTER')
    for cluster in missing_clusters:
        clusterToPatient.pop(cluster)
            
    return clusterToPatient



def pmm_cross_validate(num_components, patientToGenes, num_folds, kf_random_state=None, max_iter=30, rand_num=5, far_rand_num=5, min_cluster_size=0):
    """
    :return: The average likelihood of the model when applied to a new test set, and its BIC
    """

    kf = KFold(len(patientToGenes), n_folds=num_folds, random_state=kf_random_state)

    lls = []
    missing_patients = []
    bics = []
    for train_index, test_index in kf:

        train_patientToGenes = dict([patientToGenes.items()[x] for x in train_index])
        test_patientToGenes = dict([patientToGenes.items()[x] for x in test_index])
        best_data, _ = best_pmm(train_patientToGenes, num_components, max_iter=max_iter, rand_num=rand_num,
                                               far_rand_num=far_rand_num, min_cluster_size=min_cluster_size)

        test_stats, test_cluster = pmm_fit_to_data(test_patientToGenes, best_data['Classes'], best_data['Means'], best_data['Probabilities'])

        plot_pmm_clusters(test_patientToGenes, test_cluster, num_components, title='Test clusters size ' + str(num_components))

        lls.append(test_stats['Likelihood'])
        missing_patients.append(test_stats['Missing'])
        bics.append(test_stats['BIC'])

    return sum(lls) * 1.0/len(lls), sum(missing_patients) * 1.0 / len(missing_patients), sum(bics) * 1.0/ len(bics)





def pmm_fit_to_data(patientToGenes, classes, lam, p_k, data=None, min_cluster_size=0):
    """
    :param patientToGenes:
    :param lam:
    :param p_k:
    :param data:
    :return: data, clusterToPatient
    """

    if not data:
        data = collections.OrderedDict()


    D = [len(patientToGenes[p]) for p in patientToGenes]
    numCases = len(D)
    num_components = len(lam)

    ll_kd = np.array([ [np.log(p_k[k]) + np.log(poisson(lam[k]).pmf(d)) for d in D] for k in classes])
    likelihood_sums = np.zeros(numCases)

    for i in range(numCases):
        likelihood_sums[i] = sum([(np.exp(ll_kd[k][i]) if ll_kd[k][i] > -np.inf else 0) for k in range(num_components)] )

    # complete log likelihood

    ll = sum(np.log(np.array([ls for ls in likelihood_sums if ls > 0])))

    clusterToPatient = pmm_to_cluster(patientToGenes, classes, lam, p_k)

    print "LL:", np.round(ll), "Missing patients: ", len(clusterToPatient[-1]) if -1 in clusterToPatient else 0

    data['Number'] = num_components
    data['OriginalNumber'] = num_components
    mp = zip(*sorted(zip(list(np.round(lam, 1)), list(np.round(p_k, 2))), key = lambda entry: entry[0]))

    data['Means'], data['Probabilities'] =  list(mp[0]), list(mp[1])   
    data['Likelihood'] = np.round(ll)
    data['Classes'] = classes
    data['AIC'] = np.round(2 * (len(p_k) + len(lam)) - 2 * ll)
    data['BIC'] = np.round(-2 * ll + (len(p_k) + len(lam)) * np.log(numCases))
    data['Missing'] = len(clusterToPatient[-1]) if -1 in clusterToPatient else 0
    data['MinClusterSize'] = min([len(clusterToPatient[c]) if c != -1 else np.inf  for c in clusterToPatient])
    data['MoreThanMin'] = 1 if data['MinClusterSize'] > min_cluster_size else 0
    data['Merged'] = False
    data['MergeHistory'] = set()

    return data, clusterToPatient




def partition_pmm(patientToGenes, num_components, diff_thresh=0.01, num_bins=50, max_iter=100, by_iter=True,
                  rand_init=False, far_rand_init=False, do_plot=False, get_best=True, min_cluster_size=0,
                 verbose=True):


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

    if verbose:
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
    missing_best = numCases

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

            
        missing_new = len([x for x in likelihood_sums if x == 0])
        # complete log likelihood

        ll_new = sum(np.log(np.array([ls for ls in likelihood_sums if ls > 0])))

        if num_iter == 0:
            data['Initial LL'] = np.round(ll_new)

        if verbose:
            print "ll_new is ", ll_new, "missing is ", missing_new


        if ll_new > ll_best or missing_new < missing_best:
            ll_best = ll_new
            p_k_best = p_k
            lam_best = lam
            missing_best = missing_new

        # When we break out of the loop, take previous value since it might have jumped out
        if (by_iter):
            if num_iter > max_iter:
                break
            elif abs(ll_new - ll) < diff_thresh:
                break
        else:
            if abs(ll_new - ll) < diff_thresh:

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





    data, clusterToPatient = pmm_fit_to_data(patientToGenes, classes, lam, p_k, data=data, min_cluster_size=min_cluster_size)
    # plot patients in each cluster

    if do_plot:
        plot_pmm_clusters(patientToGenes, clusterToPatient, num_components, num_bins=100)


    # clusterToPatient = pmm_to_cluster(patientToGenes, classes, lam, p_k)

    #
    #
    #
    #
    # data['Number'] = num_components
    # data['Means'] = np.round(lam, 1)
    # data['Probabilities'] = np.round(p_k, 2)
    # data['Likelihood'] = np.round(ll)
    # data['Classes'] = classes
    # data['AIC'] = np.round(2 * (len(p_k) + len(lam)) - 2 * ll)
    # data['BIC'] = np.round(-2 * ll + (len(p_k) + len(lam)) * np.log(numCases))
    # data['Missing'] = len(clusterToPatient[-1]) if -1 in clusterToPatient else 0
    # data['MinClusterSize'] = min([len(clusterToPatient[c]) if c != -1 else np.inf  for c in clusterToPatient])
    # data['MoreThanMin'] = 1 if data['MinClusterSize'] > min_cluster_size else 0

    return data, lls



def sort_data_by_means(data):
    """ Sort in ascending order. Don't need to change cluster labels"""
    data_items = data.items()
    mean_indices = ((i, data['Means'][i]) for i in range(len(data['Means'])))
    mean_indices = sorted(mean_indices, key=lambda entry: min(entry[1]) if isinstance(entry[1], list)
                         else entry[1])
    
    conversion_array = [m[0] for m in mean_indices] # this should map to the correct index now. these are new clusters
    
    new_data = collections.OrderedDict()
    
    for key in data:
        value = data[key]
        if isinstance(value, np.ndarray):
            new_value = np.zeros(len(value))
            for i in range(len(conversion_array)):
                new_value[i] = value[conversion_array[i]]
            new_data[key] = new_value
        if isinstance(value, list):
            new_value = [value[conversion_array[i]] for i in range(len(conversion_array))]
            new_data[key] = new_value
            
        else:
            new_data[key] = value
    
    return new_data
    

def merge_clusters(data, clusterToPatient, patientToGenes,
                  missing_limit=0.5, min_cluster_size=30):
    """Merge adjacent clusters. Choosse to merge those clusters that
    are the most similar, as measured by the likelihood of one within
    another.
    missing_limit is the limit on number of patients that can't
    be explained by one cluster. Clusters will be sorted first
    by those who are below the minimum cluster size,
    less missing patients in their merging
    cluster, then by those that have the highest likelihood
    """
    # get the likelihood of each cluster rel. to other ones
    # only look at adjacent clusters! sort them
    
    data = sort_data_by_means(data)
    
    print "****************************************"
    print "Begin merging."
    # first go forward

    
    classes = data['Classes']
    p_k = data['Probabilities']
    lam = data['Means']
    
    
    all_list = []
    
    for i in range(len(lam) - 1):
        from_index, to_index = i, i + 1
        from_class, to_class = classes[from_index], classes[to_index]
        patients = clusterToPatient[from_class]
        p = [len(patientToGenes[patient]) for patient in patients]
        
        #check if we're dealing with merged clusters. if so... add the likelihoods of the individual
        # underlying poissons?
        if isinstance(p_k[from_index], list):
            clust_probs = p_k[from_index]
            clust_means = lam[from_index]
            clust_size = len(clust_means)
            
            from_ll = [max([np.log(clust_probs[x]) + 
                           np.log(poisson(clust_means[x]).pmf(d)) for x in range(clust_size)])
                          for d in p]
        else:
            from_ll = [np.log(p_k[from_index]) + np.log(poisson(lam[from_index]).pmf(d)) for d in p]
            
        if isinstance(p_k[to_index], list):
            clust_probs = p_k[to_index]
            clust_means = lam[to_index]
            clust_size = len(clust_means)
            
            to_ll = [max([np.log(clust_probs[x]) + 
                           np.log(poisson(clust_means[x]).pmf(d)) for x in range(clust_size)])
                          for d in p]
        else:
            to_ll = [np.log(p_k[to_index]) + np.log(poisson(lam[to_index]).pmf(d)) for d in p]
            
        missing = np.isinf(from_ll) ^ np.isinf(to_ll)
        
        missing_indices = np.where(missing)[0]
        good_indices = np.where(~missing)[0]
        
        missing_num = len(missing_indices)
        
        ll_diffs = [to_ll[j] - from_ll[j] for j in good_indices]
        
        ll_diffs_total = sum(ll_diffs)
        
        all_list.append([(from_index, to_index), missing_num, ll_diffs_total, missing_num > missing_limit * len(p),
                        len(patients) < min_cluster_size])
        
    # now go backwards
    for i in reversed(range(1, len(lam))):
        from_index, to_index = i, i - 1
        from_class, to_class = classes[from_index], classes[to_index]
        patients = clusterToPatient[from_class]
        p = [len(patientToGenes[patient]) for patient in patients]
        
                #check if we're dealing with merged clusters. if so... add the likelihoods of the individual
        # underlying poissons?
        if isinstance(p_k[from_index], list):
            clust_probs = p_k[from_index]
            clust_means = lam[from_index]
            clust_size = len(clust_means)
            
            from_ll = [max([np.log(clust_probs[x]) + 
                           np.log(poisson(clust_means[x]).pmf(d)) for x in range(clust_size)])
                          for d in p]
        else:
            from_ll = [np.log(p_k[from_index]) + np.log(poisson(lam[from_index]).pmf(d)) for d in p]
            
        if isinstance(p_k[to_index], list):
            clust_probs = p_k[to_index]
            clust_means = lam[to_index]
            clust_size = len(clust_means)
            
            to_ll = [max([np.log(clust_probs[x]) + 
                           np.log(poisson(clust_means[x]).pmf(d)) for x in range(clust_size)])
                          for d in p]
        else:
            to_ll = [np.log(p_k[to_index]) + np.log(poisson(lam[to_index]).pmf(d)) for d in p]
        
        
        missing = np.isinf(from_ll) ^ np.isinf(to_ll)
        
        missing_indices = np.where(missing)[0]
        good_indices = np.where(~missing)[0]
        
        missing_num = len(missing_indices)
        
        ll_diffs = [to_ll[j] - from_ll[j] for j in good_indices]
        
        ll_diffs_total = sum(ll_diffs)
        
        
        all_list.append([(from_index, to_index), missing_num, ll_diffs_total, missing_num < missing_limit * len(p),
                        len(patients) < min_cluster_size])
        
    
    # sort by the cluster that's below the min size, then byminimum missing, then by maximum likelihood ratio
    all_list = sorted(all_list, key=lambda entry: (entry[4], entry[3], entry[2]), reverse=True)
    
    print "Possible merged clusters is ", all_list
    print "Best cluster is ", all_list[0]
    

    (from_index, to_index), missing_num, ll_diffs_total, more_than_missing, cluster_too_small = all_list[0]

    # calculate the new AIC, BIC, make new cluster to patient, make new classes..new means? update probabilities
    
    # Record merge history
    new_data = data
    if 'MergeHistory' not in new_data:
        new_data['MergeHistory'] = set()
    
    new_data['MergeHistory'].add((str([lam[from_index], lam[to_index]]),
                  str([p_k[from_index], p_k[to_index]]),
                  (len(clusterToPatient[classes[from_index]]), len(clusterToPatient[classes[to_index]])),
                  missing_num, ll_diffs_total, ('Num classes befpre', len(classes), ('Cluster too small?', cluster_too_small))))
        
    new_clusterToPatient = clusterToPatient
    moved_patients = new_clusterToPatient[classes[from_index]]
    new_clusterToPatient[classes[to_index]] = new_clusterToPatient[classes[to_index]].union(moved_patients)
    new_clusterToPatient.pop(classes[from_index])

    
    print "MERGING the probs and likelihoods"
    if not isinstance(p_k[from_index], list):
        p_k[from_index] = [p_k[from_index]]
        lam[from_index] = [lam[from_index]]
    if not isinstance(p_k[to_index], list):
        p_k[to_index] = [p_k[to_index]]
        lam[to_index] = [lam[to_index]] 
    p_k[to_index].extend(p_k[from_index])
    lam[to_index].extend(lam[from_index])
    new_data['Probabilities'] = p_k
    new_data['Means'] = lam
    
    
    print "MERGING: HERE ARE OLD VALUES", new_data
    #remove all the old values
    new_data['Merged'] = True
    new_data['Number'] -= 1
    for key in new_data:
        value = new_data[key]
        if isinstance(value, np.ndarray):
            value = list(value)
            value = value[0: from_index] + value[from_index + 1 :]
            value = np.array(value)
            new_data[key] = value
        elif isinstance(value, list):
            value = value[0: from_index] + value[from_index + 1 :]
            new_data[key] = value

    print "New classe:", new_data['Classes'], "VS NEW KEYS", new_clusterToPatient.keys()
            
    # integrate the old patients to the new ones

    
    
    new_data['MinClusterSize'] = min(len(new_clusterToPatient[c]) for c in new_clusterToPatient)
    
    print "MERGING: HERE ARE NEW VALUES", new_data
    
    plot_pmm_clusters(patientToGenes, clusterToPatient, new_data['Number'], title='Merging')
    
    print "End merging."
    print "****************************************"    
    
    return new_data, new_clusterToPatient
 
    
#     data['Number'] = num_components
#     data['Means'], data['Probabilities'] = zip(*sorted(zip(list(np.round(lam, 1)), list(np.round(p_k, 2))), key = lambda entry: entry[0]))
#     data['Likelihood'] = np.round(ll)
#     data['Classes'] = classes
#     data['AIC'] = np.round(2 * (len(p_k) + len(lam)) - 2 * ll)
#     data['BIC'] = np.round(-2 * ll + (len(p_k) + len(lam)) * np.log(numCases))
#     data['Missing'] = len(clusterToPatient[-1]) if -1 in clusterToPatient else 0
#     data['MinClusterSize'] = min([len(clusterToPatient[c]) if c != -1 else np.inf  for c in clusterToPatient])
#     data['MoreThanMin'] = 1 if data['MinClusterSize'] > min_cluster_size else 0

def backward_selection(data, clusterToPatient, patientToGenes, min_cluster_size = 30,
                       max_components = 10):
    """Merge clusters until a criterion is satisfied. Missing patients are assumed to
    be assigned already.
    """
    

    merged_data = data
    merged_cluster = clusterToPatient
    
    while (merged_data['Number'] > max_components or merged_data['MinClusterSize'] < min_cluster_size):
        merged_data, merged_cluster = merge_clusters(merged_data, merged_cluster, patientToGenes,
                                                    min_cluster_size = min_cluster_size)
    
    return merged_data, merged_cluster
    







def plot_pmm_clusters(patientToGenes, clusterToPatient, num_components, num_bins=100, title=None):
    D = [len(patientToGenes[p]) for p in patientToGenes]

    bins = range(0, max(list(D)), max(list(D))/num_bins)
    plt.figure()
    for cluster in clusterToPatient:
        plt.hist([len(patientToGenes[p]) for p in clusterToPatient[cluster]], bins=bins, label=str(cluster), alpha = 1.0/num_components)
    plt.xlabel('# Somatic Mutations In Tumor', fontsize=20)
    plt.ylabel('Number of Samples', fontsize=20)
    plt.legend()
    if not title:
        plt.title("Cluster size " + str(num_components), fontsize=20)
    else:
        plt.title(title, fontsize=20)
    plt.show()

def plot_likelihoods(ll_record):
    plt.figure()
    for i in range(len(ll_record)):
        plt.plot(ll_record[i], label=str(i))
    plt.title("Log-likelihood change in EM", fontsize=20)
    plt.legend(loc=4)
    plt.show()

# If there are any patients that aren't assigned, i.e. in cluster -1
# Throw them out?
def load_patient_cohorts(partitionfile, patientToGenes, add_to_closest=True, delimiter='\t'):
    clusterToProp = {}

    with open(partitionfile, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        for row in reader:
            if (row[0] == 'Cluster'): break
        # reader = csv.DictReader(csvfile, delimiter=delimiter)
        # print "begun dict reader\n"
        for row in reader:
            c = eval(row[0])
            print c
            clusterToProp[c] = {}
            clusterToProp[c]['Mean'] = eval(row[1]) if row[1] else 0
            clusterToProp[c]['Probability'] = eval(row[2]) if row[2] else 0
            clusterToProp[c]['Patients'] = set(row[3:]) if row[3] else set()


    if -1 in clusterToProp:
        if add_to_closest:
            other_cs = clusterToProp.keys()
            other_cs.remove(-1)
            print "Removed ", clusterToProp[-1]
            for patient in clusterToProp[-1]:
                sims = [(abs(len(patientToGenes[patient]) - clusterToProp[c]['Mean']), c) for c in other_cs]
                sims = sorted(sims, key = lambda entry: entry[0])
                best_c = sims[0][1]
                clusterToProp[best_c]['Patients'].add(patient)
            print "completed"

        clusterToProp.pop(-1)

    sorted_clusters = sorted(clusterToProp.keys(), key = lambda entry: clusterToProp[entry]['Mean'])
    
    oldclusterToProp = clusterToProp.copy()
    clusterToProp = {}
    cohort_dict = {}
    
    for i in range(len(sorted_clusters)):
        cohort_dict[i] = oldclusterToProp[sorted_clusters[i]]['Patients']
        clusterToProp[i] = oldclusterToProp[sorted_clusters[i]]
    
    min_cohort = cohort_dict[0]
    
    
    
    

#     for c in clusterToProp:
#         cohort_dict[c] = clusterToProp[c]['Patients']
#     min_cohort = cohort_dict[sorted(clusterToProp.keys(), key=lambda entry: clusterToProp[entry]['Mean'])[0]]

    return cohort_dict, clusterToProp, min_cohort

# INDEX BY LOSSES

def run_partitions(mutationmatrix = None, #'/Users/jlu96/maf/new/OV_broad/OV_broad-cna-jl.m2',
        patientFile = None, #'/Users/jlu96/maf/new/OV_broad/shared_patients.plst',
        out_file = None, #'/Users/jlu96/conte/jlu/Analyses/CancerMutationDistributions/OV_broad-cna-jl-PMM-crossval.txt',
        partition_file = None, #'/Users/jlu96/maf/new/OV_broad/OV_broad-cna-jl.ppf',
        load_pmm_file = None, #'/Users/jlu96/conte/jlu/Analyses/CancerMutationDistributions/OV_broad-cna-jl-PMM.txt',
        dna_pmm_comparison_file = None, #'/Users/jlu96/conte/jlu/Analyses/CancerMutationDistributions/OV_broad-cna-jl-PMM-dnacomp.txt',
        cluster_matrix = None, # '/Users/jlu96/maf/new/OV_broad/OV_broad-cna-jl-cluster.m2',
        min_cluster_size = 15,
        num_init = 9,
        minComp = 2,
        maxComp = 5,
        do_plot = True,
        do_gmm = False,
        do_dna = False,
        num_integrated = 4,
        do_kmeans = False,
        do_pmm = True,
        do_cross_val = False,
        do_pmm_dna = True,
        do_back_selection = True,
        write_cluster_matrices = True,
        rand_num = 3,
        far_rand_num = 3,
        kf_random_state = 1,
        kf_num_folds = 5,

        geneFile = None,
        minFreq = 0,
        dna_gene_file = '/Users/jlu96/conte/jlu/Analyses/CancerGeneAnalysis/DNADamageRepair_loss.txt',
       out_dir = '/Users/jlu96/conte/jlu/Analyses/CancerMutationDistributions/',
        write_all_partitions = True):
    
    mutationmatrix_list = mutationmatrix.split('/')
    matrix_dir = '/'.join(mutationmatrix_list[:-1]) + '/'
    prefix = (mutationmatrix_list[-1]).split('.m2')[0]
    

    if not patientFile:
        patientFile = matrix_dir + 'shared_patients.plst'
        
    if not out_file:
        if do_cross_val:
            out_file = out_dir + prefix + '-PMM-crossval-kf' + str(kf_num_folds) + '.txt'
        else:
            out_file = out_dir + prefix + '-PMM-comparisons.txt'
    
    if not partition_file:
        partition_file = matrix_dir + prefix + '.ppf'
        
    
    if not load_pmm_file:
        load_pmm_file = out_dir + prefix + '-PMM.txt'
    
    if not dna_pmm_comparison_file:
        dna_pmm_comparison_file = out_dir + prefix + '-PMM-dnacomp.txt'
        
    if not cluster_matrix:
        cluster_matrix = matrix_dir + prefix + '-cluster.m2'

    
    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq)

    p_gene_list = []

    with open(dna_gene_file, 'rU') as row_file:
        reader = csv.reader(row_file, delimiter='\t')
        for row in reader:
            p_gene_list.append(row[0])
        dna_cohort_dict = partition_gene_list(patientToGenes, p_gene_list, binary=not bool(num_integrated))


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

        partition_stats_list = []
        for num_components in np.arange(minComp, maxComp, 1):
            best_data, clusterToPatient = best_pmm(patientToGenes, num_components, rand_num=rand_num, far_rand_num=far_rand_num,
                                                   min_cluster_size=min_cluster_size)

            if do_back_selection:
                # assign the missing data
                clusterToPatient = assign_missing(clusterToPatient, patientToGenes)
                best_data, clusterToPatient = backward_selection(best_data, clusterToPatient, patientToGenes, min_cluster_size = min_cluster_size,
                       max_components = maxComp)
            
            if do_pmm_dna:
                print "cfirst lasses are ", best_data['Classes'], "clusterToPatient is ", clusterToPatient.keys()
                pmm = PMM(lam=best_data['Means'], p_k=best_data['Probabilities'], patientToGenes=patientToGenes,
                         data=best_data, clusterToPatient=clusterToPatient, classes=best_data['Classes'],
                          do_fit=False)

                partition_stats_list.extend(pmm.compare_dna(dna_cohort_dict))

                best_data = pmm.data


            if do_cross_val:
            #cross validate each of the components
                print "*******************************************************************************************************"
                print "BEGINNING CROSS VALIDATION for ", num_components
                print "*******************************************************************************************************"
                best_data['TestLL'], best_data['TestMissing'], best_data['TestBIC'] = pmm_cross_validate(num_components, patientToGenes,
                                                                                                         num_folds=kf_num_folds,
                                                                                                     kf_random_state=kf_random_state,
                                                                                   rand_num=rand_num, far_rand_num=far_rand_num,
                                                                                   min_cluster_size=min_cluster_size)
                best_data['TestFolds'] = kf_num_folds

                print "*******************************************************************************************************"
                print "EMDING CROSS VALIDATION  for ", num_components
                print "*******************************************************************************************************"

            datas.append(best_data)
            clusters.append(clusterToPatient)
            
            if write_all_partitions:
                with open(partition_file + str(num_components), 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter='\t')

                    writer.writerow(['Likelihood', best_data['Likelihood']])
                    writer.writerow(['BIC', best_data['BIC']])
                    writer.writerow(['NumComponents', best_data['Number']])
                    writer.writerow(['Cluster', 'Mean', 'Probability', 'Patients'])
                    if 'Merged' in best_data and best_data['Merged']:
                        for k in range(len(clusterToPatient)):
                            lam = best_data['Means'][k]
                            p_k = best_data['Probabilities'][k]
                            writer.writerow([best_data['Classes'][k] , lam, p_k]  + list(clusterToPatient[best_data['Classes'][k]]))
                        
                    else:
                        for k in clusterToPatient:
                            if k != -1:
                                lam = best_data['Means'][k]
                                p_k = best_data['Probabilities'][k]
                            else:
                                lam = None
                                p_k = None
                            writer.writerow([k, lam, p_k] + list(clusterToPatient[k]))

        # get the best BIC
        combined = zip(datas, clusters)
        if do_cross_val:
            combined = sorted(combined, key=lambda entry: ( -1 * entry[0]['MoreThanMin'], np.round(entry[0]['TestMissing']), -1 * entry[0]['TestLL'], entry[0]['TestBIC'], entry[0]['BIC']))
        else:
            combined = sorted(combined, key=lambda entry: ( -1 * entry[0]['MoreThanMin'], entry[0]['BIC']))

        datas, clusters = zip(*combined)




        with open(out_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=datas[-1].keys(), delimiter='\t', extrasaction='ignore')
            print datas
            writer.writeheader()
            for row in datas:
                writer.writerow(row)


        best_data = datas[0]
        clusterToPatient = clusters[0]

        # code to parition by best clusters
        with open(partition_file, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')

            writer.writerow(['Likelihood', best_data['Likelihood']])
            writer.writerow(['BIC', best_data['BIC']])
            writer.writerow(['NumComponents', best_data['Number']])
            writer.writerow(['Cluster', 'Mean', 'Probability', 'Patients'])
            if 'Merged' in best_data and best_data['Merged']:
                for k in range(len(clusterToPatient)):
                    lam = best_data['Means'][k]
                    p_k = best_data['Probabilities'][k]
                    writer.writerow([best_data['Classes'][k] , lam, p_k]  + list(clusterToPatient[best_data['Classes'][k]]))
                        
            else:
                for k in clusterToPatient:
                    if k != -1:
                        lam = best_data['Means'][k]
                        p_k = best_data['Probabilities'][k]
                    else:
                        lam = None
                        p_k = None
                    writer.writerow([k, lam, p_k] + list(clusterToPatient[k]))

        if write_cluster_matrices:
            for cluster in clusterToPatient:
                with open(cluster_matrix + str(cluster), 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter='\t')
                    for patient in clusterToPatient[cluster]:
                        writer.writerow('\t'.join([patient] + list(patientToGenes[patient])))


        if do_pmm_dna:
            with open(dna_pmm_comparison_file, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=partition_stats_list[0].keys(), delimiter='\t')
                writer.writeheader()
                print "header written"
                for row in partition_stats_list:
                    writer.writerow(row)
