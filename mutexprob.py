__author__ = 'jlu96'

import scipy.stats as stats
import numpy as np
import time
import matplotlib.pyplot as plt
import bisect

class Memoize:
    def __init__(self, f):
        self.f = f
        self.memo = {}
    def __call__(self, *args):
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]

class TableLookup:
    def __init__(self, f):
        self.f = f

        #A dictionary of mut_average to overlap,
        #WHere mut_average is the average of the mutation frequencies,
        #And overlap is the amount of overlap needed to reach significance
        self.n = 0
        self.p = 0.0
        self.trials = 0
        self.p_buffer = 0.0
        self.set_size = 0
        self.interval = 0
        self.max_percentage = 0
        self.min_frequencies = []
        self.necessary_overlaps = []


    def __call__(self, *args, **kwargs):
        # # First check the table of computed overlaps necessary for significance.
        # # If the number of overlaps you have is less than the number required,
        # # Then simply return False.
        #
        #Unzip the commands to issigoverlap_multi
        n, mut_nums, p, overlap, trials = args

        #
        # p_buffer = kwargs['p_buffer'] if ('p_buffer' in kwargs) else 0.01
        # precision = kwargs['precision'] if ('precision' in kwargs) else 0.5
        # interval = kwargs['interval'] if ('interval' in kwargs) else 10
        # max_frequency = kwargs['max_frequency'] if ('max_frequency' in kwargs) else 0.5 * n
        #
        # #In case the function has been called with other parameters, reset the table.
        # if (self.n != n or self.p != p or self.set_size != len(mut_nums) or self.trials != trials):
        #     self.n = n
        #     self.p = p
        #     self.trials = trials
        #     self.p_buffer = p_buffer
        #     self.set_size = len(mut_nums)
        #     self.interval = interval
        #     self.max_frequency = max_frequency
        #
        #     #set the minimum frequencies and table
        #
        #     self.min_frequencies = np.round(np.arange(1, max_frequency, interval))
        #     self.necessary_overlaps = np.zeros(len(self.min_frequencies))
        #
        #     p_lax = p + p_buffer
        #
        #     print "Precomputing thresholds for significant overlap for ", len(self.min_frequencies), "frequencies"
        #     t = time.time()
        #     for i in range(len(self.min_frequencies)):
        #         self.necessary_overlaps[i] = precomputesigoverlap(n, trials, p_lax, precision,
        #                                                    *[self.min_frequencies[i] for j in range(self.set_size)])
        #
        #     print "Total time needed to precompute thresholds", time.time() - t
        #
        # # Check for the lowest mutation number in the set, and find the overlap needed for significance
        # min_mut_num = min(mut_nums)
        # index = bisect.bisect_left(self.min_frequencies, min_mut_num)
        # if overlap < self.necessary_overlaps[index]:
        #     return False
        #
        # else:
        #     return cooccurprob_approximate(n, overlap, trials, *mut_nums) <= p

        if overlap < min(mut_nums) * 0.1:
            return False
        else:
            return cooccurprob_approximate(n, overlap, trials, *mut_nums) <= p


@TableLookup
def issigcooccur_multi(n, mut_nums, p, overlap, trials):
    return


def issigmutex(n, mut_nums, p, overlap = 0):
    """Determines the statistical significance of the mutually exclusive mutations,
    which occur at mut_nums frequencies in the n samples.
    **n** (*int*) - number of samples
    **mut_nums** (*list*) - list of each mutation's frequency in the samples
    **p** (*float*) - p-value
    """
    surprise = mutexprob(n, mut_nums, overlap)
    return surprise < p

def issigcooccur(n, mut_nums, p, overlap = 0):
    """Determines the statistical significance of the mutually exclusive mutations,
    which occur at mut_nums frequencies in the n samples.
    **n** (*int*) - number of samples
    **mut_nums** (*list*) - list of each mutation's frequency in the samples
    **p** (*float*) - p-value
    """
    surprise = cooccurprob(n, mut_nums, overlap)
    return surprise < p

def mutexprob(n, mut_nums, overlap = 0):
    """Calculates the tail probability of mutually exclusive
    mutations occurring in at mut_nums frequencies in the
    n samples."""

    if (len(mut_nums) == 2): return mutexprob2(n, mut_nums[0], mut_nums[1], overlap)

def cooccurprob(n, mut_nums, overlap = 0):
    """Calculates the tail probability of cooccurring
    mutations occurring in at mut_nums frequencies in the
    n samples."""
    if (len(mut_nums) == 2): return cooccurprob2(n, mut_nums[0], mut_nums[1], overlap)

@Memoize
def mutexprob2(n, a, b, overlap):
    oddsratio, pvalue = stats.fisher_exact([[overlap, a - overlap], [b - overlap, n - a - b + overlap]],
                              alternative= 'less')
    return pvalue * 0.5


@Memoize
def cooccurprob2(n, a, b, overlap):
    oddsratio, pvalue = stats.fisher_exact([[overlap, a - overlap], [b - overlap, n - a - b + overlap]],
                              alternative= 'greater')
    return pvalue * 0.5
    # p = 1.0
    # prev_m = 0
    # mut_nums = [a, b]
    # for m in mut_nums:
    #     p *= ratio(n, prev_m, m)
    #     prev_m += m
    # return p

@Memoize
def prob3(n, a, b, c):
    p = 1.0
    prev_m = 0
    mut_nums = [a, b, c]
    for m in mut_nums:
        p *= ratio(n, prev_m, m)
        prev_m += m
    return p

@Memoize
def ratio(n, a, b):
    """Calculates (n-a_C_b)/(n_C_b): the number of ways of mutating samples not already mutated,
    divided by the number of ways of mutating samples, period."""
    i = 1.0
    for c in range(b):
        i *= (n-a-c)
        i /= (n-c)
    return i

@Memoize
def p_exclusive(n, a, b, pA, pB):
    return ratio(n, a, b) * stats.binom.pmf(a, n, pA) * stats.binom.pmf(b, n, pB)

@Memoize
def p_exclusive_total(n, pA, pB):
    a_range = range(n + 1)
    b_range = range(n + 1)
    p_value = 0.0
    for a in a_range:
        for b in b_range:
            p_value += p_exclusive(n, a, b, pA, pB)
    return p_value





# class Memoize_CooccurProb:
#     def __init__(self, f):
#         self.f = f
#         self.memo = {}
#     def __call__(self, *args, **kwargs):
#         keyword_tuple = tuple(sorted(kwargs.items()))
#         if not (args, keyword_tuple) in self.memo:
#             n = args[0]
#             overlap = args[1]
#             trials = kwargs['trials']
#             self.memo[(args, keyword_tuple)] = self.f(n, overlap, trials, *args[2:])
#         return self.memo[(args, keyword_tuple)]

#@Memoize_CooccurProb

@Memoize
def cooccurprob_approximate(n, overlap, trials, *mut_nums):
    t = time.time()
    overlap_limit = len(mut_nums) - 1
    cooccur_num = 0
    for i in range(trials):
        # Create a random mutation matrix and sum across columns.
        row = np.zeros(n)
        numOverlaps = 0
        for m in mut_nums:
            random_row = np.concatenate([np.ones(m), np.zeros(n - m)])
            np.random.shuffle(random_row)
            row += random_row


        if (len([entry for entry in row if entry > overlap_limit]) >= overlap): cooccur_num += 1

    print "Time used", time.time() - t, 'for mut_nums', mut_nums, ' and overlap ', overlap, ' p equals ', cooccur_num * 1.0/trials

    return cooccur_num * 1.0/trials

#This needs to be adjusted like below for cooccur
#mut_nums is a key word arg: add the *
#chang
def mutexprob_approximate(n, num_mutex, trials, *mut_nums):
    t = time.time()
    overlap_limit = len(mut_nums) - 1
    mutex_num = 0
    for i in range(trials):
        # Create a random mutation matrix and sum across columns.
        row = np.zeros(n)
        numOverlaps = 0
        for m in mut_nums:
            random_row = np.concatenate([np.ones(m), np.zeros(n - m)])
            np.random.shuffle(random_row)
            row += random_row

        if len([entry for entry in row if entry == 1.0]) >= num_mutex: mutex_num += 1

    print "Time used", time.time() - t, 'for mut_nums', mut_nums, ' and number of mutually exclusive alterations ', \
        num_mutex, ' p equals ', mutex_num * 1.0/trials

    return mutex_num * 1.0/trials

    # t = time.time()
    # print "Number of trials", trials
    # mutex_num = 0
    # for i in range(trials):
    #     # Create a random mutation matrix and sum across columns.
    #     row = np.zeros(n)
    #     numOverlaps = overlap
    #     for m in mut_nums:
    #         row += np.random.permutation(np.concatenate([np.ones(m), np.zeros(n - m)]))
    #
    #     #Check if it's mutually exclusive
    #
    #     for entry in row:
    #         if (entry > 1.1): #in case of floating point error
    #             numOverlaps -= 1
    #             if (numOverlaps < 0):
    #                 break
    #     if (numOverlaps >= 0): mutex_num += 1
    # print "Time used", time.time() - t
    #
    # return mutex_num * 1.0/trials

class PrecomputeSigOverlap:
    def __init__(self, f):
        self.f = f
        self.n = 0
        self.p = 0.0
        self.trials = 0
        self.set_size = 0
        self.min_frequencies = []
        self.necessary_overlaps = []

    def __call__(self, *args):
        n, trials, p, precision, mut_nums = args

        #In case the function has been called with other parameters, reset the table.
        if (self.n != n or self.p != p or self.set_size != len(mut_nums) or self.trials != trials):
            self.n = n
            self.p = p
            self.trials = trials
            self.set_size = len(mut_nums)
            self.min_frequencies = []
            self.necessary_overlaps = []

        # This better be equal among all of the indices...
        min_mut_num = min(mut_nums)
        index = bisect.bisect_left(self.min_frequencies, min_mut_num)
        min_overlap = self.necessary_overlaps[index]

        self.min_frequencies.append(min_mut_num)
        new_necessary_overlap = self.f(*args, min_overlap=min_overlap)
        self.necessary_overlaps.append(new_necessary_overlap)

        return new_necessary_overlap

def precomputesigoverlap(n, trials, p, precision, *mut_nums, **kwargs):
    # For a given number of samples n and mutation numbers mut_nums,
    # returns the minimum amount of overlap for it to be significant.

    #Return -1 if can't get significance..

    min_overlap = kwargs['min_overlap'] if ('min_overlap' in kwargs) else 1

    #Check all overlaps up to the minimum.
    overlap_range_length = int(round((min(mut_nums) - min_overlap) * precision))
    overlap_range = np.round(np.linspace(min_overlap, min(mut_nums), overlap_range_length))

    for overlap in overlap_range:
        pvalue = cooccurprob_approximate(n,  overlap, trials, *mut_nums)

        if pvalue <= p:
            return overlap

    return -1



def getsigoverlap(n, trials, p, precision, *mut_nums, **kwargs):
    # For a given number of samples n and mutation numbers mut_nums,
    # returns the minimum amount of overlap for it to be significant.

    #Return -1 if can't get significance..

    min_overlap = kwargs['min_overlap'] if ('min_overlap' in kwargs) else 1

    #Check all overlaps up to the minimum.
    overlap_range_length = int(round((min(mut_nums) - min_overlap) * precision))
    overlap_range = np.round(np.linspace(min_overlap, min(mut_nums), overlap_range_length))

    for overlap in overlap_range:
        pvalue = cooccurprob_approximate(n,  overlap, trials, *mut_nums)

        if pvalue <= p:
            return overlap

    return -1



def plotmutnumsvsoverlap(n, constants=[1, 1, 1], p=0.05, mut_precision=0.1, overlap_precision=0.5, minFraction = 0.1,
                         maxFraction = 0.5, trials=10000):
    t = time.time()
    global mut_range
    global overlaps

    length = int(round(mut_precision * n))

    constants = np.array(constants)
    constants /= (sum(constants) / len(constants)) #normalize the constants

    mut_range = np.round(np.linspace(n * minFraction/min(constants), n * maxFraction/max(constants), length))
    overlaps = np.zeros(length)
    numGenes = 3

    for i in range(length):
        m = mut_range[i]
        mut_nums = [round(m * constants[g]) for g in range(numGenes)]
        # overlaps[i] = getsigoverlap(n, mut_nums, p, overlap_precision)
        overlaps[i] = precomputesigoverlap(n, trials, p, overlap_precision, *mut_nums)

    plt.figure()
    plt.title('Alteration Number versus Overlap Needed for Significance. n = ' + str(n) + ', p = ' + str(p))
    plt.xlabel("Average Alteration Number. mut_nums = " + str(constants))
    plt.ylabel("Overlap for significance")
    plt.plot(mut_range, overlaps)
    plt.show()

    # plt.figure()
    # plt.title('Alteration Number versus Cooccurrence Ratio Needed for Significance. n = ' + str(n) + ', p = ' + str(p))
    # plt.xlabel("Average Alteration Number. mut_nums = " + str(constants))
    # plt.ylabel("Ratio for significance")
    # plt.plot(mut_range, sigratios)
    # plt.show()

    print "Time used: ", time.time() - t, " for n = ", n, ", length = ", length, \
        "mutation precision = ", mut_precision, ", overlap precision = ", overlap_precision
    return mut_range, overlaps



def multiplotmutnumsvsoverlap(n,  constants_list=[[1, 1, 1], [1, 2, 2]], p=0.05, mut_precision=0.1, overlap_precision=0.5,
                              minFraction = 0.1, maxFraction = 0.5):

    mut_range_list = []
    overlap_list = []
    for constants in constants_list:
        mut_range, overlap = plotmutnumsvsoverlap(n, constants, p, mut_precision, overlap_precision, minFraction, maxFraction)
        mut_range_list.append(mut_range)
        overlap_list.append(overlap)

    return mut_range_list, overlap_list


def ratiopvalue(n, mut_nums):
    length = 20
    overlap_range = np.linspace(0, min(mut_nums) * 2 / 3, length)
    ratios = np.zeros(length)
    pvalues = np.zeros(length)
    for i in range(length):
        overlap = overlap_range[i]
        ratios[i] = cooccurrence_ratio(n, mut_nums, overlap)
        pvalues[i] = cooccurprob_approximate(n, mut_nums, overlap)
    return ratios, pvalues

def plotratiovspvalue(n, mut_nums):
    ratios, pvalues = ratiopvalue(n, mut_nums)
    plt.figure()
    plt.plot(ratios, pvalues, 'bo')
    plt.xlabel('Cooccurrence Ratio')
    plt.ylabel('pvalue')
    plt.title('Cooccurrence Ratio vs pvalue for n = ' + str(n) + ' and mutation frequencies' + str(mut_nums))
    plt.show()


#
# def random_matrix_row(n, mut_nums):
#     """Returns an array of sums across mutations of a random mutation matrix;
#     i.e. the row of the total number of mutations each sample has"""
#     row = np.zeros(n)
#     for m in mut_nums:
#         row += np.random.permutation(np.concatenate([np.ones(m), np.zeros(n - m)]))
#     return row
#
# def matrix_row_mutex(row, overlap = 0):
#     for entry in row:
#         if (entry > 1.1): #in case of floating point error
#             overlap -= 1
#             if (overlap < 0):
#                 return False
#     return True

#
# def random_mut_matrix(n, mut_nums):
#     matrix = []
#     for m in mut_nums:
#         random_row = np.random.permutation(np.concatenate([np.ones(m), np.zeros(n - m)]))
#         matrix.append(random_row)
#     return matrix

# def matrix_mutex(matrix, overlap = 0):
#     column_sum = sum(matrix)
#     for entry in column_sum:
#         if (entry > 1.1): #in case of floating point error
#             overlap -= 1
#             if (overlap < 0):
#                 return False
#     return True

# def mutexprob_approximate(n, mut_nums, overlap = 0, trials = 10000):
#     mutex_num = 0
#     for i in range(trials):
#         row = random_matrix_row(n, mut_nums)
#         mutex_num += matrix_row_mutex(row, overlap)
#     return mutex_num * 1.0/trials