__author__ = 'jlu96'

import numpy as np

def BH(results, pvalues, Q=0.05, use_dependencies=False):
    """Benjamini-Hochberg

    Given a list of results and list of their corresponding p-values:
    return those that control the FDR at the given level Q

    returns: results, pvalues
    """
    if len(results) != len(pvalues):
        raise ValueError("Results and p-values must have same length")

    result_pvalue_sorted = sorted(zip(results, pvalues), key=lambda entry: entry[1])
    m = len(pvalues)

    base_threshold = Q * 1.0 / m

    if use_dependencies:
        # change the base threshold to Q * 1.0/m * sum_i=1^m 1/i
        whole_sum = np.sum(1.0/ np.arange(1, m + 1))
        base_threshold /= whole_sum


    # i will be the cutoff index for discoveries
    for i in range(len(result_pvalue_sorted)):
        result, pvalue = result_pvalue_sorted[i]
        threshold = (i + 1) * base_threshold
        if pvalue > threshold:
            discovered = result_pvalue_sorted[:i]
            if len(discovered) == 0:
                return [], [], threshold
            discovered_results, discovered_pvalues = zip(*discovered)[0], zip(*discovered)[1]
            return discovered_results, discovered_pvalues, threshold

    return results, pvalues, threshold


def main():

    # Test code
    Pvalues = [0.001, 0.008, 0.027, 0.039, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216, 0.222, 0.251, 0.269, 0.275, 0.34,
          0.341, 0.569, 0.594, 0.696, 0.762, 0.94, 0.942, 0.975, 0.986, 1.0]
    results = [str(p) for p in Pvalues]
    
    reject_r, p_r = BH(results, Pvalues, Q=0.25)
    print reject_r, p_r
    
    reject = BH(results, Pvalues, Q=0.25, use_dependencies=True)    
    print reject


if __name__ == '__main__':
    main()
