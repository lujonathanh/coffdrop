#COFFDROP
Coffdrop is a comprehensive analysis toolkit for mutually exclusive and co-occurring pairs  of mutations in sequenced tumor samples, with False Discovery Rate control. Coffdrop is written in Python. It was developed by Jonathan Lu, Jason Pitt, and Lorenzo Pesce at the University of Chicago.

Coffdrop implements a binomial statistical model to assess for significance of the mutual exclusivity/co-occurrence of a pair. To control false discoveries, one can limit the tested pairs by performing an initial screen of the pairs over the patients with the least mutations, then choosing only the most significant ones to test across the whole distribution. Coffdrop uses the Benjamini-Hochberg procedure to control False Discoveries.

After detecting significant pairs, Coffdrop 
    1. searches for enriched genes and chromosomal regions
    2. searches for enriched pairs.
    3. plots the mutual-exclusivity and co-occurrence networks and finds genes with the highest degree centrality
    4. Searches for triplets of mixed mutually exclusive and co-occurring pairs.

Furthermore, it has a flexible preprocessing feature to allow for:
    1. handling various mutation types, particularly Copy Number Alterations, which can create significant artefacts due to lack of independence among alterations in nearby genes. Thus, one can require genes to be a certain distance away before being run
    2. testing only those genes above a certain frequency

##Requirements
Coffdrop requires the following Python modules:
    1. NetworkX
    2. SciPy
    3. NumPy
    4. matplotlib

#Usage
See "Coffdrop workflow" in the wiki for details.

#Input
Coffdrop provides several python scripts for processing and integrating MAF and GISTIC files into the alteration matrix format, detailed below.

Alteration matrix. This tab-separated file lists alterations in your dataset. Each row lists the alterations for a single sample. In each row, the first column lists a sample ID, and the remaining columns list genes that are altered in that sample. Note that the matrix is not necessarily symmetric, as different samples will have different numbers of alterations.
In all files, lines starting with '#' are ignored.


#Output
Output files are txt files with each identified pair or triplet as one row.

We provide example matrices (".m2") in the data folder.