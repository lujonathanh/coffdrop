#Mutex
Mutex is a software package identifying mutually exclusive and co-occurring pairs and triplets of mutations in sequenced tumor samples. Mutex is written in Python. It was developd by Jonathan Lu and Jason Pitt in the Kevin White Lab at the Institute for Genomics and Systems Biology at the University of Chicago.

Mutex implements a permutation statistical model to assess for significance of the interaction. Furthermore, it has a flexible preprocessing feature to allow for handling various mutation types, particularly Copy Number Alterations, which can create significant artefacts due to lack of independence among alterations in nearby genes.

##Requirements
mutex requires the following Python modules:
1.  NetworkX
2.  SciPy
3.  NumPy

#Usage
See "Mutex workflow" in the wiki for details.

#Input
Alteration matrix. This tab-separated file lists alterations in your dataset. Each row lists the alterations for a single sample. In each row, the first column lists a sample ID, and the remaining columns list genes that are altered in that sample. Note that the matrix is not necessarily symmetric, as different samples will have different numbers of alterations.
In all files, lines starting with '#' are ignored.

Mutex provides several python scripts for processing and integrating MAF and GISTIC files.

#Output
Output files are tsv files with each identified pair or triplet as one row.

We provide example matrices (".m2") in the data folder.