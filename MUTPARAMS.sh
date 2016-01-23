# Centralized parameter file

# ************************************************************************************
# INPUT OUTPUT
# ************************************************************************************

export MUTATIONMATRIXPATH=/Users/jlu96/maf/new/PRAD_broad/PRAD_broad-som.m2

# This will be prefix all of the output files.
export OUTPUTPREFIX=/Users/jlu96/conte/jlu/Analyses/CooccurImprovement/Analyses/LowMutatednessScore/data/PRAD_broad-som


# ************************************************************************************
# LOADING PARAMETERS
# ************************************************************************************

# Limit to genes mutated in at least this # of patients
export MINFREQ=

# Limit to genes mutated in at least this % of patients
export MINPERCENTILE=

# Limit to this # of top most frequently mutated genes
export TOP_NUMBER=

# List of patients to limit loading to
export PATIENTFILE=

# List of Genes to limit loading to
export GENEFILE=

# Blacklists
export PATIENTBLACKLIST=
export GENEBLACKLIST=



# ************************************************************************************
# LIMIT WHICH PAIRS TO TEST
# ************************************************************************************

# If GENELIST1 and GENELIST2 are both set, only pairs containing one gene in GENELIST1 and one in
# GENELIST2 will be tested. If only GENELIST1 is set, only pairs containing two genes in GENELIST1
# Will be tested.
export GENELIST1=
export GENELIST2=

# If set, only test pairs from PAIRLIST. The column header must be "Gene0   Gene1"
export PAIRLIST=

# If set, only test pairs whose genes have a minimum mutation frequency at least this.
export TESTMINFREQ=4

# If set, only test genes that are at least this distance away.
# This will only apply to genes ending in "gain" or "loss".
# E.G. Say TP53 and MYC are 30 bp away. If the threshold is set to 1000000, the pair (TP53gain, MYCloss) will not
# be tested. But the pair (TP53, MYCloss) would be tested. This allows only Copy Number pairs to be
# filtered out
export MINDISTANCETHRESHOLD=



# ************************************************************************************
# CALCULATION PARAMETERS
# ************************************************************************************

# Set to 1 to generate Permutation matrices to test hypotheses.
export USE_PERM_MATRICES=
# Number of permutations. Default: 100
export NUM_PERMUTATIONS=
# Set to 1 to use the binary permutation method, not network
export BINARY_PERM_METHOD=
# Network Permutation Method will perform Q * |A| swaps, where A is the number of mutations in the
# matrix. Default: 1.
export Q=
# An existing directory to hold the written matrices
export PERM_MATRIX_DIRECTORY=
# Set to 1 to write the matrices to the directory
export WRITE_MATRICES=



# Mutex Pair PValue Threshold. Default is 0.05
export MPROB=
# Maximum Overlap between mutually exclusive pairs.
export MAXOVERLAP=



# Cooccur Pair PValue Threshold. Default is 0.05
export CPROB=
# Minimum overlap between cooccurring pairs
export MINCOOCCUR=
# minimum co-occurrence ratio for cooccuring pairs. This is the overlap/
export MINRATIO=


# If "Triplet", search for triplets among pairs. If "Network", calculate Network attributes
export GROUPTYPE=

# Number of processes to divide calculations across.
export NUMPROCESSES=
