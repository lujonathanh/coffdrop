# Centralized parameter file.
# binpairs.sh will convert a given mutation matrix to a segmented mutation matrix
#

# ************************************************************************************
# INPUT OUTPUT
# ************************************************************************************

# The input matrix
export OLDMATRIX=

# The output, segmented matrix.
export NEWSEGMATRIX=

# The name of the file to write the concordance of each gene with its reprsentative
# segment. Default is $OLDMATRIX_binnedgenes.tsv
export GENEBINENTRIESFILE=

# The bin option. To load, set to 1. To bin by co-occurrence, set to 0. Default: 0
export LOADGENESEGMENTS=



# ************************************************************************************
# BIN OPTION 1: LOADING
# ************************************************************************************

# Name of gene segments file
export GENESEGMENTFILE=SEGMENTINFO_COSMIC

# Type of gene segments file. Use of loading an .gene2seg file output from cna2matrix.py,
# instead of a SEGMENTINFO file output from binpairs.sh
export ISGENE2SEG=



# ************************************************************************************
# BIN OPTION 2: BIN BY CO-OCCURRENCE
# ************************************************************************************
# This option will find all co-occurring pairs

# The name of output file to write all SEGMENT information to.
export SEGMENTINFO=


# Maximum distance between adjacent genes to be in the same bin. Default: 1000000 bp
export DISTANCETHRESHOLD=
# Maximum co-occurrence p-value between adjacent genes to be in the same bin. Default: 1e-20
export PVALUETHRESHOLD=
# Minimum co-occurrence ratio between adjacent genes to be in the same bin (0.0-1.0). Default: 0.0
export RATIOTHRESHOLD=


# Number of processes to start and divide co-occurring pair testing across.
export NUMPROCESSES=




# ************************************************************************************
# MISCELLANEOUS LOADING PARAMETERS
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
