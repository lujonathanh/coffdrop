__author__ = 'jlu96'

import mutex as mex
import mutex_triangles as met
import time
import csv
import sys

global missing_chromosome
missing_chromosome = 'Z'

def getgenepairs(geneToCases, genes1, genes2=None, test_minFreq=0, closer_than_distance=None, only_filter_copy_distance=True,
                 remove_missing_genes=True):
    """
    :param genes1: First list of genes.
    :param genes2: Second list of genes. If None, defaults to making pairs from the first gene list.
    :return: genepairs, a set of all the possible genepairs between genes1 and genes2
    """
    t = time.time()

    if not genes2:
        genes2 = genes1

    relevant_genes = set([gene for gene in geneToCases.keys() if len(geneToCases[gene]) >= test_minFreq])

    if remove_missing_genes:
        print "All genes that weren't found in gene_positions file will not be considered"
        for gene in list(relevant_genes):
            if get_segment_gene_info(gene)['Chromosome'] == missing_chromosome:
                relevant_genes.remove(gene)

    genepairs = set()

    for gene1 in genes1:
        for gene2 in genes2:
            if gene1 != gene2:
                if gene1 in relevant_genes and gene2 in relevant_genes:
                    genepair = frozenset(sorted([gene1, gene2]))
                    try:
                        if only_filter_copy_distance:
                            if not closer_than_distance:
                                genepairs.add(genepair)
                            elif (not is_segment(gene1) or not is_segment(gene2)): # if one of them is somatic it doesn't matter
                                genepairs.add(genepair)
                            elif not check_pair_same_segment(genepair, bin_distance_threshold=closer_than_distance):
                                    # Remove those genes that weren't found in the
                                genepairs.add(genepair)
                        else:
                            if not closer_than_distance or not check_pair_same_segment(genepair, bin_distance_threshold=closer_than_distance):
                                genepairs.add(genepair)
                    except KeyError:
                        pass


    print "Number of pairs to test: ", len(genepairs)
    return list(genepairs)


def get_gene_bins_cooccur_same_segment(pairdict, geneToCases, cratiothresh, mutfreqdiffthresh,
                                   mutfreqdiffratiothresh, coveragethresh, probabilitythresh,
                                    bincol='RoundedLogPCov', bin_distance_threshold=2000000):

    print "Binning genes by pairs...."
    print "Bin Distance threshold is ", bin_distance_threshold

    t = time.time()
    geneToBin = {}
    for gene in geneToCases:
        geneToBin[gene] = set([gene])

    sorted_pairs = sorted(pairdict.keys(), key=lambda pair: pairdict[pair][bincol])


    # Check the gene positions right here!
    for pair in sorted_pairs:
        if pairdict[pair]['Concordance'] and pairdict[pair]['Probability'] < probabilitythresh:
            if pairdict[pair]['CooccurrenceRatio'] >= cratiothresh \
                and pairdict[pair]['MutationFrequencyDifference'] <= mutfreqdiffthresh and \
                pairdict[pair]['MutationFrequencyDifferenceRatio'] < mutfreqdiffratiothresh \
                and pairdict[pair]['Coverage'] > coveragethresh \
                and check_pair_same_segment(pair, bin_distance_threshold=bin_distance_threshold):

                gene0, gene1 = tuple(pair)
                geneToBin[gene0] = geneToBin[gene0].union(geneToBin[gene1])
                for gene in geneToBin[gene0]:
                    geneToBin[gene] = geneToBin[gene0]


    # Convert the bins into the new names.
    for gene in geneToBin:
        bin_set = frozenset(geneToBin[gene])
        suffix = iter(bin_set).next()[-4:]
        bin_name = '_'.join(sorted([the_bin[:-4] for the_bin in bin_set])) + suffix
        geneToBin[gene] = bin_name

    # There's probably some unnecessary stuff below, but stability over code efficiency. I think we can just
    # jump straight to bin sets.


    print len(geneToBin), " genes reduced to ", len(set(geneToBin.values())), " bins"
    print "Ratio threshold: ", cratiothresh, ', Diff ratio threshold: ', mutfreqdiffratiothresh, ', Coverage threshold: ', coveragethresh, ', Pvalue threshold: ', probabilitythresh

    return geneToBin




def bin_sets_from_geneToBin(genes, geneToBin):
    # Return: geneToBinSet, bin_setToBin
    # Given a mapping of the gene to the segment, get the set of all genes (bin_set) within that segment, mapped to that segment (Bin)

    bins = set(geneToBin.values())

    geneToBinSet = {}
    bin_setToBin = {}


    # Load the genes already in geneToBin
    for bin in bins:
        bin_set = set()
        for gene in geneToBin:
            if geneToBin[gene] == bin:
                bin_set.add(gene)

        for gene in bin_set:
            geneToBinSet[gene] = frozenset(bin_set)

        bin_setToBin[frozenset(bin_set)] = bin



    return geneToBinSet, bin_setToBin




def update_geneToCases_patientToGenes(geneToCases, patientToGenes, bin_setToBin, at_least_half=True):



    newGeneToCases = {}

    for bin_set in bin_setToBin:
        bin = bin_setToBin[bin_set]

        newGeneToCases[bin] = set()

        patients = set.union(*[geneToCases[gene] for gene in bin_set])

        for patient in patients:
            # find out how many genes are held by patients
            if (not at_least_half) or len(patientToGenes[patient].intersection(bin_set)) >= len(bin_set)/2:
                newGeneToCases[bin].add(patient)

    newPatientToGenes = {}
    for patient in patientToGenes:
        newPatientToGenes[patient] = set()

    for gene in newGeneToCases:
        for patient in newGeneToCases[gene]:
            newPatientToGenes[patient].add(gene)

    print "Old number of genes: ", len(geneToCases.keys())
    print "New number of genes: ", len(newGeneToCases.keys())

    return newGeneToCases, newPatientToGenes



def get_gene_bin_entries(geneToCases, newGeneToCases, geneToBinSet, bin_setToBin):
    gene_bin_entries = {}

    for gene in geneToBinSet:
        bin_name = bin_setToBin[geneToBinSet[gene]]

        # consolidate info
        entry = {}
        entry['Gene'] = gene
        entry['Bin'] = bin_name
        gene_cases = geneToCases[gene]
        bin_cases = newGeneToCases[bin_name]
        shared_cases = gene_cases.intersection(bin_cases)
        total_cases = gene_cases.union(bin_cases)
        entry['GeneMutationFrequency'] = len(gene_cases)
        entry['BinMutationFrequency'] = len(bin_cases)
        entry['SharedPatients'] = len(shared_cases)
        entry['TotalPatients'] = len(total_cases)
        entry['SharedRatio'] = entry['SharedPatients'] * 1.0 / entry['TotalPatients']
        gene_bin_entries[gene] = entry

    return gene_bin_entries


def convert_genes_to_bins(genes, geneToBin):

    return set([geneToBin[gene] for gene in genes])



def check_pair_same_segment(pair, bin_distance_threshold=2000000, ):
    gene0, gene1 = tuple(pair)
    gene0info = get_segment_gene_info(gene0)
    gene1info = get_segment_gene_info(gene1)

    # return gene0info['Chromosome'] == gene1info['Chromosome'] and \
    #        (min([abs(gene0info['Start'] - gene1info['End']), abs(gene0info['End'] - gene1info['Start'])]) \
    #            < bin_distance_threshold)

    if gene0info['Chromosome'] == gene1info['Chromosome']:
        indices = [gene0info['Start'], gene0info['End'], gene1info['Start'], gene1info['End']]
        sorted_indices = sorted(indices)
        if indices[0:2] == sorted_indices[0:2] or indices[0:2] == sorted_indices[2:]:
            if abs(sorted_indices[2] - sorted_indices[1]) > bin_distance_threshold:
                return False
            else:
                return True
        else:
            return True
    else:
        return False


def get_gene_distance(gene1, gene2, usestart=True):
    try:
        geneToPosition = get_gene_distance.geneToPosition
    except AttributeError:
        get_gene_distance.geneToPosition = load_gene_positions()
        geneToPosition = get_gene_distance.geneToPosition

    if gene1 not in geneToPosition:
            #print 'Position for gene ', gene1, " not found"
            return None
    elif gene2 not in geneToPosition:
        #print 'Position for gene ', gene2, " not found"
        return None
    else:
        gene1pos = geneToPosition[gene1]
        gene2pos = geneToPosition[gene2]

    # POS- change equals to check for tuples, go for each
    if gene1pos['Chromosome'] == gene2pos['Chromosome']:
        indices = [gene1pos['Start'], gene1pos['End'], gene2pos['Start'], gene2pos['End']]
        sorted_indices = sorted(indices)
        if indices[0:2] == sorted_indices[0:2] or indices[0:2] == sorted_indices[2:]:
            return abs(sorted_indices[2] - sorted_indices[1])
        else:
            return 1
    else:
        return None


def is_segment(segment):
    return segment[-4:] in {'loss', 'gain'}


def get_segment_gene_info(segment):
    try:
        get_segment_gene_info.geneToStat
        get_segment_gene_info.segToStat
    except AttributeError:
        get_segment_gene_info.geneToStat = load_gene_positions()
        get_segment_gene_info.segToStat = {}

    try:
        return get_segment_gene_info.segToStat[segment]
    except KeyError:

        seg_info = {}


        if segment[-4:] in {'loss', 'gain'}:
            alt_type = segment[-4:]
            genes = segment[:-4].split('_')
        else:
            alt_type = 'somatic'
            genes = segment.split('_')

        # alt_type = segment[-4:]


        seg_info['Name'] = segment

        in_dictionary = False
        for gene in genes:
            if gene in get_segment_gene_info.geneToStat and get_segment_gene_info.geneToStat[gene]['Chromosome'] != missing_chromosome:
                seg_info['Chromosome'] = get_segment_gene_info.geneToStat[gene]['Chromosome']
                in_dictionary = True
                break
            else:
                if gene not in get_segment_gene_info.geneToStat:
                    get_segment_gene_info.geneToStat[gene] = {}
                    get_segment_gene_info.geneToStat[gene]['Chromosome'] = missing_chromosome
                    get_segment_gene_info.geneToStat[gene]['Start'] = 0
                    get_segment_gene_info.geneToStat[gene]['End'] = 1

        if not in_dictionary:
            seg_info['Chromosome'] = missing_chromosome
            for gene in genes:
                get_segment_gene_info.geneToStat[gene] = {}
                get_segment_gene_info.geneToStat[gene]['Chromosome'] = missing_chromosome
                get_segment_gene_info.geneToStat[gene]['Start'] = 0
                get_segment_gene_info.geneToStat[gene]['End'] = 1


        seg_info['Genes'] = genes
        seg_info['Number'] = len(genes)
        seg_info['Type'] = alt_type

        seg_info['Start'] = min([get_segment_gene_info.geneToStat[gene]['Start'] for gene in genes])
        seg_info['End'] = max([get_segment_gene_info.geneToStat[gene]['End'] for gene in genes])

        seg_info['Length'] = seg_info['End'] - seg_info['Start']

        # LEFT OFF HERE 8/4/15 -jlu

        cytobands = get_cytobands2(seg_info['Chromosome'], seg_info['Start'], seg_info['End'])
        seg_info['Cytobands'] = cytobands
        seg_info['JoinedCytobands'] = '_'.join(cytobands)

        seg_info['ID'] = ':'.join([str(entry) for entry in [seg_info['Chromosome'], seg_info['Start'],  seg_info['End'],
                                                            seg_info['Length'], seg_info['Number'], seg_info['Type']]])

        get_segment_gene_info.segToStat[segment] = seg_info

        if not in_dictionary:
            print "Gene ", segment, " was not found in gene_positions.txt. We pretend it is on Chromosome Z."

        return seg_info


def write_segment_infos(segments, filename, header=['Name',  'Number', 'Type', 'Genes', 'Chromosome', 'Start', 'End', 'Length', 'Cytobands', 'ID']):

    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, delimiter='\t', extrasaction='ignore', fieldnames=header)
        writer.writeheader()

        for segment in segments:
            writer.writerow(get_segment_gene_info(segment))


def load_cytobands2(filename='cytoBand.txt'):
    chrom_dict = {}
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            chrom = row[0].strip('chr')
            start = eval(row[1])
            end = eval(row[2])
            cytoband = row[3]

            if chrom not in chrom_dict:
                chrom_dict[chrom] = []

            chrom_dict[chrom].append((start, end, str(chrom) + cytoband))
    for chrom in chrom_dict:
        chrom_dict[chrom] = sorted(chrom_dict[chrom], key=lambda entry: entry[0]) # sort ascending

    return chrom_dict

def get_cytobands2(chrom, start, end):

    try:
        chrom_dict = get_cytobands2.chrom_dict
    except AttributeError:
        get_cytobands2.chrom_dict = load_cytobands2()
        chrom_dict = get_cytobands2.chrom_dict



    if chrom not in chrom_dict:
        return []

    cyto_tuples = chrom_dict[chrom] # these are in the form start, end, cytoband


    cytobands = [entry[2] for entry in cyto_tuples if (entry[0] <= start and entry[1] >= start) or
                 (entry[0] >= start and entry[0] <= end and entry[1] >= start and entry[1] <= end) or
                 (entry[0] <= end and entry[1] >= end) or (entry[0] <= start and entry[1] >= end)] # first intersection
    #
    #
    # for cyto in cyto_dict:
    #     if chrom == cyto_dict[cyto]['Chromosome']:
    #         if not ((cyto_dict[cyto]['End'] < start and cyto_dict[cyto]['Start'] < start)
    #             or (cyto_dict[cyto]['End'] >= end and cyto_dict[cyto]['Start'] >= end )):
    #             cytobands.append(cyto)

    return cytobands




def load_gene_positions(filename='gene_positions.txt', genecol='Associated Gene Name', chromcol='Chromosome Name', extra_genecols=['Ensembl Gene ID'],
                   startcol='Gene Start (bp)', endcol='Gene End (bp)', delimiter='\t',
                   chromosomes = set(['X', 'Y'] + [str(i) for i in range(1, 23)]),
                        include_extras=False):
    t = time.time()
    geneToStat = {}
    extra_genes = set()
    with open(filename,'rU') as geneToStatFile:
        statreader = csv.DictReader(geneToStatFile, delimiter=delimiter)
        for row in statreader:
            if row[chromcol] in chromosomes:
                gene = row[genecol]
                position = {}
                position['Chromosome'] = row[chromcol]
                position['Start'] = eval(row[startcol])
                position['End'] = eval(row[endcol])

                if gene not in geneToStat:
                    geneToStat[gene] = position
                elif include_extras:
                    print "Extra entry for ", gene
                    if not isinstance(geneToStat[gene]['Chromosome'], tuple):
                        geneToStat[gene]['Chromosome'] = (geneToStat[gene]['Chromosome'],)
                        geneToStat[gene]['Start'] = (geneToStat[gene]['Start'],)
                        geneToStat[gene]['End'] = (geneToStat[gene]['End'],)
                    geneToStat[gene]['Chromosome'] += (position['Chromosome'],)
                    geneToStat[gene]['Start'] += (position['Start'],)
                    geneToStat[gene]['End'] += (position['End'],)
                    extra_genes.add(gene)

                if extra_genecols:
                    genes = [row[extra_genecol] for extra_genecol in extra_genecols]
                    for gene in genes:
                        if gene not in geneToStat:
                            geneToStat[gene] = position

    if not include_extras:
        return geneToStat
    else:
        return geneToStat, extra_genes


def write_gene_positions(genes, filename='gene_positions.txt', genecol='Associated Gene Name', chromcol='Chromosome Name',
                   startcol='Gene Start (bp)', endcol='Gene End (bp)', delimiter='\t',
                   chromosomes = set(['X', 'Y'] + [str(i) for i in range(1, 23)])):
    # Get header and other irrelevant columns
    with open(filename, 'rU') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=delimiter)
        header = reader.fieldnames


    other_cols = set(header)
    other_cols.remove(genecol)
    other_cols.remove(chromcol)
    other_cols.remove(startcol)
    other_cols.remove(endcol)

    with open(filename, 'ab') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter=delimiter, extrasaction='ignore')
        for gene in genes:
            seg_info = get_segment_gene_info(gene)
            seg_info[genecol] = seg_info['Name']
            seg_info[chromcol] = seg_info['Chromosome']
            seg_info[startcol] = seg_info['Start']
            seg_info[endcol] = seg_info['End']

            for col in other_cols:
                seg_info[col] = 'NA'

            writer.writerow(seg_info)



def load_gene_to_bin(seg_info, geneToCases, no_throw_out_extras=False, is_gene2seg=False):
    # Load the segments from another file and map the genes to them
    # is_segment_report means that the each
    relevant_genes = set(geneToCases.keys())

    geneToBin = {}

    # Left off here. Use bgbp to load and rewrite AllCancersmutation matrix. Then segment again on this new mutation matrix.
    if is_gene2seg:
        suffixes = ['gain', 'loss', '']
        with open(seg_info, 'rU') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                gene = row['Gene']
                segment = row['Segment']

                gene_segment_pairs = [(gene + suffix, segment + suffix) for suffix in suffixes]

                for gene_segment_pair in gene_segment_pairs:
                    if gene_segment_pair[0] in relevant_genes:
                        geneToBin[gene_segment_pair[0]] = gene_segment_pair[1]

    else:
        with open(seg_info, 'rU') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                segment = row['Name']
                suffix = segment[-4:]

                segment_genes = [gene + suffix for gene in segment[:-4].split('_')]

                for gene in segment_genes:
                    if gene in relevant_genes:
                        geneToBin[gene] = segment

    extra_genes = relevant_genes.difference(set(geneToBin.keys()))
    print len(extra_genes), " genes out of ", len(relevant_genes), " were not found in loading segment file"

    if no_throw_out_extras:
        print "Setting these extras to their own bin."
        for gene in relevant_genes:
            if gene not in geneToBin:
                    geneToBin[gene] = gene
                    extra_genes.remove(gene)
    else:
        print "Removing these extras from geneToCases"



    return geneToBin, extra_genes







# def update_geneToCases(geneToCases, bin_to_ID):
#
#     return newGeneToCases
#
#
# def update_patientToGenes(patientToGenes, bin_to_ID):
#
#     return newPatientToGenes



def writemutationmatrix(patientToGenes, filename):
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')

        for patient in patientToGenes:
            writer.writerow([patient] + list(patientToGenes[patient]))




def get_parser():
    # Parse arguments
    import argparse

    description = 'Given an input number of samples n, and probability of mutual exclusivity p, ' \
                  'plots the number of mutations that each sample' \
                  'must have in order to reach that probability p.'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-m', '--mutation_matrix', required=True,
                        help='File name for mutation data.')
    parser.add_argument('-nm', '--newmutationmatrix', required=True)


    parser.add_argument('-o', '--output_prefix', default=None,
                        help='Output path prefix (TSV format). Otherwise just prints.')
    parser.add_argument('-mf', '--min_freq', type=int, default=0,
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-tp', '--top_percentile', type=float, default=100,
                        help='Limit to this percentage of mutations of greatest abundance.')
    parser.add_argument('-tn', '--top_number', type=int, default=0,
                        help='Limit to this number of mutations of greatest abundance.')

    parser.add_argument('-pf', '--patient_file', default=None,
                        help='File of patients to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None,
                        help='File of genes to be included (optional).')
    parser.add_argument('-pbf', '--patient_blacklist_file', default=None,
                        help='File of patients to be excluded (optional).')
    parser.add_argument('-gbf', '--gene_blacklist_file', default=None,
                        help='File of genes to be excluded (optional).')


    parser.add_argument('-gl1', '--gene_list_1', default=None,
                        help='First sets of genes to draw from')

    parser.add_argument('-gl2', '--gene_list_2', default=None,
                        help='Second set of genes to draw from')




    #
    # parser.add_argument('-s', '--set_size', type=int, default=2, help='Number of genes per set')
    # parser.add_argument('-t', '--type', default='m',
    #                     help='Use m for mutual exclusivity, c for cooccuring')

    parser.add_argument('-mo', '--max_overlaps', type=int, default=400,
                        help='Maximum allowed number of overlapping mutations per mutually exclusive set')
    parser.add_argument('-mc', '--min_cooccur', type=int, default=1,
                        help='Minimum number of cooccurrences per cooccuring set')
    parser.add_argument('-mcr', '--min_cooccurrence_ratio', type=float, default=0.0,
                        help='The minimum cooccurrence ratio for a set to be deemed significant.')
    parser.add_argument('-cdt', '--cooccur_distance_threshold', type=float, default=None,
                    help='The maximum distance between two cooccuring genes to be deemed significant.')


    parser.add_argument('-bdt', '--bin_distance_threshold', type=float, default=10000000,
                help='The maximum distance between two cooccuring genes to be deemed significant.')



    parser.add_argument('-mp', '--mprob', type=float, default=0.05,
                        help='Significance Threshold for mutual exclusivity.')
    parser.add_argument('-cp', '--cprob', type=float, default=0.05,
                    help='Significance Threshold for mutual exclusivity.')


    parser.add_argument('-mdf', '--mdictfile', default=None, help='Mutual Exclusivity pair file')
    parser.add_argument('-cdf', '--cdictfile', default=None, help='Cooccurring pair file')


    parser.add_argument('-pcn', '--parallel_compute_number', default=None, type=int, help='Number of parallel computations'
                                                                                          'to use.')


    parser.add_argument('-fcoding', '--filter_coding', type=int, default=0, help='Filter only by those genes that do code.')

    parser.add_argument('-fmgl', '--filter_mutex_gain_loss', type=int, default=1, help='Filter the same-gene gainloss mutual exclusivity')
    parser.add_argument('-fmss', '--filter_mutex_same_segment', type=int, default=0, help='Pairs with gain and loss filtered.')



    parser.add_argument('-fcss', '--filter_cooccur_same_segment', type=int, default=0, help='Pairs filtered as below.')
    parser.add_argument('-fcss_rt', '--fcss_cratiothresh', type=float, default=0.0, help='Maximum ratio for cooccurring pairs,'
                                                                                         'which are probably the same at this point')
    parser.add_argument('-fcss_mfdrt', '--fcss_mutfreqdiffratiothresh', type=float, default=1, help='Maximum mutation frequency'
                                                                                                       'difference ratio for cooccurring pairs,'
                                                                                     'which are probably the same at this point')
    parser.add_argument('-fcss_mfdt', '--fcss_mutfreqdiffthresh', type=float, default=100, help='Maximum mutation frequency difference.')
    parser.add_argument('-fcss_ct', '--fcss_coveragethresh', type=float, default=0, help='Minimum coverage for cooccurring pairs,'
                                                                                     'which are probably the same at this point')
    parser.add_argument('-fcss_pt', '--fcss_probabilitythresh', type=float, default=1e-20, help='Maximum ratio for cooccurring pairs,'
                                                                             'which are probably the same at this point')

    parser.add_argument('-bgbp', '--bin_genes_by_pairs', default=False, help='Bin the genes by the pairs in a second iteration.')
    parser.add_argument('-rpt', '--ridiculous_pvalue_thresh', type=float, default=1e-10, help="Filter ridiculous pvalues")


    parser.add_argument('-glc', '--gainlosscombine', type=int, default=0, help='Combine gain and loss in file')
    parser.add_argument('-fc', '--filter_cooccur', default=False, help='Filter cooccurring pairs')
    parser.add_argument('-ud', '--use_downstream', default=True, help='Use the cooccurring pairs downstream')

    parser.add_argument('-r_p', '--ratioperc', type=float, default=50, help='Cooccurrence Ratio Percentile')
    parser.add_argument('-s_p', '--scoreperc', type=float, default=70, help='Set Score Percentile')
    parser.add_argument('-c_p', '--coverageperc', type=float, default=30, help='Coverage Percentile')
    parser.add_argument('-cs_p', '--combscoreperc', type=float, default=0.0, help='Combined Score Percentile')

    parser.add_argument('-leb', '--local_edge_bet', default=None, help='Calculate the local edge betweenness.')

    parser.add_argument('-gt', '--group_type', default='TripletNetwork', help='Type of group to find')
    parser.add_argument('-pt', '--pair_type', default=None, help='Type of pair to limit to')


    parser.add_argument('-ctst', '--calc_triple_stats', default=True, help='Calculate the cooccurrence ratio and overlap'
                                                                          'and coverage of Triplet.')
    parser.add_argument('-ctsc', '--calc_triple_scores', default=False, help='Calculate the set score, etc. of Triplet.')
    parser.add_argument('-om', '--only_mixed', default=True, help="Only search for mixed triplets.")

    # parser.add_argument('-v', '--cytoscape_output', type=int, default=1, help='Set to 1 to produce cytoscape-compatible'
    #                                                                           'outputs')
    # parser.add_argument('-co', '--complete_output', type=int, default=1, help='Set to 1 to produce a file with complete information')

    #parser.add_argument('-cif', '--check_in_file', default=None, help='Name of file to check if the set is in')


    parser.add_argument('-js', '--just_sets', type=int, default=0, help='Just find the possible cooccurring sets,'
                                                                        'without checking for significance')
    parser.add_argument('-pd', '--plot_distribution', type=int, default=0,
                        help='Plot distribution of mutation frequencies '
                             'with this number of bins before running')
    parser.add_argument('-ppd', '--plot_patient_distribution', type=int, default=0,
                        help='Plot distribution of mutation number'
                             'per gene '
                             'with this number of bins before running')

    parser.add_argument('-lgs', '--load_gene_segments', default=None)
    parser.add_argument('-gsf', '--gene_segment_file', default=None)
    parser.add_argument('-ig2s', '--is_gene2seg', type=int, default=0, help="If the file is a gene2seg fie output from cna2matrix.py, as opposed to SEGMENTINFO file.")
    parser.add_argument('-gbef', '--gene_bin_entries_file', default=None)
    parser.add_argument('-sgif', '--segment_info_file', default=None, help='Name of the file to write SEGMENTINFO to.')

    parser.add_argument('-ntoe', '--no_throw_out_extras', type=int, default=-0)



    return parser

def run(args):
    mdictfile = args.mdictfile
    cdictfile = args.cdictfile

    mprob = args.mprob
    cprob = args.cprob
    cooccur_distance_threshold = args.cooccur_distance_threshold
    bin_distance_threshold = args.bin_distance_threshold

    mutationmatrix = args.mutation_matrix
    newmutationmatrix = args.newmutationmatrix

    file_prefix = args.output_prefix
    if not file_prefix:
        file_prefix = newmutationmatrix

    geneFile = args.gene_file
    patientFile = args.patient_file
    gene_blacklist = args.gene_blacklist_file
    patient_blacklist = args.patient_blacklist_file
    minFreq = args.min_freq
    minCooccur = args.min_cooccur
    min_cooccurrence_ratio = args.min_cooccurrence_ratio
    top_percentile = args.top_percentile
    top_number = args.top_number

    parallel_compute_number = args.parallel_compute_number


    filter_cooccur_same_segment = args.filter_cooccur_same_segment
    fcss_cratiothresh = args.fcss_cratiothresh
    fcss_mutfreqdiffthresh = args.fcss_mutfreqdiffthresh
    fcss_mutfreqdiffratiothresh = args.fcss_mutfreqdiffratiothresh
    fcss_coveragethresh = args.fcss_coveragethresh
    fcss_probabilitythresh = args.fcss_probabilitythresh



    gene_segment_file = args.gene_segment_file
    load_gene_segments = args.load_gene_segments
    is_gene2seg = args.is_gene2seg
    gene_bin_entries_file = args.gene_bin_entries_file
    no_throw_out_extras = args.no_throw_out_extras
    segment_info_file = args.segment_info_file


    if not gene_bin_entries_file:
        gene_bin_entries_file = file_prefix + '_binnedgenes.tsv'

    if not segment_info_file:
        segment_info_file = file_prefix + '_SEGMENTINFO.tsv'

    #-----------------------------------------------------






    mutations = mex.remove_blacklists(gene_blacklist, patient_blacklist,
                                  *mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq))
    numGenes, numCases, genes, patients, geneToCases, patientToGenes = mutations
    print 'Filtered Mutation data: %s genes x %s patients' % (numGenes, numCases)


    # Load segment info
    if load_gene_segments:

        # extra_genes is the genes not found in the segment file.
        # If throw_out_extras is False, extra_genes will be empty.
        geneToBin, extra_genes = load_gene_to_bin(gene_segment_file, geneToCases, no_throw_out_extras=no_throw_out_extras, is_gene2seg=is_gene2seg)

        numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.remove_extra_genes(extra_genes, numGenes, numCases, genes, patients, geneToCases, patientToGenes)


    else:
        print "Beginning bin genes by co-occurring pairs. "
        genepairs = getgenepairs(geneToCases, genes, closer_than_distance=bin_distance_threshold)
        print "Pairs retrieved. Calculating cooccurring pairs to make bins."

        cpairsdict, cgenedict = met.complete_cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, fcss_probabilitythresh, minCooccur,
                      cooccur_distance_threshold, fcss_cratiothresh, parallel_compute_number,
                      filter_cooccur_same_segment, fcss_cratiothresh, fcss_mutfreqdiffratiothresh,
                      fcss_coveragethresh, fcss_probabilitythresh)

        print "Cooccurring pairs calculated."
        geneToBin = get_gene_bins_cooccur_same_segment(cpairsdict, geneToCases, fcss_cratiothresh, fcss_mutfreqdiffthresh,
                           fcss_mutfreqdiffratiothresh, fcss_coveragethresh, fcss_probabilitythresh, bin_distance_threshold=bin_distance_threshold)
        # Write these new bins out
        new_bins = convert_genes_to_bins(genes, geneToBin)
        write_segment_infos(new_bins, filename=segment_info_file)
        print "New SEGMENTINFO written to ", segment_info_file

        write_gene_positions(new_bins)
        print "New segment positions appended to gene_positions.txt"


    # Update to the new mutation matrix.
    geneToBinSet, bin_setToBin = bin_sets_from_geneToBin(genes, geneToBin)

    newGeneToCases, newPatientToGenes = update_geneToCases_patientToGenes(geneToCases, patientToGenes, bin_setToBin, at_least_half=True)

    gene_bin_entries = get_gene_bin_entries(geneToCases, newGeneToCases, geneToBinSet, bin_setToBin)

    if gene_bin_entries_file:
        met.writeanydict(gene_bin_entries, gene_bin_entries_file)
        print "Gene bin entries written to ", gene_bin_entries_file


    # Write the new mutation matrix out.
    writemutationmatrix(newPatientToGenes, filename=newmutationmatrix)



def main():
    run(get_parser().parse_args(sys.argv[1:]))
    # Write the new genes/patients out in a matrix

if __name__ == '__main__':
    main()
