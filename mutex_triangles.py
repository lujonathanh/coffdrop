__author__ = 'jlu96'

import mutex as mex
import csv
import collections
import time
from scipy import stats
import edgereader as edg
import networkx as nx
import sys
import parallel_compute_working as pac
import bingenesbypairs as bgbp


class Triplet:
    def __init__(self, number, pairdict=None):
        assert len(pairdict) == 3
        self.number = number
        self.pairdict = pairdict.copy()
        self.genes = set.union(*[set(pair) for pair in pairdict.keys()])

        types = [self.pairdict[pair]['Type'] for pair in self.pairdict]
        self.type = ''.join(sorted(types))

        self.stats = collections.OrderedDict()
        self.stats['ID'] = self.number
        self.calc_double_stats()

    def calc_double_stats(self):
        # Left off here, move writable dict to here? Quickly.
        gene0, gene1, gene2 = tuple(self.genes)
        self.stats['Gene0'], self.stats['Gene1'], self.stats['Gene2'] = gene0, gene1, gene2
        pair01, pair02, pair12 = frozenset([gene0, gene1]), frozenset([gene0, gene2]), frozenset([gene1, gene2]),
        self.stats['TripletType'] = self.type

        score_names = ['Coverage', 'SetScore', 'CooccurrenceRatio', 'CombinedScore', 'TotalEdgeBetweennessCentrality',
                       'AverageNodeClosenessCentrality', 'AverageNodeClusteringCoefficient', 'RoundedLogPCov']
        for score_name in score_names:
            pair_scores = [self.pairdict[pair][score_name] for pair in self.pairdict if score_name in self.pairdict[pair]]
            if pair_scores:
                self.stats['Average' + score_name] = sum(pair_scores)/len(pair_scores)
            else:
                self.stats['Average' + score_name] = 0

        self.stats['TripletScore'] = self.stats['AverageCombinedScore'] * 1.0 / self.stats['AverageRoundedLogPCov']

        self.stats['01Type'] = self.pairdict[pair01]['Type']
        self.stats['02Type']= self.pairdict[pair02]['Type']
        self.stats['12Type'] = self.pairdict[pair12]['Type']

        self.stats['01Concordance'] = self.pairdict[pair01]['Concordance']
        self.stats['02Concordance']= self.pairdict[pair02]['Concordance']
        self.stats['12Concordance'] = self.pairdict[pair12]['Concordance']

        self.stats['Somatic'] = 0
        for gene in tuple(self.genes):
            if gene[-4:] not in {'loss', 'gain'}:
                self.stats['Somatic'] += 1


        self.stats['Probability'] = mex.prod([self.pairdict[pair]['Probability'] for pair in self.pairdict])


    def calc_triple_stats(self, numCases, geneToCases, patientToGenes, compute_prob=False, calc_triple_scores=False):


        # self.overlap = mex.numoverlaps(self.genes, geneToCases)
        # self.coverage = mex.numcoverage(self.genes, geneToCases)
        # self.
        # # print "Genes are ", self.genes
        # # print "Overlap is", self.overlap
        # # print "Coverage is", self.coverage
        # # print "GeneToCases is", geneToCases

        if calc_triple_scores and self.type == 'CooccurringCooccurringCooccurring':
            triple_stats = mex.analyze_cooccur_set_new(numCases, geneToCases, patientToGenes, self.genes, compute_prob=compute_prob)

            triple_stats.update(self.stats)

            self.stats = triple_stats.copy()

            self.overlap =  self.stats['Overlap']
            self.coverage = self.stats['Coverage']
            self.ratio = self.stats['CooccurrenceRatio']
        else:
            self.overlap = mex.numoverlaps(self.genes, geneToCases)
            self.coverage = mex.numcoverage(self.genes, geneToCases)
            self.ratio = self.overlap * 1.0 / self.coverage
            self.stats['Overlap'] = self.overlap
            self.stats['Coverage'] = self.coverage
            self.stats['CooccurrenceRatio'] = self.ratio


    def getwritabledict(self):
        #[pair1, pair2, pair3] = self.pairdict.keys()

        return self.stats











def getTriplets(pairsdict, genesdict, pairstotest, numCases=None, geneToCases=None, patientToGenes = None, compute_prob=False, name='Triplets',
                calc_triple_stats=True, calc_triple_scores=False, only_mixed=False):
    """
    :param pairsdict:
    :param genesdict:
    :param numCases:
    :param geneToCases:
    :param patientToGenes:
    :param compute_prob:
    :param name:
    :param calc_triple_stats: Calc cooccurrence ratio, etc.
    :param calc_triple_scores: Calc the set scores, etc.
    :return:
    """


    t0 = time.time()
    Triplets = []
    num_Triplets = 0

    # Number of triplets each pair has
    print "number of pairs ", len(pairsdict)
    pairsdict_Triplets = pairsdict.copy()


    for pair in pairsdict:
        pairsdict_Triplets[pair][name] = set()
        pairsdict_Triplets[pair]['num' + name] = 0
        pairsdict_Triplets[pair]['type' + name] = set()

    print "Pair info of triplets initialized"


    # Number of triplets each gene is in
    genesdict_Triplets = {}
    for gene in genesdict:
        genesdict_Triplets[gene] = {}
        genesdict_Triplets[gene]['Gene'] = gene
        genesdict_Triplets[gene][name] = set()
        genesdict_Triplets[gene]['num' + name] = 0
        #genesdict_Triplets[gene]['type' + name] = set()


    gene_triples = set()

    # t = time.time()
    #
    # time
    for pair in pairstotest:
        # print "Getting pairs: ", pair
        gene1, gene2 = tuple(pair)
        if gene1 in genesdict and gene2 in genesdict:
            genes1 = genesdict[gene1]
            genes2 = genesdict[gene2]

            # The set of genes that cooccur with both genes.
            both = genes1.intersection(genes2)
            if both:
                for bothgene in both:
                    gene_triple = frozenset([gene1, gene2, bothgene])

                    bothpair1 = frozenset([gene1, bothgene])
                    bothpair2 = frozenset([gene2, bothgene])
                    # print "Both pairs: ", bothpair1, bothpair2

                    if (not only_mixed) or 'MutuallyExclusive' in [pairsdict[pair]['Type'], pairsdict[bothpair1]['Type'], pairsdict[bothpair2]['Type']]:

                        if gene_triple not in gene_triples:
                            gene_triples.add(gene_triple)




                            pairdict = {}
                            pairdict[pair] = pairsdict[pair]
                            pairdict[bothpair1] = pairsdict[bothpair1]
                            pairdict[bothpair2] = pairsdict[bothpair2]

                            newTriplet = Triplet(num_Triplets, pairdict=pairdict)
                            Triplets.append(newTriplet)

                            # pairsdict_Triplets[pair][name].add(num_Triplets)
                            # pairsdict_Triplets[pair]['num' + name] += 1
                            # pairsdict_Triplets[bothpair1][name].add(num_Triplets)
                            # pairsdict_Triplets[bothpair1]['num' + name] += 1
                            # pairsdict_Triplets[bothpair2][name].add(num_Triplets)
                            # pairsdict_Triplets[bothpair2]['num' + name] += 1


                            if calc_triple_stats:
                                newTriplet.calc_triple_stats(numCases, geneToCases, patientToGenes, compute_prob=compute_prob,
                                                 calc_triple_scores=calc_triple_scores)


                        for tpair in pairdict:
                            pairsdict_Triplets[tpair][name].add(num_Triplets)
                            pairsdict_Triplets[tpair]['num' + name] += 1
                            pairsdict_Triplets[tpair]['type' + name].add(newTriplet.type)

                        num_Triplets += 1

                        genesdict_Triplets[gene1][name].add(num_Triplets)
                        genesdict_Triplets[gene1]['num' + name] += 1
                        genesdict_Triplets[gene2][name].add(num_Triplets)
                        genesdict_Triplets[gene2]['num' + name] += 1
                        genesdict_Triplets[bothgene][name].add(num_Triplets)
                        genesdict_Triplets[bothgene]['num' + name] += 1

                #print len(both), "triplets calculated in ", time.time() - t0

    # for pair in pairsdict:
    #     if pair not in pairstotest:
    #         pairsdict_Triplets.pop(pair)

    print len(Triplets), " triplets calculated in ", time.time() - t0
    sorted_pairs = sorted(pairsdict_Triplets.keys(), key=lambda entry: pairsdict_Triplets[entry]['num' + name], reverse=True)
    sorted_genes = sorted(genesdict_Triplets.keys(), key=lambda entry: genesdict_Triplets[entry]['num' + name], reverse=True)

    print "Including sorting time ", time.time() - t0

    return Triplets, pairsdict_Triplets, sorted_pairs, genesdict_Triplets, sorted_genes


def combine_pairsdict_Triplets(pairsdict0, pairsdict1):
    # pairsdict_Triplets = pairsdict.copy()
    # for pair in pairsdict:
    #     pairsdict_Triplets[pair][name] = set()
    #     pairsdict_Triplets[pair]['num' + name] = 0
    #     pairsdict_Triplets[pair]['type' + name] = set()

    pairsdict = pairsdict0.copy()

    for pair in pairsdict1:
        if pair not in pairsdict:
            pairsdict[pair] = pairsdict1[pair]
        else:
            pairsdict[pair]['Triplets'] = pairsdict[pair]['Triplets'].union(pairsdict1[pair]['Triplets'])
            pairsdict[pair]['numTriplets'] += pairsdict1[pair]['numTriplets']
            pairsdict[pair]['typeTriplets'] = pairsdict[pair]['typeTriplets'].union(pairsdict1[pair]['typeTriplets'])

    return pairsdict


def combine_genesdict_Triplets(genesdict0, genesdict1):
    # genesdict_Triplets = {}
    # for gene in genesdict:
    #     genesdict_Triplets[gene] = {}
    #     genesdict_Triplets[gene]['Gene'] = gene
    #     genesdict_Triplets[gene][name] = set()
    #     genesdict_Triplets[gene]['num' + name] = 0

    genesdict = genesdict0.copy()

    for gene in genesdict1:
        if gene not in genesdict:
            genesdict[gene] = genesdict1[gene]
        else:
            genesdict[gene]['Triplets'] = genesdict[gene]['Triplets'].union(genesdict1[gene]['Triplets'])
            genesdict[gene]['numTriplets'] += genesdict1[gene]['numTriplets']

    return genesdict




def writeTriplets(Triplets, file_prefix, delimiter='\t', fieldnames=None):
    # Dictwriter
    # Get writabledict from each Triplet, make header from
    Tripletfile = file_prefix + '_Triplets.tsv'
    genefile = file_prefix + '_Tripletgenes.txt'

    with open(Tripletfile, 'w') as csvfile:
        if not fieldnames:
            fieldnames = Triplets[0].getwritabledict().keys()
        writer = csv.DictWriter(csvfile, delimiter=delimiter, fieldnames=fieldnames, extrasaction='ignore')

        writer.writeheader()
        for Triplet in Triplets:
            writer.writerow(Triplet.getwritabledict())

    print "Triplets written to ", file_prefix + '_Triplets.tsv'






def writeanydict(anydict, filename, delimiter='\t', fieldnames=None, orderedkeys=None):

    with open(filename, 'w') as csvfile:
        if not fieldnames:
            fieldnames = anydict.values()[0].keys()

        writer = csv.DictWriter(csvfile, delimiter=delimiter, fieldnames=fieldnames, extrasaction='ignore')

        writer.writeheader()

        if orderedkeys:
            for key in orderedkeys:
                writer.writerow(anydict[key])
        else:
            for entry in anydict.values():
                writer.writerow(entry)

def writegenedict(genedict, filename, delimiter='\t', fieldnames=None):
     with open(filename, 'w') as csvfile:
        if not fieldnames:
            fieldnames = ["Gene", "Entries"]

        writer = csv.writer(csvfile, delimiter=delimiter)

        writer.writerow(fieldnames)

        for gene in genedict:
            writer.writerow([gene] + [genedict[gene]])





def getgenepairs(geneToCases, genes1, genes2=None):
    """
    :param genes1: First list of genes.
    :param genes2: Second list of genes. If None, defaults to making pairs from the first gene list.
    :return: genepairs, a set of all the possible genepairs between genes1 and genes2
    """

    if not genes2:
        genes2 = genes1

    relevant_genes = set(geneToCases.keys())
    genepairs = set()

    for gene1 in genes1:
        for gene2 in genes2:
            if gene1 != gene2:
                if gene1 in relevant_genes and gene2 in relevant_genes:
                    genepair = frozenset([gene1, gene2])
                    genepairs.add(genepair)
    return list(genepairs)

def loadgenepairs(pair_list_file, geneToCases):
    relevant_genes = set(geneToCases.keys())
    genepairs = set()
    with open(pair_list_file, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if row[0] in relevant_genes and row[1] in relevant_genes:
                genepairs.add(frozenset(row))

    return list(genepairs)


def mutexpairs(numCases, geneToCases, patientToGenes, genepairs, p=1.0, maxOverlap=200, compute_prob=True):
    """
    :param numCases: The number of patients, total.
    :param geneToCases: Mapping of genes to cases.
    :param patientToGenes: Mapping of patients to genes
    :param p: Threshold for probability
    :param genepairs: The list of genepairs.
    :return: mpairsdict, mgenedict
        mpairsdict: A dictionary of sets which are significant, according to p, mapped to the statistics of that set.
        mgenedict: A dictionary of genes, mapped to the genes with which they are mutually exclusive.
    """

    mpairsdict = {}
    mgenedict = {}

    for genepair in genepairs:

        mstats = mex.analyze_mutex_set_new(numCases, geneToCases, patientToGenes, genepair, compute_prob=compute_prob)
        mprob = mstats['Probability']
        moverlap = mstats['Overlap']
        if mprob < p and moverlap <= maxOverlap:
            mpairsdict[genepair] = mstats
            gene1, gene2 = tuple(genepair)
            if gene1 not in mgenedict:
                mgenedict[gene1] = set()
                mgenedict[gene1].add(gene2)
            else:
                mgenedict[gene1].add(gene2)

            if gene2 not in mgenedict:
                mgenedict[gene2] = set()
                mgenedict[gene2].add(gene1)
            else:
                mgenedict[gene2].add(gene1)

    return mpairsdict, mgenedict







def cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, p=1.0, minCooccur=1, compute_prob=True, cooccur_distance_threshold=None,
                 min_cooccurrence_ratio=0.0):
    """
    :param numCases: The number of patients, total.
    :param geneToCases: Mapping of genes to cases.
    :param patientToGenes: Mapping of patients to genes
    :param p: Threshold for probability
    :param genepairs: The list of genepairs.
    :return: cpairsdict, cgenedict
        cpairsdict: A dictionary of sets which are significant, according to p, mapped to the statistics of that set.
        cgenedict: A dictionary of genes, mapped to the genes with which they are mutually exclusive.
    """

    cpairsdict = {}
    cgenedict = {}

    for genepair in genepairs:

        # Preliminary Distance Filter
        cratio = mex.cooccurrence_ratio(tuple(genepair), geneToCases)
        if cratio >= min_cooccurrence_ratio:

            cstats = mex.analyze_cooccur_set_new(numCases, geneToCases, patientToGenes, genepair, compute_prob=compute_prob,
                                                 getdistance=cooccur_distance_threshold)
            cprob = cstats['Probability']
            coverlap = cstats['Overlap']
            cratio = cstats['CooccurrenceRatio']
            if cprob < p and coverlap >= minCooccur and cratio >= min_cooccurrence_ratio:
                if not (cooccur_distance_threshold) or cstats['Distance'] > cooccur_distance_threshold:
                    cpairsdict[genepair] = cstats
                    gene1, gene2 = tuple(genepair)
                    if gene1 not in cgenedict:
                        cgenedict[gene1] = set()
                        cgenedict[gene1].add(gene2)
                    else:
                        cgenedict[gene1].add(gene2)

                    if gene2 not in cgenedict:
                        cgenedict[gene2] = set()
                        cgenedict[gene2].add(gene1)
                    else:
                        cgenedict[gene2].add(gene1)

    return cpairsdict, cgenedict


def filterpairs_new(pairsdict, genesdict, entryToFilter):
    """
    :param pairsdict:
    :param genesdict:
    :param entry: Entry of pairsdict.
    :param pair_filter: Function that returns True if want to keep.
    :return:
    """
    new_pairsdict = pairsdict.copy()
    new_genesdict = genesdict.copy()
    for pair in pairsdict:
        for entry in entryToFilter:
            if (entry not in pairsdict[pair]) or not entryToFilter[entry](pairsdict[pair][entry]):
                if (entry not in pairsdict[pair]):
                    print "Attribute not in pair"
                else:
                    pass
                    # print "Entry: ", entry
                    # print pairsdict[pair][entry]

                gene0, gene1 = tuple(pair)

                new_pairsdict.pop(pair)

                if gene0 in new_genesdict:
                    new_genesdict[gene0].remove(gene1)
                    if not new_genesdict[gene0]:
                        new_genesdict.pop(gene0)
                if gene1 in new_genesdict:
                    new_genesdict[gene1].remove(gene0)
                    if not new_genesdict[gene1]:
                        new_genesdict.pop(gene1)
                break

    return new_pairsdict, new_genesdict




def pairs_limittogenes(pairsdict, genesdict, genes):
    newpairsdict = {}
    newgenesdict = {}
    #geneset = set(genes)

    for gene in genes:
        if gene in genesdict:
            othergenes = genesdict[gene]
            for othergene in othergenes:
                pair = frozenset([gene, othergene])
                newpairsdict[pair] = pairsdict[pair]
                if gene not in newgenesdict:
                    newgenesdict[gene] = set()
                    newgenesdict[gene].add(othergene)
                else:
                    newgenesdict[gene].add(othergene)
                if othergene not in newgenesdict:
                    newgenesdict[othergene] = set()
                    newgenesdict[othergene].add(gene)
                else:
                    newgenesdict[othergene].add(gene)
            # bothgenes = geneset.intersection(othergenes)
            # if bothgenes:
            #     for bothgene in bothgenes:
            #         pair = frozenset([gene, bothgene])
            #         newpairsdict[pair] = pairsdict[pair]
            #         if gene not in newgenesdict:
            #             newgenesdict[gene] = {bothgene}
            #         else:
            #             newgenesdict[gene].add(bothgene)
            #         if bothgene not in newgenesdict:
            #             newgenesdict[bothgene] = {gene}
            #         else:
            #             newgenesdict[bothgene].add(gene)
    return newpairsdict, newgenesdict






def complete_mutexpairs(numCases, geneToCases, patientToGenes, genepairs, mprob, maxOverlap, parallel_compute_number,
                        filter_mutex_gain_loss, filter_mutex_same_segment):
    """
    Complete Wrapper function for mutexpairs.
    :return: mpairsdict, mgenedict
    """

    t = time.time()

    if parallel_compute_number:
        mpairsdict, mgenedict = pac.parallel_compute_new(mutexpairs, [numCases, geneToCases, patientToGenes, genepairs, mprob, maxOverlap, True],
                                                         list(genepairs), 3, pac.partition_inputs, {0: pac.combine_dictionaries},
                                                         number=parallel_compute_number,
                                                         procnumber=parallel_compute_number)
        mgenedict = edg.get_gene_dict(mpairsdict)

        # old_mpairsdict, old_mgenedict = mutexpairs(numCases, geneToCases, patientToGenes, genepairs, p=mprob, maxOverlap=maxOverlap)
        #
        # print len(set(mpairsdict.keys()).difference(set(old_mpairsdict.keys())))
        # print len(set(old_mpairsdict.keys()).difference(set(mpairsdict.keys())))

    else:
        mpairsdict, mgenedict = mutexpairs(numCases, geneToCases, patientToGenes, genepairs, p=mprob, maxOverlap=maxOverlap)




    if filter_mutex_gain_loss:
        prevlen = len(mpairsdict)
        mpairsdict, mgenedict = filter_mutex_gainloss(mpairsdict, mgenedict)
        print prevlen - len(mpairsdict), " mutex pairs removed by same-gene filter."

    if filter_mutex_same_segment:
        prevlen = len(mpairsdict)
        mpairsdict, mgenedict = filter_mutex_samesegment(mpairsdict, mgenedict)
        print prevlen - len(mpairsdict), " mutex pairs removed by same-segment mutual exclusivity filter."

    print len(mpairsdict), " mutex pairs found in ", time.time() - t

    return mpairsdict, mgenedict


def calc_Network(pairsdict, genesdict, geneToCases, local_edge_bet, pair_filename, gene_filename, weight=None):
    pairsdict_Network = pairsdict.copy()
    genesdict_Network = {}
    for gene in genesdict:
        genesdict_Network[gene] = {}

    G = nx.Graph()
    for pair in pairsdict_Network:
        if weight:
            G.add_edge(*tuple(pair), weight=pairsdict_Network[pair][weight])
            # print "weight is ", pairsdict_Network[pair][weight]
        else:
            G.add_edge(*tuple(pair))

    t0 = time.time()
    # print G
    # print
    node_bet_cent = nx.algorithms.centrality.betweenness_centrality(G, weight='weight')
    node_close_cent = nx.algorithms.centrality.closeness_centrality(G, distance='weight')


    node_clust_coef = nx.algorithms.clustering(G, weight='weight')
    edge_bet_cent = nx.algorithms.centrality.edge_betweenness_centrality(G, weight='weight')
    print "Network attributes calc'ed in ", time.time() - t0


    for gene in genesdict_Network:
        genesdict_Network[gene]['Gene'] = gene
        genesdict_Network[gene]['MutationFrequency'] = len(geneToCases[gene])
        genesdict_Network[gene]['NumPairs'] = len(genesdict[gene] if gene in genesdict else [])
        genesdict_Network[gene]['NodeBetweennessCentrality'] = node_bet_cent[gene]
        genesdict_Network[gene]['NodeClosenessCentrality'] = node_close_cent[gene]
        genesdict_Network[gene]['NodeClusteringCoefficient'] = node_clust_coef[gene]


    numpairs = 0
    t = time.time()

    for pair in pairsdict_Network:

        gene0, gene1 = tuple(pair)
        # print "Edge bet cent", edge_bet_cent
        # print "Pair is ", pair
        # print "Tuple pair is ", tuple(pair)
        #print edge_bet_cent[tuple(pair)]

        try :
            pairsdict_Network[pair]['TotalEdgeBetweennessCentrality'] = edge_bet_cent[(gene1, gene0)]
        except KeyError:
            pairsdict_Network[pair]['TotalEdgeBetweennessCentrality'] = edge_bet_cent[(gene0, gene1)]
        pairsdict_Network[pair]['AverageNodeClosenessCentrality'] = (node_close_cent[gene0] + node_close_cent[gene1]) / 2.0
        pairsdict_Network[pair]['AverageNodeClusteringCoefficient'] = (node_clust_coef[gene0] + node_clust_coef[gene1]) / 2.0
        # pairsdict_Network[pair]['HarmonicNodeClosenessCentrality'] = math.sqrt(node_close_cent[gene0] * node_close_cent[gene1]))
        # pairsdict_Network[pair]['NumTriplets'] = len(nx.common_neighbors(G, *tuple(pair)))
        # Local Edge Betweenness Centrality

        if local_edge_bet:
            # Get subgraph
            t1 = time.time()
            neighbors0 = set([neighbor for neighbor in G.neighbors(gene0)])
            neighbors1 = set([neighbor for neighbor in G.neighbors(gene1)])
            neighbors = set.union(neighbors0, neighbors1)
            S = G.subgraph(neighbors)
            # get the edge betweenness centrality
            local_edge_bet_cent = nx.algorithms.centrality.edge_betweenness_centrality(S)
            try:
                pairsdict_Network[pair]['LocalEdgeBetweennessCentrality'] = local_edge_bet_cent[(gene0, gene1)]
            except KeyError:
                pairsdict_Network[pair]['LocalEdgeBetweennessCentrality'] = local_edge_bet_cent[(gene1, gene0)]
            print "Local edge calculated in ", time.time() - t1

        numpairs += 1
        # print numpairs, "Number of pairs calculated in ", time.time() - t

    writeanydict(pairsdict_Network, pair_filename)
    writeanydict(genesdict_Network, gene_filename)


    # print "Network calculated in ", time.time() - t
    print "Network pairs written to ",  pair_filename
    print "Network genes written to ", gene_filename



    # print "Network pairs written to " + file_prefix + "_Networkpairs.tsv"
    # print "Network genes written to " + file_prefix + "_Networkgenes.tsv"





def complete_cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, cprob, minCooccur,
                          cooccur_distance_threshold, min_cooccurrence_ratio, parallel_compute_number,
                          filter_cooccur_same_segment, fcss_cratiothresh, fcss_mutfreqdiffratiothresh,
                          fcss_coveragethresh, fcss_probabilitythresh):
    t1 = time.time()





    if parallel_compute_number:
        cpairsdict, cgenedict = pac.parallel_compute_new(cooccurpairs, [numCases, geneToCases, patientToGenes, genepairs, cprob, minCooccur,
                                                                        True, cooccur_distance_threshold, min_cooccurrence_ratio],
                                                         list(genepairs), 3, pac.partition_inputs, {0: pac.combine_dictionaries},
                                                         number=parallel_compute_number,
                                                         procnumber=parallel_compute_number)
        cgenedict = edg.get_gene_dict(cpairsdict)


        # old_cpairsdict, old_cgenedict = cooccurpairs(numCases, geneToCases, patientToGenes, cgenepairs, p=cprob, minCooccur=minCooccur, min_cooccurrence_ratio=min_cooccurrence_ratio)
        #
        # print len(set(cpairsdict.keys()).difference(set(old_cpairsdict.keys())))
        # print len(set(old_cpairsdict.keys()).difference(set(cpairsdict.keys())))



    else:
        cpairsdict, cgenedict = cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, p=cprob, minCooccur=minCooccur,
                                             compute_prob=True, cooccur_distance_threshold=cooccur_distance_threshold,
                                             min_cooccurrence_ratio=min_cooccurrence_ratio)

    if filter_cooccur_same_segment:
        prevlen = len(cpairsdict)
        cpairsdict, cgenedict = filter_cooccur_samesegment(cpairsdict, cgenedict, fcss_cratiothresh, fcss_mutfreqdiffratiothresh, fcss_coveragethresh,
                                                           fcss_probabilitythresh)

        print prevlen - len(cpairsdict), " cooccur pairs removed by same-segment cooccur filter"


    t2 = time.time()
    print len(cpairsdict), " cooccurring pairs found in ", t2 - t1

    return cpairsdict, cgenedict


def remove_extra_triplets(Triplets):
    new_Triplets = []

    gene_sets = set()

    for Triplet in Triplets:
        gene_set = frozenset(Triplet.genes)
        if gene_set not in gene_sets:
            gene_sets.add(gene_set)
            new_Triplets.append(Triplet)

    return new_Triplets


def sort_triplets_by_type(Triplets):
    Triplet_dict = {}

    for Triplet in Triplets:
        type = Triplet.type
        if type not in Triplet_dict:
            Triplet_dict[type] = []
        Triplet_dict[type].append(Triplet)

    return Triplet_dict



def get_parser():
    # Parse arguments
    import argparse

    description = 'Given an input number of samples n, and probability of mutual exclusivity p, ' \
                  'plots the number of mutations that each sample' \
                  'must have in order to reach that probability p.'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-m', '--mutation_matrix', default=None,
                        help='File name for mutation data.')
    parser.add_argument('-o', '--output_prefix', default=None,
                        help='Output path prefix (TSV format). Otherwise just prints.')
    parser.add_argument('-mf', '--min_freq', type=int, default=0,
                        help='Minimum gene mutation frequency.')
    parser.add_argument('-mperc', '--min_percentile', type=float, default=0)

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
    parser.add_argument('-fcss_rt', '--fcss_cratiothresh', type=float, default=0.9, help='Maximum ratio for cooccurring pairs,'
                                                                                         'which are probably the same at this point')
    parser.add_argument('-fcss_mfdrt', '--fcss_mutfreqdiffratiothresh', type=float, default=1, help='Maximum mutation frequency'
                                                                                                       'difference ratio for cooccurring pairs,'
                                                                                     'which are probably the same at this point')
    parser.add_argument('-fcss_mfdt', '--fcss_mutfreqdiffthresh', type=float, default=20, help='Maximum mutation frequency difference.')
    parser.add_argument('-fcss_ct', '--fcss_coveragethresh', type=float, default=50, help='Minimum coverage for cooccurring pairs,'
                                                                                     'which are probably the same at this point')
    parser.add_argument('-fcss_pt', '--fcss_probabilitythresh', type=float, default=1e-20, help='Maximum ratio for cooccurring pairs,'
                                                                             'which are probably the same at this point')

    parser.add_argument('-bgbp', '--bin_genes_by_pairs', default=False, help='Bin the genes by the pairs in a second iteration.')
    parser.add_argument('-rpt', '--ridiculous_pvalue_thresh', type=float, default=0, help="Filter ridiculous pvalues")


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
    parser.add_argument('-ctsc', '--calc_triple_scores', default=True, help='Calculate the set score, etc. of Triplet.')
    parser.add_argument('-om', '--only_mixed', default=False, help="Only search for mixed triplets.")

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
    parser.add_argument('-plf', '--pair_list_file', default=None, help='File from which to load pair lists.')



    return parser



def main():
    run(get_parser().parse_args(sys.argv[1:]))

def run(args):


    tstart = time.time()


    # ------------------------------------------------------------------------------------------------------------
    # Turn command line arguments to shorter variable handles
    # ------------------------------------------------------------------------------------------------------------
    mdictfile = args.mdictfile
    cdictfile = args.cdictfile

    mprob = args.mprob
    cprob = args.cprob
    cooccur_distance_threshold = args.cooccur_distance_threshold
    mutationmatrix = args.mutation_matrix
    file_prefix = args.output_prefix
    geneFile = args.gene_file
    patientFile = args.patient_file
    gene_blacklist = args.gene_blacklist_file
    patient_blacklist = args.patient_blacklist_file
    gene_file_1 = args.gene_list_1
    gene_file_2 = args.gene_list_2

    pair_list_file = args.pair_list_file

    minFreq = args.min_freq
    maxOverlap = args.max_overlaps
    minCooccur = args.min_cooccur
    min_cooccurrence_ratio = args.min_cooccurrence_ratio
    min_percentile = args.min_percentile
    top_percentile = args.top_percentile
    top_number = args.top_number

    parallel_compute_number = args.parallel_compute_number

    filter_coding = args.filter_coding
    filter_mutex_gain_loss = args.filter_mutex_gain_loss
    filter_mutex_same_segment = args.filter_mutex_same_segment

    filter_cooccur_same_segment = args.filter_cooccur_same_segment
    fcss_cratiothresh = args.fcss_cratiothresh
    fcss_mutfreqdiffthresh = args.fcss_mutfreqdiffthresh
    fcss_mutfreqdiffratiothresh = args.fcss_mutfreqdiffratiothresh
    fcss_coveragethresh = args.fcss_coveragethresh
    fcss_probabilitythresh = args.fcss_probabilitythresh
    ridiculous_pvalue_thresh = args.ridiculous_pvalue_thresh

    filter_cooccur = args.filter_cooccur
    use_downstream = args.use_downstream
    scoreperc = args.scoreperc
    ratioperc = args.ratioperc
    coverageperc = args.coverageperc
    combscoreperc = args.combscoreperc


    calc_triple_stats = args.calc_triple_stats
    calc_triple_scores = args.calc_triple_scores
    only_mixed = args.only_mixed

    group_type = args.group_type
    local_edge_bet = args.local_edge_bet
    pair_type = args.pair_type
    # ------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------







    if not mdictfile and not cdictfile and not mutationmatrix:
        print "No input file specified"
        return








<<<<<<< Updated upstream

    # LOAD MUTATION DATA
    # --------------------------------

=======
    # Load mutation data
>>>>>>> Stashed changes
    if mutationmatrix:

        mutations = mex.remove_blacklists(gene_blacklist, patient_blacklist,
                                  *mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq))
        numGenes, numCases, genes, patients, geneToCases, patientToGenes = mutations

        print 'Pre-percentile filter Mutation data: %s genes x %s patients' % (numGenes, numCases)

        if min_percentile != 0:
            print "Loading using minimum percentile: ", min_percentile
            minFreq = int(min_percentile * 0.01 * len(patients))
            print "Minimum mutation frequency: ", minFreq
            mutations = mex.remove_blacklists(gene_blacklist, patient_blacklist,
                                  *mex.load_mutation_data(mutationmatrix, patientFile, geneFile, minFreq))
            numGenes, numCases, genes, patients, geneToCases, patientToGenes = mutations
            print 'Post-top-percentile filter Mutation data: %s genes x %s patients' % (numGenes, numCases)

        if top_number:
            numGenes, numCases, genes, patients, geneToCases, patientToGenes = mex.filtertopnumbermatrix(top_number, *mutations)


        print "Post-filter Mutation data: %s genes x %s patients" % (numGenes, numCases)


    else:
        print "No mutation matrix specified"
        return





    # LOAD OR CALCULATE PAIRS
    if mdictfile or cdictfile:
        if mdictfile:
            medges = edg.EdgeReader(mdictfile, include_gene_cols=True)
            mpairsdict = medges.get_pair_dict()
            mgenedict = medges.get_gene_dict()
            if not file_prefix:
                file_prefix = mdictfile[:-4]

            writegenedict(mgenedict, file_prefix + 'Mgenes.tsv')

        if cdictfile:
            cedges = edg.EdgeReader(cdictfile, include_gene_cols=True)
            cpairsdict = cedges.get_pair_dict()
            cgenedict = cedges.get_gene_dict()
            if not file_prefix:
                file_prefix = cdictfile[:-4]

    else:

        if not file_prefix:
            file_prefix = mutationmatrix

        file_prefix += '.g' + str(numGenes) + '.p' + str(numCases) + ('.mf' + str(minFreq) if minFreq else '') + '.mp' + str(mprob) + '.mo' + str(maxOverlap) + \
                      '.cp' + str(cprob) +  '.mc' + str(minCooccur) + '.mcr' + str(min_cooccurrence_ratio) + ('.cdt' + str(cooccur_distance_threshold) if cooccur_distance_threshold else '') \
                       + '.pcn' + ('0' if not parallel_compute_number else str(parallel_compute_number))



        print "Getting gene pairs to test..."

        if gene_file_1:
            with open(gene_file_1) as f:
                gene_list_1 = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
        else:
            gene_list_1 = None

        if gene_file_2:
            with open(gene_file_2) as f:
                gene_list_2 = set(l.rstrip().split()[0] for l in f if not l.startswith("#"))
        else:
            gene_list_2 = None



        if gene_list_1 and gene_list_2:
            genepairs = getgenepairs(geneToCases, genes1=gene_list_1, genes2=gene_list_2)
        elif gene_list_1:
            genepairs = getgenepairs(geneToCases, genes1=gene_list_1)
        elif pair_list_file:
            genepairs = loadgenepairs(pair_list_file, geneToCases)
        else:
            genepairs = getgenepairs(geneToCases, genes1=genes)

        print "Gene pairs finished"

<<<<<<< Updated upstream








        # MUTEX BEGINNING HERE
=======
>>>>>>> Stashed changes




        # MUTEX BEGINNING HERE


        if cooccur_distance_threshold:
            print "Before distance threshold: ", len(genepairs), " pairs to test"
            genepairs = filter_pair_distances(genepairs, cooccur_distance_threshold)
            print "After distance threshold: ", len(genepairs), " pairs to test"
            cooccur_distance_threshold = None


        cpairsdict, cgenedict = complete_cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, cprob, minCooccur,
                          cooccur_distance_threshold, min_cooccurrence_ratio, parallel_compute_number,
                          filter_cooccur_same_segment, fcss_cratiothresh, fcss_mutfreqdiffratiothresh,
                          fcss_coveragethresh, fcss_probabilitythresh)

        mpairsdict, mgenedict = complete_mutexpairs(numCases, geneToCases, patientToGenes, genepairs, mprob, maxOverlap, parallel_compute_number,
                filter_mutex_gain_loss, filter_mutex_same_segment)


        if ridiculous_pvalue_thresh:
            entryToFilter = {}

            entryToFilter['Probability'] = lambda prob: prob >= ridiculous_pvalue_thresh

            mpairsdict, mgenedict = filterpairs_new(mpairsdict, mgenedict, entryToFilter)
            cpairsdict, cgenedict = filterpairs_new(cpairsdict, cgenedict, entryToFilter)
            print "After filtering ridiculous pvalues below ", ridiculous_pvalue_thresh
            print len(cpairsdict), "cooccur pairs remain"
            print len(mpairsdict), "mutex pairs remain"



        file_prefix += '.mn' + str(len(mpairsdict)) + '.cn' + str(len(cpairsdict))

        if mpairsdict:

            writeanydict(mpairsdict, file_prefix + 'Mpairs.tsv')
            writegenedict(mgenedict, file_prefix + 'Mgenes.tsv')
            print len(mpairsdict), "Mutex pairs written to ", file_prefix + 'Mpairs.tsv'
        if cpairsdict:
            writeanydict(cpairsdict, file_prefix + 'Cpairs.tsv')
            writegenedict(cgenedict, file_prefix + 'Cgenes.tsv')
            print len(cpairsdict), 'Cooccur pairs written to ', file_prefix + 'Cpairs.tsv'




    if filter_cooccur:
        entryToFilter = {}

        if scoreperc:
            score_thresh = stats.scoreatpercentile([cpairstat['SetScore'] for cpairstat in cpairsdict.values()], scoreperc)
            entryToFilter['SetScore'] = lambda score: score <= score_thresh

        if combscoreperc:
            combscore_thresh = stats.scoreatpercentile([cpairstat['CombinedScore'] for cpairstat in cpairsdict.values()], combscoreperc)
            entryToFilter['CombinedScore'] = lambda score: score <= combscore_thresh

        if ratioperc:
            ratio_thresh = stats.scoreatpercentile([cpairstat['CooccurrenceRatio'] for cpairstat in cpairsdict.values()], ratioperc)
            entryToFilter['CooccurrenceRatio'] = lambda score: score >= ratio_thresh

        if coverageperc:
            coverage_thresh = stats.scoreatpercentile([cpairstat['Coverage'] for cpairstat in cpairsdict.values()], coverageperc)
            entryToFilter['Coverage'] = lambda score: score >= coverage_thresh

        cpairsdict_filtered, cgenedict_filtered = filterpairs_new(cpairsdict, cgenedict, entryToFilter)

        filtered_suffix = '_' + '.sp' + str(scoreperc) + '.csp' + str(combscoreperc) + '.rp' + str(ratioperc) + '.cp' + str(coverageperc) + '.' + \
                          str(len(cpairsdict_filtered)) + 'Cpairs_filtered'

        writeanydict(cpairsdict_filtered, file_prefix + filtered_suffix + '.tsv')

        print 'Filtered cooccurring pairs written to ', file_prefix + filtered_suffix + '.tsv'

        if use_downstream:
            cpairsdict = cpairsdict_filtered
            cgenedict = cgenedict_filtered
            file_prefix += filtered_suffix


    if pair_type == 'm':
        pairsdict = mpairsdict.copy()
        genesdict = mgenedict.copy()
    elif pair_type == 'c':
        pairsdict = cpairsdict.copy()
        genesdict = cgenedict.copy()
    else:
        pairsdict = cpairsdict.copy()
        pairsdict.update(mpairsdict)

        genesdict = cgenedict.copy()
        for gene in mgenedict:
            if gene in genesdict:
                genesdict[gene] = genesdict[gene].union(mgenedict[gene])
            else:
                genesdict[gene] = mgenedict[gene]

        print "Num mutex ", len(mpairsdict)
        print "Num cooccur ", len(cpairsdict)

    print "Number of pairs total ", len(pairsdict)





    if group_type == 'Network' or group_type == 'TripletNetwork':

        if mpairsdict:
            print "------------------------------------------------------------------"
            print "MUTUAL EXCLUSIVITY NETWORK"
            calc_Network(mpairsdict, mgenedict, geneToCases, local_edge_bet, file_prefix + "_NetworkMpairs.tsv",
                 file_prefix + "_NetworkMgenes.tsv")
            print "MUTUAL EXCLUSIVITY LOGPCOV WEIGHT"
            calc_Network(mpairsdict, mgenedict, geneToCases, local_edge_bet, file_prefix + "_NetworkRoundedLogPCovMpairs.tsv",
                 file_prefix + "_NetworkRoundedLogPCovMgenes.tsv", weight='RoundedLogPCov')
        if cpairsdict:
            print "------------------------------------------------------------------"
            print "COOCCURRING NETWORK"
            calc_Network(cpairsdict, cgenedict, geneToCases, local_edge_bet, file_prefix + "_NetworkCpairs.tsv",
                 file_prefix + "_NetworkCgenes.tsv")
            print "COOCCURRING NETWORK COMBINED SCORE WEIGHT"
            calc_Network(cpairsdict, cgenedict, geneToCases, local_edge_bet, file_prefix + "_NetworkCombScoreCpairs.tsv",
                 file_prefix + "_NetworkCombScoreCgenes.tsv", weight="CombinedScore")
            print "COOCCURRING NETWORK LOGPCOV WEIGHT"
            calc_Network(cpairsdict, cgenedict, geneToCases, local_edge_bet, file_prefix + "_NetworkRoundedLogPCovCpairs.tsv",
                 file_prefix + "_NetworkRoundedLogPCovCgenes.tsv", weight='RoundedLogPCov')
        if pairsdict:
            print "------------------------------------------------------------------"
            print "COMBINED NETWORK"
            calc_Network(pairsdict, genesdict, geneToCases, local_edge_bet, file_prefix + "_Networkpairs.tsv",
                         file_prefix + "_Networkgenes.tsv")
            print "COMBINED NETWORK LOGPCOV WEIGHT"
            calc_Network(pairsdict, genesdict, geneToCases, local_edge_bet, file_prefix + "_NetworkRoundedLogPCovpairs.tsv",
             file_prefix + "_NetworkRoundedLogPCovgenes.tsv", weight='RoundedLogPCov')
        print "------------------------------------------------------------------"



    if group_type == 'Triplet' or group_type == 'TripletNetwork':

        print "BEGINNING TRIPLETS"

        name = "Triplets"
        t3 = time.time()


        if parallel_compute_number:

            Triplets, pairsdict_Triplets, sorted_pairs, genesdict_Triplets, sorted_genes = pac.parallel_compute_new(getTriplets, [pairsdict, genesdict, pairsdict.keys(), numCases, geneToCases, patientToGenes, False, 'Triplets',
                calc_triple_stats, calc_triple_scores, only_mixed], pairsdict.keys(), 2, pac.partition_inputs, {0: pac.combine_lists,
                                                                                             1: combine_pairsdict_Triplets,
                                                                                            2: pac.combine_lists,
                                                                                            3: combine_genesdict_Triplets,
                                                                                            4: pac.combine_lists}, number=parallel_compute_number,
                                                                                            procnumber=parallel_compute_number)

            Triplets = remove_extra_triplets(Triplets)


        else:
            Triplets, pairsdict_Triplets, sorted_pairs, genesdict_Triplets, sorted_genes = getTriplets(pairsdict, genesdict, pairsdict.keys(), numCases=numCases, geneToCases=geneToCases, patientToGenes = patientToGenes, compute_prob=False, name='Triplets',
                calc_triple_stats=calc_triple_stats, calc_triple_scores=calc_triple_scores, only_mixed=only_mixed)

        print len(Triplets), " triplets found in", time.time() - t3

        Triplet_dict = sort_triplets_by_type(Triplets)

        for type in Triplet_dict:
            print len(Triplet_dict[type]), " of type ",  type
            writeTriplets(Triplet_dict[type], file_prefix + '.n' + str(len(Triplet_dict[type])) + type)

        # file_prefix += '.n' + str(len(Triplets))


        pairs_header = pairsdict_Triplets.values()[0].keys()
        pairs_header.remove(name)
        genes_header = genesdict_Triplets.values()[0].keys()
        genes_header.remove(name)

        writeanydict(pairsdict_Triplets, file_prefix + "_Tripletspairs.tsv", orderedkeys=sorted_pairs, fieldnames=pairs_header)
        writeanydict(genesdict_Triplets, file_prefix + "_Tripletsgenes.tsv", orderedkeys=sorted_genes, fieldnames=genes_header)

        print "Triplet pairs written to " + file_prefix + "_Tripletspairs.tsv"
        print "Triplet genes written to " + file_prefix + "_Tripletsgenes.tsv"


    print "Time used ", time.time() - tstart

if __name__ == '__main__':
    main()









































































































# ----------------------------------------------------------------------------------------------------------------------
# BELOW FOLLOWS UNUSED CODE
# ----------------------------------------------------------------------------------------------------------------------




















































































































































































def getTriangles(mpairsdict, mgenedict, cpairsdict, cgenedict, name="Triangles"):
    """
    :return: Triangles:
        A list of all initialized triangles.
        mpairsdict_Triangles:
        All the mpairs with Triangles assigned.
        cpairsdict_Triangles:
        All the cpairs with Triangles assigned.
    """
    # Get mutexpairs, mutexdict

    # mpairsdict, mgenedict = mutexpairs(numCases, geneToCases, patientToGenes, genepairs, p=mprob, maxOverlap=maxOverlap)
    # cpairsdict, cgenedict = cooccurpairs(numCases, geneToCases, patientToGenes, genepairs, p=cprob, minCooccur=minCooccur)

    Triangles = []
    num_Triangles = 0
    mpairsdict_Triangles = mpairsdict.copy()
    for mpair in mpairsdict:
        mpairsdict_Triangles[mpair][name] = set()

    cpairsdict_Triangles = cpairsdict.copy()
    for cpair in cpairsdict:
        cpairsdict_Triangles[cpair][name] = set()


    for mpair in mpairsdict:
        gene1, gene2 = tuple(mpair)
        if gene1 in cgenedict and gene2 in cgenedict:
            cgenes1 = cgenedict[gene1]
            cgenes2 = cgenedict[gene2]

            # The set of genes that cooccur with both genes.
            cboth = cgenes1.intersection(cgenes2)
            if cboth:
                for cbothgene in cboth:
                    # The cooccurring pairs of the triangle.
                    cbothpair1 = frozenset([gene1, cbothgene])
                    cbothpair2 = frozenset([gene2, cbothgene])


                    mpairdict = {}
                    mpairdict[mpair] = mpairsdict[mpair]

                    cpairdict = {}
                    cpairdict[cbothpair1] = cpairsdict[cbothpair1]
                    cpairdict[cbothpair2] = cpairsdict[cbothpair2]

                    newTriangle = Triangle(num_Triangles, mpairdict=mpairdict, cpairdict=cpairdict)
                    Triangles.append(newTriangle)

                    mpairsdict_Triangles[mpair][name].add(num_Triangles)
                    cpairsdict_Triangles[cbothpair1][name].add(num_Triangles)
                    cpairsdict_Triangles[cbothpair2][name].add(num_Triangles)


                    # Add to the mpairsdict_Triangles and cpairsdict_Triangles
                    # if mpair not in mpairsdict_Triangles:
                    #     mpairsdict_Triangles[mpair] = mpairsdict[mpair]
                    #     mpairsdict_Triangles[mpair][name] = {num_Triangles}
                    # else:
                    #     mpairsdict_Triangles[mpair][name].add(num_Triangles)
                    #
                    # if cbothpair1 not in cpairsdict_Triangles:
                    #     cpairsdict_Triangles[cbothpair1] = cpairsdict[cbothpair1]
                    #     cpairsdict_Triangles[cbothpair1][name] = {num_Triangles}
                    # else:
                    #     cpairsdict_Triangles[cbothpair1][name].add(num_Triangles)
                    #
                    # if cbothpair2 not in cpairsdict_Triangles:
                    #     cpairsdict_Triangles[cbothpair2] = cpairsdict[cbothpair2]
                    #     cpairsdict_Triangles[cbothpair2][name] = {num_Triangles}
                    # else:
                    #     cpairsdict_Triangles[cbothpair2][name].add(num_Triangles)


                    num_Triangles += 1

    sorted_mpairs = sorted(mpairsdict_Triangles.keys(), key=lambda entry: len(mpairsdict_Triangles[entry]), reverse=True)
    sorted_cpairs = sorted(mpairsdict_Triangles.keys(), key=lambda entry: len(mpairsdict_Triangles[entry]), reverse=True)

    return Triangles, mpairsdict_Triangles, cpairsdict_Triangles, sorted_mpairs, sorted_cpairs






def writeTriangles(Triangles, file_prefix, delimiter='\t'):
    # Dictwriter
    # Get writabledict from each triangle, make header from
    Trianglefile = file_prefix + '_triangles.tsv'
    genefile = file_prefix + '_genes.txt'

    with open(Trianglefile, 'w') as csvfile:
        fieldnames = Triangles[0].getwritabledict().keys()
        writer = csv.DictWriter(csvfile, delimiter=delimiter, fieldnames=fieldnames)

        writer.writeheader()
        for Triangle in Triangles:
            writer.writerow(Triangle.getwritabledict())

    with open(genefile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        # writer.writerow(['Genes'])
        for Triangle in Triangles:
            nodes = Triangle.nodes
            for node in nodes:
                writer.writerow([node])


class Triangle:
    def __init__(self, number, mpairdict=None, cpairdict=None):
        assert len(mpairdict) == 1
        assert len(cpairdict) == 2

        self.number = number

        self.mpairdict = {}
        # Turn into our own sets.
        for inputmpair in mpairdict:
            # mpair = set(inputmpair)
            self.mpairdict[inputmpair] = mpairdict[inputmpair]


        self.cpairdict = {}
        for inputcpair in cpairdict:
            # cpair = set(inputcpair)
            self.cpairdict[inputcpair] = cpairdict[inputcpair]

        self.nodes = set.union(*[set(cpair) for cpair in cpairdict.keys()])
        self.medge = self.mpairdict.keys()[0]
        self.vertex = next(iter(self.nodes.difference(self.medge)))


        self.statsdict = {}




    def getwritabledict(self):
        wdict = collections.OrderedDict()
        wdict['Vertex'] = self.vertex
        wdict['Edge1'], wdict['Edge2'] = tuple(self.medge)
        return wdict
        # Later possibly add mutation frequencies, gene information, probabilities, etc.

    #def getGOwritabledict(self):
        # when the time comes, look up GO annotation to optimize how convenient this is.

    def getmpairs(self):
        return self.mpairdict

    def getcpairs(self):
        return self.cpairdict


class CoMPair:
    def __init__(self, mpairdict=None):
        assert len(mpairdict) == 2
        self.mpairdict = mpairdict.copy()
        self.genes = set.union(*[set(mpair) for mpair in mpairdict.keys()])
    def getwritabledict(self):
        [mpair1, mpair2] = self.mpairdict.keys()
        wdict = collections.OrderedDict()
        wdict['AGene1'], wdict['AGene2'] = tuple(mpair1)
        wdict['BGene1'], wdict['BGene2'] = tuple(mpair2)
        return wdict

def getCoMPairs(mpairsdict, mgenedict, cpairsdict, cgenedict):
    """

    :return: CoMPairs:
        A list of all initialized CoMPairs.
    """
    # Get mutexpairs, mutexdict

    CoMPairs = []
    for mpair in mpairsdict:
        gene1, gene2 = tuple(mpair)
        if gene1 in cgenedict and gene2 in cgenedict:
            cgenes1 = cgenedict[gene1]
            cgenes2 = cgenedict[gene2]

            # The set of genes that cooccur with either in the pair.
            ceither = cgenes1.union(cgenes2)

            if ceither:
                for mpair2 in mpairsdict:
                    mpair2_unfrozen = iter(mpair2)
                    if mpair2_unfrozen.next() in ceither and mpair2_unfrozen.next() in ceither:
                        mpairdict = {}
                        mpairdict[mpair] = mpairsdict[mpair]
                        mpairdict[mpair2] = mpairsdict[mpair2]
                        newCoMPair = CoMPair(mpairdict=mpairdict)
                        CoMPairs.append(newCoMPair)

    return CoMPairs


def writeCoMPairs(CoMPairs, file_prefix, delimiter='\t'):
    # Dictwriter
    # Get writabledict from each CoMPair, make header from
    CoMPairfile = file_prefix + '_CoMPairs.tsv'
    genefile = file_prefix + '_genes.txt'

    with open(CoMPairfile, 'w') as csvfile:
        fieldnames = CoMPairs[0].getwritabledict().keys()
        writer = csv.DictWriter(csvfile, delimiter=delimiter, fieldnames=fieldnames)

        writer.writeheader()
        for CoMPair in CoMPairs:
            writer.writerow(CoMPair.getwritabledict())

    with open(genefile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        # writer.writerow(['Genes'])
        for CoMPair in CoMPairs:
            genes = CoMPair.genes
            for gene in genes:
                writer.writerow([gene])


def filter_pair_distances(genepairs, cooccur_distance_threshold):
    newgenepairs = genepairs[:]

    for genepair in genepairs:
        gene0, gene1 = tuple(genepair)
        if gene0[-4:] in {'loss', 'gain'} and gene1[-4:] in {'loss', 'gain'}:
            if bgbp.get_gene_distance(*genepair) and (bgbp.get_gene_distance(*genepair) < cooccur_distance_threshold):
                # print genepair
                newgenepairs.remove(genepair)
    return newgenepairs



def filter_pair_coding(genepairs, filename='/Users/jlu96/conte/jlu/geneToLength_all_firstseen.txt'):
    # Load the new genes.
    genes = set()
    with open(filename, 'rU') as f:
        arrs = [l.rstrip().split("\t") for l in f if not l.startswith("#")]
        for arr in arrs:
            gene = arr[0]
            genes.add(gene)

    newgenepairs = genepairs.copy()
    for genepair in genepairs:
        gene0, gene1 = tuple(genepair)
        if gene0 not in genes or gene1 not in genes:
            newgenepairs.remove(genepair)

    return newgenepairs


def filter_mutex_gainloss(mpairdict, mgenedict):
    newmpairdict = mpairdict.copy()
    newmgenedict = mgenedict.copy()
    for pair in mpairdict:
        genes = tuple(pair)
        if genes[0][:-4] == genes[1][:-4]:
            newmpairdict.pop(pair)
            newmgenedict[genes[0]].remove(genes[1])
            if not newmgenedict[genes[0]]:
                newmgenedict.pop(genes[0])
            newmgenedict[genes[1]].remove(genes[0])
            if not newmgenedict[genes[1]]:
                newmgenedict.pop(genes[1])
    return newmpairdict, newmgenedict

def filter_mutex_samesegment(mpairdict, mgenedict):
    newmpairdict = mpairdict.copy()
    newmgenedict = mgenedict.copy()
    for pair in mpairdict:
        genes = tuple(pair)
        if not mpairdict[pair]['Concordance']:
            newmpairdict.pop(pair)
            newmgenedict[genes[0]].remove(genes[1])
            if not newmgenedict[genes[0]]:
                newmgenedict.pop(genes[0])
            newmgenedict[genes[1]].remove(genes[0])
            if not newmgenedict[genes[1]]:
                newmgenedict.pop(genes[1])
    return newmpairdict, newmgenedict

def filter_cooccur_samesegment(cpairdict, cgenedict, cratiothresh, mutfreqdiffratiothresh, coveragethresh, probabilitythresh):
    newcpairdict = cpairdict.copy()
    newcgenedict = cgenedict.copy()
    for pair in cpairdict:
        genes = tuple(pair)
        if cpairdict[pair]['Concordance']:
            if cpairdict[pair]['Probability'] < probabilitythresh \
                and cpairdict[pair]['CooccurrenceRatio'] > cratiothresh and \
                cpairdict[pair]['MutationFrequencyDifferenceRatio'] < mutfreqdiffratiothresh \
                and cpairdict[pair]['Coverage'] > coveragethresh:
                newcpairdict.pop(pair)
                newcgenedict[genes[0]].remove(genes[1])
                if not newcgenedict[genes[0]]:
                    newcgenedict.pop(genes[0])
                newcgenedict[genes[1]].remove(genes[0])
                if not newcgenedict[genes[1]]:
                    newcgenedict.pop(genes[1])
    return newcpairdict, newcgenedict


def bin_genes_cooccur_same_segment(pairdict, geneToCases, patientToGenes, cratiothresh, mutfreqdiffthresh,
                                   mutfreqdiffratiothresh, coveragethresh, probabilitythresh,
              bincol='RoundedLogPCov', filename=None):

    print "Binning genes by pairs...."

    t = time.time()
    geneToBin = {}
    for gene in geneToCases:
        geneToBin[gene] = set([gene])

    sorted_pairs = sorted(pairdict.keys(), key=lambda pair: pairdict[pair][bincol])


    for pair in sorted_pairs:
        if pairdict[pair]['Concordance'] and pairdict[pair]['Probability'] < probabilitythresh:
            if pairdict[pair]['CooccurrenceRatio'] >= cratiothresh \
                and pairdict[pair]['MutationFrequencyDifference'] <= mutfreqdiffthresh and \
                pairdict[pair]['MutationFrequencyDifferenceRatio'] < mutfreqdiffratiothresh \
                and pairdict[pair]['Coverage'] > coveragethresh:
                # print "Mutation Frequencies: ", pairdict[pair]['MutationFrequencies']
                # print "Cooccurrence Ratio: ", pairdict[pair]['CooccurrenceRatio']
                gene0, gene1 = tuple(pair)
                geneToBin[gene0] = geneToBin[gene0].union(geneToBin[gene1])
                for gene in geneToBin[gene0]:
                    geneToBin[gene] = geneToBin[gene0]

    t1 = time.time()
    # print "Time to bin genes by pairs", t1 - t

    for gene in geneToBin:
        geneToBin[gene] = frozenset(geneToBin[gene])


    # Update GeneToCases by majority ratios
    # Iterate across values


    newGeneToCases = {}

    bin_sets = set(geneToBin.values())

    print len(geneToBin), " genes reduced to ", len(bin_sets), " bins"
    print "Ratio threshold: ", cratiothresh, ', Diff ratio threshold: ', mutfreqdiffratiothresh, ', Coverage threshold: ', coveragethresh, ', Pvalue threshold: ', probabilitythresh


    bin_setToBin = {}


    for bin_set in bin_sets:
        # Fix the bin name
        suffix = iter(bin_set).next()[-4:]
        bin_name = '_'.join(sorted([the_bin[:-4] for the_bin in bin_set])) + suffix
        bin_setToBin[bin_set] = bin_name

        newGeneToCases[bin_name] = set()

        patients = set.union(*[geneToCases[gene] for gene in bin_set])
        for patient in patients:
            # find out how many genes are held by patients
            if len(patientToGenes[patient].intersection(bin_set)) >= len(bin_set)/2:
                newGeneToCases[bin_name].add(patient)



    # Make the new patientToGenes

    newPatientToGenes = {}
    for patient in patientToGenes:
        newPatientToGenes[patient] = set()

    for gene in newGeneToCases:
        for patient in newGeneToCases[gene]:
            newPatientToGenes[patient].add(gene)

    t2 = time.time()
    #print "Time to update dictionaries", t2- t1



    # Consolidate information
    # Change geneToBin to bin names

    gene_bin_entries = {}

    for gene in geneToBin:
        bin_name = bin_setToBin[geneToBin[gene]]
        geneToBin[gene] = bin_name

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

    if filename:
        writeanydict(gene_bin_entries, filename)
        print "Gene bin entries written to ", filename

    t3 = time.time()
    # print "Time to get gene bin info ", t3 - t2

    bins = set(geneToBin.values())

    # print "Bins are ", iter(bins).next()
    # print "Bin sets are ", iter(bin_sets).next()

    return gene_bin_entries, geneToBin, bins, bin_sets, bin_setToBin, newGeneToCases, newPatientToGenes


def convert_genes_to_bins(genes, geneToBin):

    return set([geneToBin[gene] for gene in genes])


def genedict_to_bin_pairs(genedict, geneToBin, bin_sets, bin_setToBin, max_tries=3):

    binpairs = set()

    # print "Genedict keys are ", genedict.keys()

    for bin_set in bin_sets:
        bin = bin_setToBin[bin_set]
        paired_bins = []

        iter_bin = iter(bin_set)

        # Try to find the corresponding bins

        for attempt in range(min(max_tries, len(bin_set))):
            bin_in_dict = False
            random_gene = iter_bin.next()
            try:
                paired_genes = genedict[random_gene]
                paired_bins = convert_genes_to_bins(paired_genes, geneToBin)
                bin_in_dict = True
            except KeyError:
                pass
            if bin_in_dict:
                break

        for paired_bin in paired_bins:
            if paired_bin != bin:
                binpair = frozenset([bin, paired_bin])
                binpairs.add(binpair)

    return list(binpairs)



def filterpairs(pairsdict, genesdict, entry, pair_filter=None):
    """
    :param pairsdict:
    :param genesdict:
    :param entry: Entry of pairsdict.
    :param pair_filter: Function that returns True if want to keep.
    :return:
    """
    new_pairsdict = pairsdict.copy()
    new_genesdict = genesdict.copy()
    for pair in pairsdict:
        if not pair_filter(pairsdict[pair][entry]):
            gene0, gene1 = tuple(pair)

            new_pairsdict.pop(pair)

            if gene0 in new_genesdict:
                new_genesdict[gene0].remove(gene1)
                if not new_genesdict[gene0]:
                    new_genesdict.pop(gene0)
            if gene1 in new_genesdict:
                new_genesdict[gene1].remove(gene0)
                if not new_genesdict[gene1]:
                    new_genesdict.pop(gene1)
    return new_pairsdict, new_genesdict