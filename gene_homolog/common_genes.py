__author__ = 'dtgillis'

import actgraph
import numpy as np
import pysam
import os
import swalign


class Exon():

    def __init__(self, seq, chrm, start, end, weight):

        self.seq = seq
        self.chromosome = chrm
        self.start = start
        self.end = end
        self.weight = weight



class GeneMultiGraph():

    def __init__(self, exon_list, intron_list, splice_list, start_node_list, end_node_list, novel_node_list ,
                 exon_wt_list, intron_wt_list, splice_wt_list):

        self.exon_list = exon_list
        self.intron_list = intron_list
        self.splice_list = splice_list
        self.start_node_list = start_node_list
        self.end_node_list = end_node_list
        self.novel_node_list = novel_node_list
        self.exon_wt_list = exon_wt_list
        self.intron_wt_list = intron_wt_list
        self.splice_wt_list = splice_wt_list
        self.gene = None


class Gene():
    """
    Class representing a gene and its name.
    """
    def __init__(self, name, organism, chromosome):

        self.name = name
        self.organism = organism
        self.chromosome = chromosome


class GeneDictionary():
    """
    This class will construct a homologous gene dictionary from
    http://www.informatics.jax.org/homology.shtml.
    The file will have a format similar to


    """

    def __init__(self, homolog_file):

        self.gene_homolog_file = homolog_file
        self.organism_list = []
        self.gene_dict = dict()

    def parse_gene_file(self):

        for line in open(self.gene_homolog_file, 'r'):

            fields = line.split('\t')

            if fields[0] == 'HomoloGene ID':
                continue
            #TODO change this to dict of dict?
            if fields[0] not in self.gene_dict:
                self.gene_dict[fields[0]] = []

            #create gene object
            gene = Gene(fields[3], fields[1], fields[8].split()[0])

            self.gene_dict[fields[0]].append(gene)


def levenshtein(source_string, target):
    """
    http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
    """
    if len(source_string) < len(target):
        return levenshtein(target, source_string)

    # So now we have len(source) >= len(target).
    if len(target) == 0:
        return len(source_string)

    # We call tuple() to force strings to be used as sequences
    # ('c', 'a', 't', 's') - numpy uses them as values by default.
    source_string = np.array(tuple(source_string))
    target = np.array(tuple(target))

    # We use a dynamic programming algorithm, but with the
    # added optimization that we only need the last two rows
    # of the matrix.
    previous_row = np.arange(target.size + 1)
    for s in source_string:
        # Insertion (target grows longer than source):
        current_row = previous_row + 1

        # Substitution or matching:
        # Target and source items are aligned, and either
        # are different (cost of 1), or are the same (cost of 0).
        current_row[1:] = np.minimum(
                current_row[1:],
                np.add(previous_row[:-1], target != s))

        # Deletion (target grows shorter than source):
        current_row[1:] = np.minimum(
                current_row[1:],
                current_row[0:-1] + 1)

        previous_row = current_row

    return previous_row[-1]


class GeneCompare():
    """
    Will do lots of work for us.
    """

    def __init__(self, mouse_act_file_list, human_act_file_list, gene_dict, stat_file_name, mouse_ref_dir,
                 human_ref_dir):

        self.mouse_act_file_list = mouse_act_file_list
        self.human_act_file_list = human_act_file_list
        self.mouse_act_list = []
        self.human_act_list = []
        self.gene_dict = gene_dict
        self.stat_file_name = stat_file_name
        self.create_stats_file()
        self.mouse_ref_dir = mouse_ref_dir
        self.human_ref_dir = human_ref_dir

    def load_act_graphs(self):
        """
        This function will load up the act graphs for the mouse and human
        act files passed in.
        :return:
        """
        for mouse_act_file in self.mouse_act_file_list:
            act_graph = actgraph.actFile(mouse_act_file)
            self.mouse_act_list.append(act_graph)

        for human_act_file in self.human_act_file_list:
            act_graph = actgraph.actFile(human_act_file)
            self.human_act_list.append(act_graph)

    def get_exon_splice_out(self, graph, node):

        out_degree = 0
        for splice_node in graph.splice_list:

            if splice_node[0] == node:
                out_degree += 1

        return out_degree

    def get_exon_splice_in(self, graph, node):

        in_degree = 0

        for splice_node in graph.splice_list:

            if splice_node[1] == node:
                in_degree += 1

        return in_degree

    def compare_graphs(self, human_graph, mouse_graph):

        stat_file = open(self.stat_file_name ,'a')

        # get reference for the exon sequences
        human_chrm = human_graph.gene.chromosome.lower()
        mouse_chrm = mouse_graph.gene.chromosome.lower()
        human_exons_seq = []
        mouse_exons_seq = []
        human_fasta = pysam.Fastafile(self.human_ref_dir + os.sep + human_chrm + ".fa")
        mouse_fasta = pysam.Fastafile(self.mouse_ref_dir + os.sep + mouse_chrm + ".fa")

        for exon in human_graph.exon_list:
            human_exons_seq.append(human_fasta.fetch(region="{0:s}:{1:d}-{2:d}".format(human_chrm, exon[0], exon[1])))

        for exon in mouse_graph.exon_list:
            mouse_exons_seq.append(mouse_fasta.fetch(region="{0:s}:{1:d}-{2:d}".format(mouse_chrm, exon[0], exon[1])))
        overall_levenshtein_scores = -1 * np.ones((len(human_exons_seq), len(mouse_exons_seq)))

        for i in range(len(human_exons_seq)):
            for j in range(len(mouse_exons_seq)):
                #if i <= j:
                human_string = human_exons_seq[i].lower()
                mouse_string = mouse_exons_seq[j].lower()

                max_length = float(max([len(human_string),len(mouse_string)]))
                overall_levenshtein_scores[i, j] = levenshtein(human_string, mouse_string)/max_length
                #else:
                 #   continue

        exons_of_interest = np.where(overall_levenshtein_scores <= .1)
        print exons_of_interest
        if len(exons_of_interest) != 0:
            similar_exons = []
            for i in range(len(exons_of_interest[0])):
                exon_dict = {}
                human_index = exons_of_interest[0][i]
                mouse_index = exons_of_interest[1][i]
                human_exon = Exon(seq=human_exons_seq[human_index],
                                  chrm=human_chrm, start=human_graph.exon_list[human_index][0],
                                  end=human_graph.exon_list[human_index][1], weight=human_graph.exon_wt_list[human_index])
                mouse_exon = Exon(seq=mouse_exons_seq[mouse_index],
                                  chrm=mouse_chrm, start=mouse_graph.exon_list[mouse_index][0],
                                  end=mouse_graph.exon_list[mouse_index][1], weight=mouse_graph.exon_wt_list[mouse_index])

                human_out_degree = self.get_exon_splice_out(human_graph, human_exon.end)
                human_out_degree += self.get_exon_splice_out(human_graph, human_exon.start)
                human_in_degree = self.get_exon_splice_in(human_graph, human_exon.start)
                human_in_degree += self.get_exon_splice_in(human_graph, human_exon.end)
                mouse_out_degree = self.get_exon_splice_out(mouse_graph, mouse_exon.end)
                mouse_out_degree += self.get_exon_splice_out(mouse_graph, mouse_exon.start)
                mouse_in_degree = self.get_exon_splice_in(mouse_graph, mouse_exon.start)
                mouse_in_degree += self.get_exon_splice_in(mouse_graph, mouse_exon.end)
                exon_dict['human'] = human_exon
                exon_dict['mouse'] = mouse_exon
                similar_exons.append(exon_dict)



        #
        # human_splice_wt = .0
        # mouse_splice_wt = .0
        #
        # for splice_wt in human_graph.splice_wt_list:
        #     human_splice_wt += float(splice_wt)
        # for splice_wt in mouse_graph.splice_wt_list:
        #     mouse_splice_wt += float(splice_wt)
        #
        # mouse_exon_wt = .0
        # human_exon_wt = .0
        #
        # for exon_wt in human_graph.exon_wt_list:
        #     human_exon_wt += float(exon_wt)
        # for exon_wt in mouse_graph.exon_wt_list:
        #     mouse_exon_wt += float(exon_wt)
        #
        # mouse_intron_wt = .0
        # human_intron_wt = .0
        #
        # for intron_wt in human_graph.intron_wt_list:
        #     human_intron_wt += float(intron_wt)
        # for intron_wt in mouse_graph.intron_wt_list:
        #     mouse_intron_wt += float(intron_wt)
        #
        # gene_stats = "{0:d} {1:f} {2:d} {3:f} {4:d} {5:f} {6:d} {7:f} {8:d} {9:f} {10:d} {11:f}\n".format(
        #     len(human_graph.exon_list), human_exon_wt, len(mouse_graph.exon_list), mouse_exon_wt,
        #     len(human_graph.intron_list), human_intron_wt, len(mouse_graph.intron_list), mouse_intron_wt,
        #     len(human_graph.splice_list), human_splice_wt, len(mouse_graph.splice_list), mouse_splice_wt)

        #stat_file.write(gene_stats)

    def create_stats_file(self):

        stat_file = open(self.stat_file_name, 'w')
        stat_file.write("human_exon_count human_exon_wt mouse_exon_count mouse_exon_wt "
                        "human_intron_count human_intron_wt mouse_intron_count mouse_intron_wt "
                        "human_splice_count human_splice_wt mouse_splice_count mouse_splice_wt\n")
        stat_file.close()

    def cross_wise_simple_stats(self):

        gene_dict = self.gene_dict.gene_dict

        for gene_id in gene_dict:

            mouse_gene = None
            human_gene = None

            for gene in gene_dict[gene_id]:

                if "mouse" in gene.organism:
                    mouse_gene = gene
                else:
                    human_gene = gene

            if mouse_gene is None or human_gene is None:
                continue

            for human_act in self.human_act_list:

                human_graph = GeneMultiGraph(*human_act.getwtgraphstruct(human_gene.name))
                human_graph.gene = human_gene

            for mouse_act in self.mouse_act_list:

                mouse_graph = GeneMultiGraph(*mouse_act.getwtgraphstruct(mouse_gene.name))
                mouse_graph.gene = mouse_gene

            self.compare_graphs(human_graph,  mouse_graph)

    def rebuild_index(self):

        for act_file in self.mouse_act_list:
            act_file.build_index()


































