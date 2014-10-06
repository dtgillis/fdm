__author__ = 'dtgillis'

import getopt
import sys
import os
from graph_utils import graph_constructor
import networkx as nx
import numpy as np
import matplotlib as plt
import Bio.pairwise2 as pairwise
from Bio.SubsMat.MatrixInfo import pam120
from Bio.Data import CodonTable
from Bio.Seq import translate
from large_blast import blast_wrapper
import cPickle

def print_help():
    print 'graph_compare_stage_2'


def get_gene_stats():
    query_pickle_prefix = None
    subject_pickle_prefix = None
    proj_temp_dir = None
    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "", ["proj_temp_dir=", "query_prefix=", "subject_prefix="])
    except getopt.GetoptError:
        print 'graph_compare_stage_2.py --proj_temp_dir=temp_dir ' \
              '--query_prefix=prefix_for_query --subject_prefix=prefix_for_subject'
        sys.exit(2)

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            print_help()
        elif opt in ("--proj_temp_dir"):
            proj_temp_dir = arg
        elif opt in ("--query_prefix"):
            query_pickle_prefix = arg
        elif opt in ("--subject_prefix"):
            subject_pickle_prefix = arg
        else:
            assert False, "invalid option"

    # build up and check that our input files now exist.
    assert os.path.exists(proj_temp_dir), "Temp directory does not exist."
    query_pickle_gene_dict = proj_temp_dir + os.sep + query_pickle_prefix + 'gene_dict.pickle'
    query_pickle_gene_dict_exons = proj_temp_dir + os.sep + query_pickle_prefix + 'gene_dict_exons.pickle'
    assert os.path.isfile(query_pickle_gene_dict), "Query Gene Dict non existent"
    subject_pickle_gene_dict = proj_temp_dir + os.sep + subject_pickle_prefix + 'gene_dict.pickle'
    subject_pickle_gene_dict_exons = proj_temp_dir + os.sep + subject_pickle_prefix + 'gene_dict_exons.pickle'

    # load all pickle files
    query_gene_dict = cPickle.load(query_pickle_gene_dict)
    query_gene_dict_exons = cPickle.load(query_pickle_gene_dict_exons)
    subject_gene_dict = cPickle.load(subject_pickle_gene_dict)
    subject_gene_dict_exons = cPickle.load(subject_pickle_gene_dict_exons)

    # read in blast out files
    blast_out = sorted(blast_out, key=lambda blast_line: float(blast_line[10]))

    for line in blast_out:
        print line
        human_graph_key = int(line[0].split("|")[-1])
        mouse_graph_key = int(line[1].split("|")[2].strip("graph-"))
        human_graph = human_gene_dict[human_graph_key]
        mouse_graph = mouse_gene_dict[mouse_graph_key]

        #get exon counts for both species
        human_exon_list = []
        for human_edge in human_graph.edges_iter(data=True):
            human_attributes = human_edge[2]
            if human_attributes['type'] == "exon":
                human_exon_list.append(human_edge)
        mouse_exon_list = []
        for mouse_edge in mouse_graph.edges_iter(data=True):
            mouse_attributes = mouse_edge[2]
            if mouse_attributes['type'] == "exon":
                mouse_exon_list.append(mouse_edge)

        score_mat = np.zeros((len(mouse_exon_list)*3, len(human_exon_list)))
        for j in range(len(human_exon_list)):
            human_attributes = human_exon_list[j][2]
            if human_attributes['type'] != "exon":
                continue
            human_codon = translate(human_attributes['seq']).replace('*', '')

            for i in range(len(mouse_exon_list)):
                mouse_attributes = mouse_exon_list[i][2]
                if mouse_attributes['type'] != "exon":
                    continue
                if len(mouse_attributes['seq']) > 3:
                    mouse_seq = [mouse_attributes['seq'],
                                    mouse_attributes['seq'][1:], mouse_attributes['seq'][2:]]
                else:
                    mouse_seq = [mouse_attributes['seq']*3]
                mouse_codons = [translate(seq).replace('*', '') for seq in mouse_seq]
                k = 0
                for mouse_codon in mouse_codons:
                    if mouse_codon =="" or human_codon == "":
                        pam_score=0.0
                    else:
                        pam_score = pairwise.align.globaldx(human_codon, mouse_codon, pam120, score_only=True)
                        #print pam_score
                        #print pairwise.format_alignment(*pairwise.align.globaldx(human_codon, mouse_codon, pam120, one_alignment_only=True)[0])

                    score_mat[i*3 + k, j] = pam_score
                    k += 1

        matches = np.where(score_mat == np.max(score_mat, axis=0))

        for i in range(len(human_exon_list)):
            print human_exon_list[matches[1][i]]
            mouse_exon = matches[0][i]
            mouse_exon_shift = mouse_exon % 3
            mouse_range = mouse_exon/3
            print mouse_exon_list[mouse_range][2]['seq'][mouse_exon_shift:]

if __name__ == '__main__':

    get_gene_stats()
