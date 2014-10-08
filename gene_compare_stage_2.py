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
import actgraph

import cPickle

def print_help():
    print 'graph_compare_stage_2'

def create_act_file_from_nx(graph, gene_name):

    out_act_file = []

    for node in graph.nodes_iter():
        chrm = node.split(":")[0]
        loc = node.split(":")[1]
        out_act_file.append('{0:s}\tnovel\tnode\t{1:s}\t{1:s}\t0.00\t0\t.\t{2:s}\n'.format(chrm, loc, gene_name))

    for edge in graph.edges_iter(data=True):
        chrm = edge[0].split(":")[0]
        start = edge[0].split(":")[-1]
        end = edge[1].split(":")[-1]
        out_act_file.append('{0:s}\tnovel\t{1:s}\t{2:s}\t{3:s}\t{4:f}\t0\t.\t{5:s}\n'.format(
          chrm,edge[2]['type'], start, end, edge[2]['depth'], gene_name))

    return out_act_file

def get_gene_stats():
    query_pickle_prefix = None
    subject_pickle_prefix = None
    proj_temp_dir = None
    blast_out_file = None
    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "", ["proj_temp_dir=", "query_prefix=", "subject_prefix=", "sorted_blastn_out="])
    except getopt.GetoptError:
        print 'gene_compare_stage_2.py --proj_temp_dir=temp_dir ' \
              '--query_prefix=prefix_for_query --subject_prefix=prefix_for_subject --sorted_blastn_out'
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
        elif opt in ("--sorted_blastn_out"):
            blast_out_file = arg
        else:
            assert False, "invalid option"

    # build up and check that our input files now exist.
    assert os.path.exists(proj_temp_dir), "Temp directory does not exist."
    query_pickle_gene_dict = proj_temp_dir + os.sep + query_pickle_prefix + '.gene_dict.pickle'
    query_pickle_gene_dict_exons = proj_temp_dir + os.sep + query_pickle_prefix + '.gene_dict_exons.pickle'
    assert os.path.isfile(query_pickle_gene_dict), "Query Gene Dict non existent"
    subject_pickle_gene_dict = proj_temp_dir + os.sep + subject_pickle_prefix + '.gene_dict.pickle'
    subject_pickle_gene_dict_exons = proj_temp_dir + os.sep + subject_pickle_prefix + '.gene_dict_exons.pickle'

    # load all pickle files
    query_gene_dict = cPickle.load(open(query_pickle_gene_dict, 'rb'))
    query_gene_dict_exons = cPickle.load(open(query_pickle_gene_dict_exons, 'rb'))
    subject_gene_dict = cPickle.load(open(subject_pickle_gene_dict, 'rb'))
    subject_gene_dict_exons = cPickle.load(open(subject_pickle_gene_dict_exons, 'rb'))

    # read in blast out files
    blast_results = [line.split('\t') for line in open(blast_out_file, 'r').readlines()]
    blast_results_sorted = sorted(blast_results, key=lambda blast_line: float(blast_line[10]))

    matches_list = []
    for line in blast_results_sorted:
        if float(line[10]) < .05:
            query_graph = line[0].split('|')[2].split('-')[-1]
            subject_graph = line[1].split('|')[-1]
            matches_list.append((query_graph, subject_graph))

    print matches_list

    for match in matches_list:

        query_graph = query_gene_dict[int(match[0])]
        subject_graph = subject_gene_dict[int(match[1])]
        query_fasta_list = []
        subject_fasta_list = []
        subject_exon_list = []
        query_exon_list = []
        for subject_edge in subject_graph.edges_iter(data=True):
            if subject_edge[2]['type'] == 'exon':
                subject_exon_list.append(subject_edge)
                attributes = subject_edge[2]
                subject_fasta_list.append('>subject|{0:s}|{1:d}\n'.format(attributes['gene'], len(subject_exon_list)-1))
                subject_fasta_list.append('{0:s}\n'.format(attributes['seq']))

        subject_graph.graph['exon_list'] = subject_exon_list
        for query_edge in query_graph.edges_iter(data=True):

            if query_edge[2]['type'] == 'exon':
                query_exon_list.append(query_edge)
                attributes = query_edge[2]
                query_fasta_list.append('>query|{0:s}|{1:d}\n'.format(attributes['gene'], len(query_exon_list) - 1))
                query_fasta_list.append('{0:s}\n'.format(attributes['seq']))
        query_graph.graph['exon_list'] = query_exon_list

        subject_protein_fasta = proj_temp_dir + os.sep + subject_pickle_prefix + '.stage2.fasta'
        open(subject_protein_fasta, 'w').writelines(subject_fasta_list)
        query_protein_fasta = proj_temp_dir + os.sep + query_pickle_prefix + '.stage2.fasta'
        open(query_protein_fasta, 'w').writelines(query_fasta_list)

        tblastx = blast_wrapper.TBlastX(query=query_protein_fasta, subject=subject_protein_fasta,tmp_dir=proj_temp_dir)
        tblastx.run_tblastx()

        tblastx_output = [line.split('\t') for line in open(tblastx.out_file, 'r')]

        sorted_tblastx_output = sorted(tblastx_output, key=lambda blast_line: float(blast_line[10]))

        matched_exons = []
        new_query_graph = nx.Graph()
        new_subject_graph = nx.Graph()
        for query_edge in query_graph.edges_iter(data=True):

            if query_edge[2]['type'] != 'exon':
                new_query_graph.add_edge(query_edge[0], query_edge[1], query_edge[2])
        for subject_edge in subject_graph.edges_iter(data=True):

            if subject_edge[2]['type'] != 'exon':
                new_subject_graph.add_edge(subject_edge[0], subject_edge[1], subject_edge[2])

        for line in sorted_tblastx_output:
            if float(line[10]) < .01:

                query_edge_num = int(line[0].split('|')[2])
                subject_edge_num = int(line[1].split('|')[2])
                query_edge = query_graph.graph['exon_list'][query_edge_num]
                subject_edge = subject_graph.graph['exon_list'][subject_edge_num]

                if (query_edge_num, subject_edge_num) not in matched_exons:

                    query_edge[2]['depth'] = float(line[2])
                    subject_edge[2]['depth'] = float(line[2])
                    matched_exons.append((query_edge_num, subject_edge_num))
                    new_query_graph.add_edge(query_edge[0], query_edge[1], query_edge[2])
                    new_subject_graph.add_edge(subject_edge[0], subject_edge[1], subject_edge[2])

        # query_act_list = create_act_file_from_nx(new_query_graph, 'query')
        # subject_act_list = create_act_file_from_nx(new_subject_graph, 'subject')
        # query_act_file = proj_temp_dir + os.sep + 'query.tmp.act'
        # subject_act_file = proj_temp_dir + os.sep + 'subject.tmp.act'
        # open(query_act_file, 'w').writelines(query_act_list)
        # open(subject_act_file, 'w').writelines(subject_act_list)
        # tmp_query_act = actgraph.actFile(query_act_file)
        # tmp_query_act.Toimage('query', proj_temp_dir)
        # tmp_subject_act = actgraph.actFile(subject_act_file)
        # tmp_subject_act.Toimage('subject', proj_temp_dir)




if __name__ == '__main__':

    get_gene_stats()
