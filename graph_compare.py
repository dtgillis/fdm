__author__ = 'dtgillis'

import getopt
import sys
import os
from graph_utils import graph_constructor
import networkx as nx
import textwrap
import networkx as nx
import numpy as np
import matplotlib as plt
import Bio.pairwise2 as pairwise
from Bio.SubsMat.MatrixInfo import pam120
from Bio.Data import CodonTable
from Bio.Seq import translate
from large_blast import blast_wrapper

def print_help():
    print 'plot_act_gen.py -a <act_file> -g <gene_list> -o <out_dir>'


def get_gene_stats():
    mouse_act_file_list = None
    human_act_file_list = None
    gene_list = None
    out_dir = None
    mouse_ref = None
    human_ref = None

    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "m:c:g:o:", ["mouse_act_file_list=", "human_act_file_list",
                                                      "gene_homolog_list=", "out_dir=", "mouse_ref=", "human_ref="])
    except getopt.GetoptError:
        print 'plot_act_gen.py -m <mouse_act_file_list> -h <human_act_file_list> -g <gene_homolog_list> -o <out_dir>'
        sys.exit(2)

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            print_help()
        elif opt in ("-m", "--mouse_act_file_list"):
            mouse_act_file_list = arg.split(',')
        elif opt in ("-c", "--human_act_file_list"):
            human_act_file_list = arg.split(',')
        elif opt in ("-o", "--out_dir"):
            out_dir = arg
        elif opt in ("-g", "--gene_homolog_list"):
            gene_list = arg
        elif opt in ("--mouse_ref_dir"):
            mouse_ref = arg
        elif opt in ("--human_ref_dir"):
            human_ref = arg
        else:
            assert False, "invalid option"

    assert mouse_act_file_list is not None, "Specify mouse act file -m <mouse_act_file>"
    assert out_dir is not None, "Specify output directory -o <out_dir>"
    assert gene_list is not None, "Specify gene homolog list -g <gene_homolog_list>"
    assert mouse_ref is not None, "Specify mouse reference directory"
    assert human_ref is not None, "Specify human reference directory"

    # make sure the files all exist
    for act_file in mouse_act_file_list:
        if not os.path.isfile(act_file):
            print "{0:s} file does not exist".format(act_file)
            sys.exit(145)
    for act_file in human_act_file_list:
        if not os.path.isfile(act_file):
            print "{0:s} file does not exist".format(act_file)
            sys.exit(145)

    if not os.path.isfile(gene_list):
        print "{0:s} file does not exist".format(gene_list)
        sys.exit(1000)

    # make sure the path to outdir exists
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except OSError:
            print "Error making dir {0:s}".format(out_dir)

    # now check gene list is not an empty list

    assert len(mouse_act_file_list) > 0, "mouse act list has size 0, you should add some mice act files?"
    assert len(human_act_file_list) > 0, "human act list has size 0, you should add some human act files?"

    human_graphs = []
    human_gene_dict = dict()
    human_gene_dict_exons = dict()
    graph_counter = 0
    create_fasta = False
    for act_file in human_act_file_list:
        human_graph = graph_constructor.MultiGraph(act_file, human_ref)
        human_graph.create_genome_wide_graph()
        human_graphs.append(human_graph)
        sub_graphs = list(nx.connected_component_subgraphs(human_graph.graph, copy=True))

        for sub_graph in sub_graphs:
            key = ""
            for edge in sub_graph.edges_iter(sub_graph, data=True):
                attributes = edge[2]
                type = attributes['type']
                if type == "exon":
                    key += attributes['seq']

            if len(key) > 0:
                human_gene_dict[graph_counter] = sub_graph
                human_gene_dict_exons[graph_counter] = key
                graph_counter += 1

    #create fasta file
    if create_fasta:
        out_put = ""
        for key in human_gene_dict.keys():
            for edge in human_gene_dict[key].edges_iter(data=True):
                attributes = edge[2]
                if 'gene' in attributes:
                    gene = attributes['gene']
                    break


            out_put+='>human|{0:s}|{1:d}\n'.format(gene, key)
            exon_sequence = human_gene_dict_exons[key]
            formated_key = textwrap.wrap(exon_sequence, 80)
            for element in formated_key:
                out_put += '{0:s}\n'.format(element)
        fasta_file_loc = '/home/dtgillis/fall_rotation/tmp/human.fasta'
        fasta_file = open(fasta_file_loc, 'w')
        fasta_file.write(out_put)
        fasta_file.close()

    mouse_graphs = []
    mouse_gene_dict = dict()
    mouse_gene_dict_exons = dict()
    for act_file in mouse_act_file_list:
        mouse_graph = graph_constructor.MultiGraph(act_file, mouse_ref)
        mouse_graph.create_genome_wide_graph()
        mouse_graphs.append(mouse_graph)
        sub_graphs = list(nx.connected_component_subgraphs(mouse_graph.graph, copy=True))
        out_put = ""
        key_num = 0
        for sub_graph in sub_graphs:
            key = ""
            for edge in sub_graph.edges_iter(sub_graph, data=True):
                attributes = edge[2]
                type = attributes['type']
                if type == "exon":
                    key += attributes['seq']
            if len(key) > 0:
                mouse_gene_dict[key_num] = sub_graph
                mouse_gene_dict_exons[key_num] = key
                key_num += 1
        if create_fasta:
            for key in mouse_gene_dict.keys():
                sub_graph = mouse_gene_dict[key]
                exon_count = 1
                for edge in sub_graph.edges_iter(data=True):
                    attributes = edge[2]
                    edge_type = attributes['type']
                    if edge_type == "exon":
                        if len(attributes['seq']) > 0:
                            out_put += '>mouse|{0:s}|graph-{1:d}|exon-{2:d}\n'.format(attributes['gene'], key, exon_count)
                            out_put += '{0:s}\n'.format(attributes['seq'].lower())
                            exon_count += 1
            fasta_file = open('/home/dtgillis/fall_rotation/tmp/mouse.fasta', 'w')
            fasta_file.write(out_put)
            fasta_file.close()

    # run blast
    blast = blast_wrapper.blast_n(query_file='/home/dtgillis/fall_rotation/tmp/human.fasta',
                                  subject_file='/home/dtgillis/fall_rotation/tmp/mouse.fasta',
                                  num_jobs=1)
    blast.break_down_subject()

    #os.system(str(cmd))

    # read and parse blast output
    #blast_out = [line.strip(os.linesep).split() for line in open('/home/dtgillis/fall_rotation/tmp/blas','r')]

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
