__author__ = 'dtgillis'

import getopt
import sys
import os
from graph_utils import graph_constructor
import networkx as nx
import textwrap
import networkx as nx
from large_blast import blast_wrapper
import cPickle

def print_help():
    print 'gene_compare_stage_1.py -a <act_file> -g <gene_list> -o <out_dir>'
    print 'Use this to sort out the intial graphs and save them and make blastn files'

def get_gene_stats():
    query_act_file = None
    subject_act_file = None
    out_dir = None
    query_ref = None
    subject_ref = None
    proj_temp_dir = None
    blast_jobs = 1
    argv = sys.argv[1:]
    #TODO make more generic ie not mouse or human
    try:
        opts, args = getopt.getopt(argv, "", ["query_act_file=", "subject_act_file=",
                                                      "out_dir=", "query_ref=", "subject_ref=",
                                                      "proj_temp_dir=", "blast_jobs="])
    except getopt.GetoptError:
        print getopt.GetoptError.msg
        print 'gene_compare_stage_1.py --query_act_file=<mouse_act_file> --query_ref=<query_ref_files> ' \
              '--subject_act_file=<human_act_file> --subject_ref=<subject_ref_files> --proj_temp_dir=<tmp_space>'
        sys.exit(2)

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            print_help()
        elif opt in ("-q", "--query_act_file"):
            query_act_file = arg
        elif opt in ("-s", "--subject_act_file"):
            subject_act_file = arg
        elif opt in ("-o", "--out_dir"):
            out_dir = arg
        elif opt in ("--query_ref"):
            query_ref = arg
        elif opt in ("--subject_ref"):
            subject_ref = arg
        elif opt in ("--proj_temp_dir"):
            proj_temp_dir = arg
        elif opt in ("--blast_jobs"):
            blast_jobs = int(arg)
        else:
            assert False, "invalid option"

    assert query_act_file is not None, "Specify mouse act file -m <mouse_act_file>"
    assert out_dir is not None, "Specify output directory -o <out_dir>"
    assert query_ref is not None, "Specify query reference directory"
    assert subject_ref is not None, "Specify subject reference directory"

    # make sure the files all exist
    if not os.path.isfile(query_act_file):
        print "{0:s} file does not exist".format(query_act_file)
        sys.exit(145)
    if not os.path.isfile(subject_act_file):
        print "{0:s} file does not exist".format(subject_act_file)
        sys.exit(145)

    # make sure the path to outdir exists
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except OSError:
            print "Error making dir {0:s}".format(out_dir)

    # get subject gene dict made
    subject_gene_dict = dict()
    subject_gene_dict_exons = dict()
    graph_counter = 0
    subject_graph = graph_constructor.MultiGraph(subject_act_file, subject_ref)
    subject_graph.create_genome_wide_graph()
    # get connected componenets
    sub_graphs = list(nx.connected_component_subgraphs(subject_graph.graph, copy=True))
    for sub_graph in sub_graphs:
        key = ""
        for edge in sub_graph.edges_iter(sub_graph, data=True):
            attributes = edge[2]
            type = attributes['type']
            if type == "exon":
                key += attributes['seq']

        if len(key) > 0:
            subject_gene_dict[graph_counter] = sub_graph
            subject_gene_dict_exons[graph_counter] = key
            graph_counter += 1

    # dump out the graphs objects
    act_file_name = subject_act_file.split(os.sep)[-1]
    cPickle.dump(subject_gene_dict, open(proj_temp_dir + os.sep + act_file_name + '.gene_dict.pickle', 'wb'))
    cPickle.dump(subject_gene_dict_exons, open(proj_temp_dir + os.sep + act_file_name + '.gene_dict_exons.pickle', 'wb'))

    # create fasta file for db construction
    out_put = []
    for key in subject_gene_dict.keys():
        for edge in subject_gene_dict[key].edges_iter(data=True):
            attributes = edge[2]
            if 'gene' in attributes:
                gene = attributes['gene']
                break


        out_put.append('>subject|{0:s}|{1:d}\n'.format(gene, key))
        exon_sequence = subject_gene_dict_exons[key].lower()
        out_put.append('{0:s}\n'.format(exon_sequence))

    fasta_file_loc = proj_temp_dir + os.sep + act_file_name + '.fasta'
    fasta_file = open(fasta_file_loc, 'w')
    fasta_file.writelines(out_put)
    fasta_file.close()

    query_gene_dict = dict()
    query_gene_dict_exons = dict()
    query_graph = graph_constructor.MultiGraph(query_act_file, query_ref)
    query_graph.create_genome_wide_graph()
    # get connected components
    sub_graphs = list(nx.connected_component_subgraphs(query_graph.graph, copy=True))
    out_put = []
    key_num = 0
    for sub_graph in sub_graphs:
        key = ""
        for edge in sub_graph.edges_iter(sub_graph, data=True):
            attributes = edge[2]
            type = attributes['type']
            if type == "exon":
                key += attributes['seq']
        if len(key) > 0:
            query_gene_dict[key_num] = sub_graph
            query_gene_dict_exons[key_num] = key
            key_num += 1
    for key in query_gene_dict.keys():
        sub_graph = query_gene_dict[key]
        exon_count = 1
        for edge in sub_graph.edges_iter(data=True):
            attributes = edge[2]
            edge_type = attributes['type']
            if edge_type == "exon":
                if len(attributes['seq']) > 0:
                    out_put.append('>query|{0:s}|graph-{1:d}|exon-{2:d}\n'.format(attributes['gene'], key, exon_count))
                    out_put.append('{0:s}\n'.format(attributes['seq'].lower()))
                    exon_count += 1

    # pickle out the graph dicts
    act_file_name_2 = query_act_file.split(os.sep)[-1]
    cPickle.dump(query_gene_dict, open(proj_temp_dir + os.sep + act_file_name_2 + '.gene_dict.pickle', 'wb'))
    cPickle.dump(query_gene_dict_exons, open(proj_temp_dir + os.sep + act_file_name_2 + '.gene_dict_exons.pickle', 'wb'))

    fasta_file = open(proj_temp_dir + os.sep + act_file_name_2 + '.fasta', 'w')
    fasta_file.writelines(out_put)
    fasta_file.close()

    # create blast files in temp directory
    query_act_file_name = query_act_file.split(os.sep)[-1]
    subject_act_file_name = subject_act_file.split(os.sep)[-1]
    blast = blast_wrapper.BlastN(query_file=proj_temp_dir + os.sep + query_act_file_name + '.fasta',
                                 subject_file=proj_temp_dir + os.sep + subject_act_file_name + '.fasta',
                                 num_jobs=blast_jobs, work_dir=proj_temp_dir)

    blast.break_down_query()
    blast.create_subject_db()
    blast.write_blast_submit_script()

    print 'Finished creating subgraphs and nucleotide level blast files for files subject {0:s} and ' \
          'query {1:s}'.format(subject_act_file, query_act_file)
    exit()
    # #os.system(str(cmd))
    #
    # # read and parse blast output
    # #blast_out = [line.strip(os.linesep).split() for line in open('/home/dtgillis/fall_rotation/tmp/blas','r')]
    #
    # blast_out = sorted(blast_out, key=lambda blast_line: float(blast_line[10]))
    #
    # for line in blast_out:
    #     print line
    #     human_graph_key = int(line[0].split("|")[-1])
    #     mouse_graph_key = int(line[1].split("|")[2].strip("graph-"))
    #     human_graph = subject_gene_dict[human_graph_key]
    #     mouse_graph = query_gene_dict[mouse_graph_key]
    #
    #     #get exon counts for both species
    #     human_exon_list = []
    #     for human_edge in human_graph.edges_iter(data=True):
    #         human_attributes = human_edge[2]
    #         if human_attributes['type'] == "exon":
    #             human_exon_list.append(human_edge)
    #     mouse_exon_list = []
    #     for mouse_edge in mouse_graph.edges_iter(data=True):
    #         mouse_attributes = mouse_edge[2]
    #         if mouse_attributes['type'] == "exon":
    #             mouse_exon_list.append(mouse_edge)
    #
    #     score_mat = np.zeros((len(mouse_exon_list)*3, len(human_exon_list)))
    #     for j in range(len(human_exon_list)):
    #         human_attributes = human_exon_list[j][2]
    #         if human_attributes['type'] != "exon":
    #             continue
    #         human_codon = translate(human_attributes['seq']).replace('*', '')
    #
    #         for i in range(len(mouse_exon_list)):
    #             mouse_attributes = mouse_exon_list[i][2]
    #             if mouse_attributes['type'] != "exon":
    #                 continue
    #             if len(mouse_attributes['seq']) > 3:
    #                 mouse_seq = [mouse_attributes['seq'],
    #                                 mouse_attributes['seq'][1:], mouse_attributes['seq'][2:]]
    #             else:
    #                 mouse_seq = [mouse_attributes['seq']*3]
    #             mouse_codons = [translate(seq).replace('*', '') for seq in mouse_seq]
    #             k = 0
    #             for mouse_codon in mouse_codons:
    #                 if mouse_codon =="" or human_codon == "":
    #                     pam_score=0.0
    #                 else:
    #                     pam_score = pairwise.align.globaldx(human_codon, mouse_codon, pam120, score_only=True)
    #                     #print pam_score
    #                     #print pairwise.format_alignment(*pairwise.align.globaldx(human_codon, mouse_codon, pam120, one_alignment_only=True)[0])
    #
    #                 score_mat[i*3 + k, j] = pam_score
    #                 k += 1
    #
    #     matches = np.where(score_mat == np.max(score_mat, axis=0))
    #
    #     for i in range(len(human_exon_list)):
    #         print human_exon_list[matches[1][i]]
    #         mouse_exon = matches[0][i]
    #         mouse_exon_shift = mouse_exon % 3
    #         mouse_range = mouse_exon/3
    #         print mouse_exon_list[mouse_range][2]['seq'][mouse_exon_shift:]

if __name__ == '__main__':

    get_gene_stats()
