__author__ = 'dtgillis'

import getopt
import sys
import os
import cPickle

import networkx as nx

from large_blast import blast_wrapper


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
    query_pickle_gene_dict = proj_temp_dir + os.sep + query_pickle_prefix + '.path_dict.pickle'
    assert os.path.isfile(query_pickle_gene_dict), "Query Gene Dict non existent"
    subject_pickle_gene_dict = proj_temp_dir + os.sep + subject_pickle_prefix + '.path_dict.pickle'

    # load all pickle files
    query_gene_dict = cPickle.load(open(query_pickle_gene_dict, 'rb'))
    subject_gene_dict = cPickle.load(open(subject_pickle_gene_dict, 'rb'))

    # read in blast out files
    #blast_results = [line.split('\t') for line in open(blast_out_file, 'r').readlines()]
    #blast_results_sorted = sorted(blast_results, key=lambda blast_line: float(blast_line[10]))
    alpha = 1
    beta = 1
    matches_dict = dict()
    log_fit_file = proj_temp_dir + os.sep + query_pickle_prefix + '.log_reg.dat'
    log_fit_out = []
    for line in open(blast_out_file, 'r'):
        line = line.strip(os.sep).split('\t')
        if float(line[10]) < .01:
            query_graph = int(line[0].split('|')[1].split('-')[-1])
            subject_graph = int(line[1].split('|')[1].split('-')[-1])
            score = float(line[2])*float(line[3])/100.0 - alpha * float(line[4]) - beta * float(line[5])
            log_fit_out.append('{0:s},{1:s},{2:s},{3:s},{4:s},{5:s}\n'.format(
                line[0].split('|')[2], line[1].split('|')[2], line[2], line[3], line[4], line[5]))
            if query_graph not in matches_dict:
                matches_dict[query_graph] = dict()
                matches_dict[query_graph][subject_graph] = score
            else:
                if subject_graph not in matches_dict[query_graph]:
                    matches_dict[query_graph][subject_graph] = score
                elif score > matches_dict[query_graph][subject_graph]:
                    matches_dict[query_graph][subject_graph] = score


    #print matches_dict
    open(log_fit_file, 'w').writelines(log_fit_out)
    print (len(matches_dict.keys()))
    count = 0
    matches_out = []
    for query_graph_num in matches_dict:
        count += 1
        if count % 100 == 0:
            sys.stdout.write('.')
        if count %1000 == 0:
            sys.stdout.write(os.linesep)
        query_path = query_gene_dict[query_graph_num]
        match_list = matches_dict[query_graph_num].items()
        sorted_match_list = sorted(match_list, key=lambda match: match[1])
        subject_path = subject_gene_dict[sorted_match_list[-1][0]]
        query_gene_name = query_path.graph['gene']
        subject_gene_name = subject_path.graph['gene']
        matches_out.append('{0:s},{1:d},{2:s},{3:d},{4:f}\n'.format(
            query_gene_name, query_graph_num, subject_gene_name, sorted_match_list[-1][0], sorted_match_list[-1][1]))

    sorted_matches_out = sorted(matches_out, key=lambda match: match[-1].strip(os.linesep), reverse=True)
    open(proj_temp_dir + os.sep + 'matching.out.a{0:1.1f}-b{1:1.1f}.dat'.format(alpha, beta), 'w').writelines(sorted_matches_out)



if __name__ == '__main__':

    get_gene_stats()
