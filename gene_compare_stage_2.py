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

    matches_dict = dict()
    for line in blast_results_sorted:
        if float(line[10]) < .001:
            query_graph = int(line[0].split('|')[2].split('-')[-1])
            subject_graph = int(line[1].split('|')[-1])
            if query_graph not in matches_dict:
                matches_dict[query_graph] = set()
            matches_dict[query_graph].add(subject_graph)

    #print matches_dict
    digraph_matches = nx.DiGraph()
    print (len(matches_dict.keys()))
    count = 0
    for query_graph_num in matches_dict:
        count += 1
        if count % 100 == 0:
            sys.stdout.write('.')
        if count %1000 == 0:
            sys.stdout.write(os.linesep)

        subject_graph_nums = matches_dict[int(query_graph_num)]
        query_graph = query_gene_dict[int(query_graph_num)]
        #one subject fasta per potential matches.
        query_fasta_list = []
        query_graph_exons = []
        for query_edge in query_graph.edges_iter(data=True):

            if query_edge[2]['type'] == 'exon':
                query_graph_exons.append(query_edge)
                attributes = query_edge[2]
                query_fasta_list.append('>subject|graph-{0:d}|{1:s}|{2:d}\n'.format(
                    query_graph_num, attributes['gene'], len(query_graph_exons)-1))
                query_fasta_list.append('{0:s}\n'.format(attributes['seq']))

        query_graph.graph['exon_list'] = query_graph_exons
        subject_fasta_list = []

        for subject_graph_num in subject_graph_nums:
            #process each query graph
            subject_graph = subject_gene_dict[subject_graph_num]
            subject_exon_list = []
            for subject_edge in subject_graph.edges_iter(data=True):
                 if subject_edge[2]['type'] == 'exon':
                    subject_exon_list.append(subject_edge)
                    attributes = subject_edge[2]
                    subject_fasta_list.append('>query|graph-{0:d}|{1:s}|{2:d}\n'.format(
                        subject_graph_num, attributes['gene'], len(subject_exon_list) - 1))
                    subject_fasta_list.append('{0:s}\n'.format(attributes['seq']))
            subject_graph.graph['exon_list'] = subject_exon_list

        subject_protein_fasta = proj_temp_dir + os.sep + subject_pickle_prefix + '.stage2.fasta'
        open(subject_protein_fasta, 'w').writelines(subject_fasta_list)
        query_protein_fasta = proj_temp_dir + os.sep + query_pickle_prefix + '.stage2.fasta'
        open(query_protein_fasta, 'w').writelines(query_fasta_list)

        tblastx = blast_wrapper.TBlastX(query=query_protein_fasta, subject=subject_protein_fasta, tmp_dir=proj_temp_dir)
        tblastx.run_tblastx()

        tblastx_output = [line.split('\t') for line in open(tblastx.out_file, 'r')]

        sorted_tblastx_output = sorted(tblastx_output, key=lambda blast_line: float(blast_line[10]))

        matched_exons = []

        weights = dict()
        for line in sorted_tblastx_output:

            if float(line[10]) < .01:

                tmp_query_graph_num = int(line[0].split('|')[1].split('-')[-1])
                tmp_subject_graph_num = int(line[1].split('|')[1].split('-')[-1])
                weight = float(line[-1].strip(os.linesep))
                if not digraph_matches.has_edge(tmp_query_graph_num, tmp_subject_graph_num):
                    digraph_matches.add_edge(tmp_query_graph_num, tmp_subject_graph_num, weight=weight)
                else:
                    digraph_matches[tmp_query_graph_num][tmp_subject_graph_num]['weight'] += weight


    overall_mapping = []
    for query_graph_num in matches_dict:
        tmp_edge = None
        for edge in digraph_matches.edges([query_graph_num], data=True):
            if edge[0] == query_graph_num:
                if tmp_edge == None:
                    tmp_edge = edge
                else:
                    if tmp_edge[2]['weight'] < edge[2]['weight']:
                        tmp_edge = edge

        if tmp_edge is not None:
            overall_mapping.append([query_graph_num, tmp_edge[1], tmp_edge[2]['weight']])
        else:
            overall_mapping.append([query_graph_num, None, 0.0])

    matching_out = []
    for mapping in overall_mapping:

        query_gene = query_gene_dict[mapping[0]].graph['exon_list'][0][2]['gene']
        if mapping[1] is not None:
            subject_graph_gene = subject_gene_dict[mapping[1]].graph['exon_list'][0][2]['gene']
        else:
            subject_graph_gene = 'None found'
            mapping[1] = 'None'

        matching_out.append('{0:s},{1:s},{2:s},{3:s},{4:f}\n'.format(
            query_gene, str(mapping[0]), subject_graph_gene, str(mapping[1]), mapping[2]
        ))

    open(proj_temp_dir + os.sep + 'matching.out.dat', 'w').writelines(matching_out)



if __name__ == '__main__':

    get_gene_stats()
