__author__ = 'dtgillis'
import networkx as nx
import pysam
import os


class GenomeDiGraphDict:

    def __init__(self, act_file, reference_dir):

        self.act_file = act_file
        self.reference_dir = reference_dir
        self.gene_struct_dict = dict()
        self.path_dict = dict()

    def create_genome_wide_graphs(self):

        act_graph_lines = [line.strip(os.linesep).split() for line in open(self.act_file, 'r')]

        for line in act_graph_lines:

            if line[2] != "node":
                if line[-1] not in self.gene_struct_dict:
                    self.gene_struct_dict[line[-1]] = nx.DiGraph()
                chrm = line[0]
                fasta_file = pysam.Fastafile(self.reference_dir + os.sep + chrm + ".fa")
                edge_type = line[2]
                weight = float(line[5])
                start = int(line[3])
                end = int(line[4])
                gene = line[-1]
                length = abs(end - start)
                if edge_type == "exon":
                    sequence = fasta_file.fetch(region="{0:s}:{1:d}-{2:d}".format(chrm, start, end))
                    self.gene_struct_dict[line[-1]].add_edge(
                        "{0:s}:{1:d}".format(chrm, start), "{0:s}:{1:d}".format(chrm, end),
                        chrm=chrm, type=edge_type, depth=weight, seq=sequence, gene=gene, length=length)
                    if 'gene' not in self.gene_struct_dict[line[-1]].graph:
                        self.gene_struct_dict[line[-1]].graph['gene'] = gene
                else:
                    self.gene_struct_dict[line[-1]].add_edge(
                        "{0:s}:{1:d}".format(chrm, start), "{0:s}:{1:d}".format(chrm, end),
                        chrm=chrm, type=edge_type, depth=weight, length=0)

        for graph in self.gene_struct_dict.values():
            max_value = 0.0
            for edge in graph.edges_iter(data=True):
                if edge[2]['type'] == 'exon':
                    max_value = max(max_value, edge[2]['length'])

            for edge in graph.edges_iter(data=True):
                if edge[2]['type'] == 'exon':
                    edge[2]['trans_weight'] = max_value - edge[2]['length']

    def make_paths_to_search(self):

        sub_graphs = self.gene_struct_dict.values()
        graph_counter = 0
        for sub_graph in sub_graphs:
            in_degrees = list(sub_graph.in_degree_iter())
            sorted_in_degree = sorted(in_degrees, key=lambda element: element[1])
            start_nodes = []
            for degree in sorted_in_degree:
                if degree[1] == 0:
                    start_nodes.append(degree[0])
                else:
                    break
            out_degrees = list(sub_graph.out_degree_iter())
            sorted_out_degree = sorted(out_degrees, key=lambda element: element[1])
            out_nodes = []
            for degree in sorted_out_degree:
                if degree[1] == 0:
                    out_nodes.append(degree[0])
                else:
                    break

            path_lengths = nx.all_pairs_dijkstra_path_length(sub_graph, weight='trans_length')
            max_distance = []
            for start_node in start_nodes:
                for end_node in out_nodes:
                    if start_node in path_lengths:
                        if end_node in path_lengths[start_node]:
                            max_distance.append((start_node, end_node, path_lengths[start_node][end_node]))

            sorted_path_length = sorted(max_distance, key=lambda key: key[-1])
            all_paths = nx.all_pairs_dijkstra_path(sub_graph,  weight='trans_length')
            paths_to_choose = []
            for path_length in sorted_path_length:
                paths_to_choose.append([all_paths[path_length[0]][path_length[1]]])

            for path in paths_to_choose:
                path_depth = 0
                path_length = 0
                for i in range(len(path[0])-1):
                    path_depth += sub_graph.get_edge_data(path[0][i], path[0][i+1])['depth']
                    path_length += sub_graph.get_edge_data(path[0][i], path[0][i+1])['length']
                path.append(path_length)
                path.append(path_depth)
                if float(path[1]) != 0:
                    path.append(path_depth/float(path[1]))
                else:
                    path.append(0.0)

            final_sort = sorted(paths_to_choose, key=lambda stuff: (stuff[1], stuff[-1]))
            graph_path = sub_graph.subgraph(final_sort[-1][0])
            self.path_dict[graph_counter] = graph_path
            graph_counter += 1

    def get_protein_fasta_from_paths(self):

        tmp_out = []
        name = self.act_file.split(os.sep)[-1]
        for graph_num in self.path_dict:
            graph = self.path_dict[graph_num]
            if 'gene' not in graph.graph:
                continue
            else:
                gene = graph.graph['gene']
            exon_count = 0
            protein = ''
            for edge in graph.edges_iter(data=True):
                if edge[2]['type'] == 'exon':
                    attributes = edge[2]
                    # tmp_out.append('>{0:s}|graph-{1:d}|{2:s}|{3:d}\n'.format(
                    #     name, graph_num, gene, exon_count))
                    # tmp_out.append('{0:s}\n'.format(attributes['seq']))
                    exon_count += 1
                    protein += attributes['seq']
            if len(protein) > 0:
                tmp_out.append('>{0:s}|graph-{1:d}|{2:s}|{3:d}\n'.format(
                    name, graph_num, gene, exon_count))
                tmp_out.append('{0:s}\n'.format(protein))

        return tmp_out





