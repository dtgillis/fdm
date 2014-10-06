__author__ = 'dtgillis'
import networkx as nx
import pysam
import os


class MultiGraph:

    def __init__(self, act_file, reference_dir):

        self.act_file = act_file
        self.reference_dir = reference_dir
        self.graph = nx.Graph()

    def create_genome_wide_graph(self):

        for line in open(self.act_file,'r'):

            fields = line.split()

            if fields[2] != "node" and fields[2] != "retint":

                chrm = fields[0]
                fasta_file = pysam.Fastafile(self.reference_dir + os.sep + chrm + ".fa")
                edge_type = fields[2]
                weight = float(fields[5])
                start = int(fields[3])
                end = int(fields[4])
                gene = fields[-1]
                if edge_type == "exon":
                    sequence = fasta_file.fetch(region="{0:s}:{1:d}-{2:d}".format(chrm, start, end))
                    self.graph.add_edge(
                        "{0:s}:{1:d}".format(chrm, start), "{0:s}:{1:d}".format(chrm, end),
                        chrm=chrm, type=edge_type, depth=weight, seq=sequence, gene=gene)
                else:
                    self.graph.add_edge(
                        "{0:s}:{1:d}".format(chrm, start), "{0:s}:{1:d}".format(chrm, end),
                        chrm=chrm, type=edge_type, depth=weight)




