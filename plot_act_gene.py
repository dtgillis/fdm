__author__ = 'dtgillis'

import actgraph
import getopt
import sys
import os


def print_help():
    print 'plot_act_gen.py -a <act_file> -g <gene_list> -o <out_dir>'


def plot_genes():
    act_file = None
    gene_list = None
    out_dir = None

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "a:gl:o:", ["act_file=", "gene_list=", "out_dir="])
    except getopt.GetoptError:
        print 'plot_act_gen.py -a <act_file> -g <gene_list> -o <out_dir>'
        sys.exit(2)

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            print_help()
        elif opt in ("-a", "--act_file"):
            act_file = arg
        elif opt in ("-o", "--out_dir"):
            out_dir = arg
        elif opt in ("-gl", "--gene_list"):
            gene_list = arg
        else:
            assert False, "invalid option"

    assert act_file is not None, "Specify act file -a <act_file>"
    assert out_dir is not None, "Specify output directory -o <out_dir>"
    assert gene_list is not None, "Specify gene list -g <gene_list>"

    # make sure the files all exist
    if not os.path.isfile(act_file):
        print "{0:s} file does not exist".format(act_file)

    # make sure the path to outdir exists
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except OSError:
            print "Error making dir {0:s}".format(out_dir)

    # now check gene list is not an empty list

    assert len(gene_list) > 0, "Gene list has size 0, you should print some genes?"

    # create the act object now
    act_graph = actgraph.actFile(act_file)
    for gene in gene_list:
        act_graph.Toimage(gene, out_dir)

if __name__ == '__main__':
    plot_genes()
