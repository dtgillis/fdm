__author__ = 'dtgillis'

import getopt
import sys
import os
import gene_homolog.common_genes as genes



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

    stat_out_file = out_dir + os.sep + 'stat_file.dat'

    gene_sets = genes.GeneDictionary(gene_list)

    gene_sets.parse_gene_file()

    gene_compare = genes.GeneCompare(mouse_act_file_list, human_act_file_list, gene_sets, stat_out_file, mouse_ref, human_ref)

    gene_compare.load_act_graphs()
    #gene_compare.rebuild_index()
    gene_compare.cross_wise_simple_stats()




if __name__ == '__main__':

    get_gene_stats()
