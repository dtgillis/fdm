__author__ = 'dtgillis'

import getopt
import sys
import os
from graph_utils import graph_constructor
from large_blast import blast_wrapper
import cPickle

def print_help():
    print 'gene_compare_stage_1.py -a <act_file> -g <gene_list> -o <out_dir>'
    print 'Use this to sort out the intial graphs and save them and make blastn files'

def get_gene_stats():
    query_act_file = None
    subject_act_file = None
    query_ref = None
    subject_ref = None
    proj_temp_dir = None
    blast_jobs = 1
    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "", ["query_act_file=", "subject_act_file=",
                                                      "query_ref=", "subject_ref=",
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
    if not os.path.exists(proj_temp_dir):
        try:
            os.makedirs(proj_temp_dir)
        except OSError:
            print "Error making dir {0:s}".format(proj_temp_dir)

    # get subject gene dict made
    subject_graph = graph_constructor.GenomeDiGraphDict(subject_act_file, subject_ref)
    subject_graph.create_genome_wide_graphs()
    # get connected componenets
    subject_graph.make_paths_to_search()
    # dump out the graphs objects
    act_file_name = subject_act_file.split(os.sep)[-1]
    cPickle.dump(subject_graph.path_dict, open(proj_temp_dir + os.sep + act_file_name + '.path_dict.pickle', 'wb'))

    # query graphs
    query_graph = graph_constructor.GenomeDiGraphDict(query_act_file, query_ref)
    query_graph.create_genome_wide_graphs()
    query_graph.make_paths_to_search()
    act_file_name_2 = query_act_file.split(os.sep)[-1]
    cPickle.dump(query_graph.path_dict, open(proj_temp_dir + os.sep + act_file_name_2 + '.path_dict.pickle', 'wb'))

    # create fast files for each path in subject / query graphs
    out_put = []
    subject_fasta_file_name = proj_temp_dir + os.sep + act_file_name + '.fasta'
    query_fasta_file_name = proj_temp_dir + os.sep + act_file_name_2 + '.fasta'
    subject_fasta = subject_graph.get_protein_fasta_from_paths()
    open(proj_temp_dir + os.sep + act_file_name + '.fasta', 'w').writelines(subject_fasta)
    query_fasta = query_graph.get_protein_fasta_from_paths()
    open(proj_temp_dir + os.sep + act_file_name_2 + '.fasta', 'w').writelines(query_fasta)
    tblastx = blast_wrapper.TBlastXBig(
        query_fasta_file_name, subject_fasta_file_name, blast_jobs, work_dir=proj_temp_dir)
    tblastx.break_down_query()
    tblastx.create_subject_db()
    tblastx.write_blast_submit_script()



    print 'Finished creating paths and nucleotide level blast files for files subject {0:s} and ' \
          'query {1:s}'.format(subject_act_file, query_act_file)
    exit()

if __name__ == '__main__':

    get_gene_stats()
