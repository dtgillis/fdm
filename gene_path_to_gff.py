__author__ = 'dtgillis'

import getopt
import sys
import os
import cPickle
import math
import matplotlib.pyplot as plt
import common
from matplotlib.patches import Arc

import networkx as nx

from large_blast import blast_wrapper


def print_help():
    print 'graph_compare_stage_2'


def to_image(ts, gene, imagedir, genomerange=[0, 10000000000], highlightnodelist=[], ext='pdf', inwgs=[]):


    #image rendering constants
    hconst=30;gwidth=100.0;exonplotdistbasis=10.0;msize=10;hoffset=10;woffset=50

    wtgraphstruct=inwgs
    exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist=wtgraphstruct

    if genomerange!=[0,10000000000]:
        glabel='%s:%s:%d-%d'%(ts,gene,genomerange[0],genomerange[1])
    else:
        glabel='%s:%s'%(ts,gene)

    imgsz=len(exonlist)+len(intronlist)
    if imgsz>100:
        message='Image too large for %s. Printing first few exons'%glabel
        common.printstatus(message,'W',common.func_name(),1)
        genomerange=[exonlist[0:30][0][0]-1,exonlist[0:30][-1][1]+1]
        glabel='%s:%s:%d-%d'%(ts,gene,genomerange[0],genomerange[1])

    allpoints=[]
    for splice in splicelist:
        if splice[0]<genomerange[0] or splice[1]>genomerange[1]:
            continue
        if splice[0] not in allpoints:
            allpoints.append(splice[0])
        if splice[1] not in allpoints:
            allpoints.append(splice[1])
    for node in startnodelist:
        if node<genomerange[0] or node>genomerange[1]:
            continue
        if node not in allpoints:
            allpoints.append(node)
    for node in endnodelist:
        if node<genomerange[0] or node>genomerange[1]:
            continue
        if node not in allpoints:
            allpoints.append(node)
    for exon in exonlist:
        if exon[0]<genomerange[0] or exon[1]>genomerange[1]:
            continue
        if (exon[0] not in allpoints) and (exon[0]-1 not in allpoints) and (exon[0]+1 not in allpoints):
            allpoints.append(exon[0])
        if (exon[1] not in allpoints) and (exon[1]-1 not in allpoints) and (exon[1]+1 not in allpoints):
            allpoints.append(exon[1])
    for exon in intronlist:
        if exon[0]<genomerange[0] or exon[1]>genomerange[1]:
            continue
        if (exon[0] not in allpoints) and (exon[0]-1 not in allpoints) and (exon[0]+1 not in allpoints):
            allpoints.append(exon[0])
        if (exon[1] not in allpoints) and (exon[1]-1 not in allpoints) and (exon[1]+1 not in allpoints):
            allpoints.append(exon[1])


    allpoints.sort()
    #print exonlist
    #print allpoints
    pointdistlist=[(allpoints[i]-allpoints[i-1]) for i in range(1,len(allpoints))]
    mindist=max(exonplotdistbasis,min(pointdistlist))
    plotdistlist=[math.log(max(x,mindist),mindist) for x in pointdistlist]
    #equidistant
    #plotdistlist=[1 for x in pointdistlist]
    #print sum(plotdistlist)
    runtotal=0
    allpointsy=[runtotal]
    for point in plotdistlist:
        runtotal+=point
        allpointsy.append(runtotal*hconst)

    #print allpointsy
    pointplotdict=dict(zip(allpoints,allpointsy))

    gheight=sum(plotdistlist)*hconst
    #gwidth=max(gheight/10,100.0)
    ax = plt.figure(figsize=(4,4*gheight/gwidth)).add_subplot(111)
    ax.axis('off')


    #plot exons
    for i in range(len(exonlist)):
        exon=exonlist[i]
        if exon[0]<genomerange[0] or exon[1]>genomerange[1]:
            continue
        if exon[0] in pointplotdict:
            y1=pointplotdict[exon[0]]
        elif exon[0]-1 in pointplotdict:
            y1=pointplotdict[exon[0]-1]
        elif exon[0]+1 in pointplotdict:
            y1=pointplotdict[exon[0]+1]
        if exon[1] in pointplotdict:
            y2=pointplotdict[exon[1]]
        elif exon[1]-1 in pointplotdict:
            y2=pointplotdict[exon[1]-1]
        elif exon[1]+1 in pointplotdict:
            y2=pointplotdict[exon[1]+1]
        ax.plot([woffset,woffset],[y1+hoffset,y2+hoffset],'b-',markersize=msize)
        ax.text(woffset-1,(y1+y2)/2+hoffset,'%5.1f'%exonwtlist[i],horizontalalignment='right',verticalalignment='center',rotation=270)

    for i in range(len(intronlist)):
        intron=intronlist[i]
        if intron[0]<genomerange[0] or intron[1]>genomerange[1]:
            continue
        if intron[0] in pointplotdict:
            y1=pointplotdict[intron[0]]
        elif intron[0]-1 in pointplotdict:
            y1=pointplotdict[intron[0]-1]
        elif intron[0]+1 in pointplotdict:
            y1=pointplotdict[intron[0]+1]
        if intron[1] in pointplotdict:
            y2=pointplotdict[intron[1]]
        elif intron[1]-1 in pointplotdict:
            y2=pointplotdict[intron[1]-1]
        elif intron[1]+1 in pointplotdict:
            y2=pointplotdict[intron[1]+1]
        ax.plot([woffset,woffset],[y1+hoffset,y2+hoffset],'c-',markersize=msize)
        ax.text(woffset-1,(y1+y2)/2+hoffset,'%5.1f'%intronwtlist[i],horizontalalignment='right',verticalalignment='center',rotation=270)

    maxawidth=0
    for i in range(len(splicelist)):
        splice=splicelist[i]
        if splice[0]<genomerange[0] or splice[1]>genomerange[1]:
            continue
        if splicewtlist[i]>0:
            y1=pointplotdict[splice[0]]
            y2=pointplotdict[splice[1]]
            center=[woffset,(y1+y2)/2+hoffset]
            awidth=abs(gwidth/gheight*(y2-y1))
            if awidth>maxawidth:
                maxawidth=awidth

    for i in range(len(splicelist)):
        splice=splicelist[i]
        if splice[0]<genomerange[0] or splice[1]>genomerange[1]:
            continue
        if splicewtlist[i]>0:
            y1=pointplotdict[splice[0]]
            y2=pointplotdict[splice[1]]
            center=[woffset,(y1+y2)/2+hoffset]
            awidth=gwidth*2/gheight*(y2-y1)*(gwidth/maxawidth)
            arcs=[Arc(xy=center, width=awidth, height=y2-y1, angle=0, theta1=270, theta2=90,lw=2,color='green',ls='dashed')] # Arc
            #arcs=[Arc(xy=center, width=awidth, height=y2-y1, angle=0, theta1=270, theta2=90,lw=2,color='green')] # Arc
            ax.add_artist(arcs[0])
            ax.text(woffset+awidth/2+1,(y1+y2)/2+hoffset,splicewtlist[i],horizontalalignment='left',verticalalignment='center', rotation=270)



    #plot nodes
    #Todo text for points
    for i in range(len(allpoints)):
        node=allpoints[i]
        if node in startnodelist:
            ax.plot(woffset,allpointsy[i]+hoffset,marker='s',markerfacecolor='k',fillstyle='bottom',markersize=msize)
        elif node in endnodelist:
            ax.plot(woffset,allpointsy[i]+hoffset,marker='s',markerfacecolor='k',fillstyle='top',markersize=msize)
        elif node in novelnodelist:
            ax.plot(woffset,allpointsy[i]+hoffset,marker='o',markerfacecolor='white',markersize=msize)
        else:
            ax.plot(woffset,allpointsy[i]+hoffset,marker='o',markerfacecolor='k',fillstyle='full',markersize=msize)
        if node in highlightnodelist:
            ax.text(woffset-msize,allpointsy[i]+hoffset,node,horizontalalignment='right',verticalalignment='center',color='red')
        else:
            ax.text(woffset-msize,allpointsy[i]+hoffset,node,horizontalalignment='right',verticalalignment='center')


    ax.text(0,max(allpointsy)+hoffset+4*msize,glabel,horizontalalignment='left',verticalalignment='center')
    ax.set_ylim(0, max(allpointsy)+hoffset+100)
    #ax.set_ylim(100, 200)
    ax.set_xlim(0, gwidth+woffset)

    #plt.show()
    ax.set_aspect('equal')
    try:
        plt.savefig('%s/%s.%s'%(imagedir,glabel.replace(':','__').replace('/','_'),ext),bbox_inches='tight', pad_inches=0.1)
    except:
        message='Image too large for %s'%glabel
        common.printstatus(message,'W',common.func_name(),1)

def create_gff_file_from_nx(graph):

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

def get_plot_input(gene_path,):
    query_nodes = gene_path.nodes()
    node_list = []
    for node in query_nodes:
        node_list.append(int(node.split(":")[-1]))

    sorted_node_list = sorted(node_list)
    exon_list = []
    exon_wt_list = []
    splice_list = []
    splice_wt_list = []
    ret_int_list = []
    ret_int_wt_list = []
    for edge in gene_path.edges_iter(data=True):

        attributes = edge[2]
        if attributes['type'] == 'exon':
            exon_list.append([int(edge[0].split(":")[-1]), int(edge[1].split(":")[-1])])
            exon_wt_list.append(float(attributes['depth']))
        elif attributes['type'] == 'splice':
            splice_list.append([int(edge[0].split(":")[-1]), int(edge[1].split(":")[-1])])
            splice_wt_list.append(float(attributes['depth']))

    end_node_list = [sorted_node_list[-1]]
    start_node_list = [sorted_node_list[0]]
    plot_list = (exon_list, [], splice_list, start_node_list, end_node_list, [], exon_wt_list, [], splice_wt_list)

    return plot_list


def write_out_gff(dir, gene_path, sample_name, score):

    out_list = []
    plot_list = get_plot_input(gene_path)
    exon_list = plot_list[0]
    for edge in gene_path.edges_iter(data=True):
        chrm = edge[0].split(":")[0]
        break
    out_list.append('track name={0:s} description="{1:s} {2:f}"\n'.format(gene_path.graph['gene'], sample_name, score))
    for exon in exon_list:
        out_list.append('{0:s}\t{1:s}\t{2:s}\t{3:d}\t{4:d}\t.\t.\t.\t.\n'.format(
                chrm, 'rna-seq', 'exon', exon[0], exon[1]))
        # elif attributes['type'] == 'splice':
        #     out_list.append('{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\t.\t.\t.\t.\n'.format(
        #         edge[0].split(":")[0], 'rna-seq', 'splice', edge[0].split(":")[-1], edge[1].split(":")[-1]))
    open(dir + os.sep + sample_name + '_' + gene_path.graph['gene'] + '.gff', 'w').writelines(out_list)


def get_gene_stats():
    query_pickle_prefix = None
    subject_pickle_prefix = None
    proj_temp_dir = None
    gene_path_file = None
    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "", ["proj_temp_dir=", "query_prefix=", "subject_prefix=", "gene_path_list="])
    except getopt.GetoptError:
        print 'gene_path_to_gff.py --proj_temp_dir=temp_dir ' \
              '--query_prefix=prefix_for_query --subject_prefix=prefix_for_subject --gene_path_list'
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
        elif opt in ("--gene_path_list"):
            gene_path_file = arg
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
    for line in open(gene_path_file, 'r'):
        line = line.strip(os.linesep).split(',')
        query_graph = query_gene_dict[int(line[1])]
        subject_graph = subject_gene_dict[int(line[3])]
        score = float(line[-1])
        query_plot_data = get_plot_input(query_graph)
        subject_plot_data = get_plot_input(subject_graph)
        write_out_gff(proj_temp_dir, query_graph, query_pickle_prefix[:-4], score)
        write_out_gff(proj_temp_dir, subject_graph, subject_pickle_prefix[:-4], score)

        to_image('{0:s}-{1:s}.{2:4.1f}'.format(
            query_pickle_prefix[:-4], subject_pickle_prefix[:-4], score),
                 query_graph.graph['gene'], proj_temp_dir, inwgs=query_plot_data)
        to_image('{0:s}-{1:s}.{2:4.1f}'.format(
            subject_pickle_prefix[:-4], query_pickle_prefix[:-4], score),
                 query_graph.graph['gene'], proj_temp_dir, inwgs=subject_plot_data)







if __name__ == '__main__':

    get_gene_stats()
