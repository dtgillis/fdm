import common
import numpy as np
import math
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

def plotscatterwithhistogram(xlist,ylist,xlabel,ylabel,title,imagefilename,logscalexflg=0,logscaleyflg=0,markertup=('.',5)):
    fig = plt.figure(figsize=(16, 16))
    ax1 = fig.add_subplot(6,2,5)
    ax2 = fig.add_subplot(2,2,3)
    ax3 = fig.add_subplot(2,6,10)
    
    ax2.set_title(title)
    if logscalexflg:
        ax2.set_xlabel('log10 %s'%xlabel)
    else:
        ax2.set_xlabel(xlabel)    
    if logscaleyflg:
        ax2.set_ylabel('log10 %s'%ylabel)
    else:
        ax2.set_ylabel(ylabel)
    

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)

    if logscalexflg:
        xlist1=[max(x,1) for x in xlist]
        plotxlist=[math.log10(x) for x in xlist1]
    else:
        plotxlist=xlist[0:]           
    if logscaleyflg:
        ylist1=[max(y,1) for y in ylist]
        plotylist=[math.log10(y) for y in ylist1]
    else:
        plotylist=ylist[0:]

    ax1.hist(plotxlist,bins=25)
    ax3.hist(plotylist,bins=25,orientation='horizontal')
    ax2.plot(plotxlist,plotylist,color='k',marker=markertup[0],markersize=markertup[1],linestyle='None')
    #plt.subplot(6,2,5)
    plt.savefig(imagefilename)

def plot_histogram(data,num_bins,xlabel,ylabel,ptitle,imagefilename,log_flg=0):
    fig = plt.figure(figsize=(16, 16))
    if log_flg==0:
        x=data[0:]
    else:
        x=[math.log10(xi) for xi in data if xi>0.0]
        xlabel='log '+xlabel
    n, bins, patches = plt.hist(x, num_bins, facecolor='green', alpha=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(ptitle)
    plt.savefig(imagefilename)

def plotFEFclusters(FEFdatadict,clusterlist,medianlist,genechredge,imagefilename):
    tslist=FEFdatadict.keys()
    goodtslist=list(itertools.chain(*clusterlist))
    fig = plt.figure(figsize=(16, 16))
    missingXY=np.array([[FEFdatadict[ts][0],math.log10(max(FEFdatadict[ts][1],0.1))] for ts in FEFdatadict if ts not in goodtslist])
    clusterpoints=[]
    for cluster in clusterlist:
        clusterXY=[[FEFdatadict[ts][0],math.log10(max(FEFdatadict[ts][1],0.1))] for ts in FEFdatadict if ts in cluster]
        clusterpoints.append(np.array(clusterXY))
    #clusterpoints=np.array(clusterpoints)
    #print clusterpoints
    
    ax1 = fig.add_subplot(6,2,5)
    ax2 = fig.add_subplot(2,2,3)
    ax3 = fig.add_subplot(2,6,10)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)        
    pxrange=[0.0,1.0]
    pyrange=[-1.0,4.0]
    
    if len(clusterpoints)==3:
        ax2.plot(clusterpoints[0][:,0],clusterpoints[0][:,1],'ob',
                 clusterpoints[1][:,0],clusterpoints[1][:,1],'or',
                 clusterpoints[2][:,0],clusterpoints[2][:,1],'og')
        ax1.hist(list(clusterpoints[0][:,0])+list(clusterpoints[1][:,0])+list(clusterpoints[2][:,0]),bins=25,range=pxrange)
        ax3.hist(list(clusterpoints[0][:,1])+list(clusterpoints[1][:,1])+list(clusterpoints[2][:,1]),bins=25,range=pyrange,orientation='horizontal')
    if len(clusterpoints)==2:
        ax2.plot(clusterpoints[0][:,0],clusterpoints[0][:,1],'ob',
                 clusterpoints[1][:,0],clusterpoints[1][:,1],'or')
        ax1.hist(list(clusterpoints[0][:,0])+list(clusterpoints[1][:,0]),bins=25,range=pxrange)
        ax3.hist(list(clusterpoints[0][:,1])+list(clusterpoints[1][:,1]),bins=25,range=pyrange,orientation='horizontal')        
    if len(missingXY)>0:
        ax2.plot(missingXY[:,0],missingXY[:,1],'ok',markerfacecolor='none')
    #medianXY=np.array([(FEFdatadict[ts][0],math.log10(max(FEFdatadict[ts][1],0.1))) for ts in FEFdatadict if ts in medianlist])
    #plt.plot(medianXY[:,0],medianXY[:,1],'sm',markersize=8)
    ax2.axis(pxrange+pyrange)
    ax1.set_xlim(pxrange)
    ax3.set_ylim(pyrange)
    ax2.set_xlabel('Fractional Edge Weight')
    ax2.set_ylabel('log Coverage')
    ax2.set_title('%s:%s:%d-%d'%(genechredge[0],genechredge[1],genechredge[2][0],genechredge[2][1]))
    plt.savefig(imagefilename)


def makejitter(y1list,y2list,pointannotation,ylabel,xlabel,title,imagefilename):
    fig = plt.figure(figsize=(8, 8))
    x1list=[0.8+0.4*x[0] for x in np.random.rand(len(y1list),1).tolist()]
    x2list=[1.8+0.4*x[0] for x in np.random.rand(len(y2list),1).tolist()]
    plt.boxplot([y1list,y2list],widths=0.4)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.plot(x1list,y1list,'b.')
    plt.plot(x2list,y2list,'r.')
    for label,x,y in pointannotation:
        plt.annotate(label,xy=(x,y),xytext=(80,-40),
                     textcoords = 'offset points', ha = 'right', va = 'bottom',
                     bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    plt.axis([0.2,3.0,0.0,1.0])
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_visible(False)
    plt.savefig(imagefilename)

def arrowplot(title,xlabellist,ylabel,y1y2data,sigflglist,wtflglist,imagefilename):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.set_xticks(range(1,len(xlabellist)+1))
    ax.set_xticklabels(xlabellist)
    for i in range(len(y1y2data)):
        if sigflglist[i]==0: 
            plt.plot(i+1,y1y2data[i][0],'bo',markerfacecolor='none')
            plt.plot(i+1,y1y2data[i][1],'rs',markerfacecolor='none')
        else:
            plt.plot(i+1,y1y2data[i][0],'bo')
            plt.plot(i+1,y1y2data[i][1],'rs')  
        if y1y2data[i][0]!=y1y2data[i][1]:
            if wtflglist[i]==1:    
                plt.annotate("",
                            xy=(i+1,y1y2data[i][0]), xycoords='data',
                            xytext=(i+1,y1y2data[i][1]), textcoords='data',
                            arrowprops=dict(arrowstyle="<|-",linewidth=3,
                            connectionstyle="arc3"),
                            )   
            else:
                plt.annotate("",
                            xy=(i+1,y1y2data[i][0]), xycoords='data',
                            xytext=(i+1,y1y2data[i][1]), textcoords='data',
                            arrowprops=dict(arrowstyle="<|-",linewidth=1,
                            connectionstyle="arc3"),
                            )   
    plt.title(title)
    plt.xlabel('')
    plt.ylabel('Splicing Fraction')
    plt.axis([0,len(xlabellist)+1,-0.05,1.05])
    try:
        plt.savefig(imagefilename)
    except:
        message='Problem in writing %s'%imagefilename
        common.printstatus(message,'W',common.func_name())

