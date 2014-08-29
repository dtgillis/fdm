"""
    This module has all the functions used by other modules of actg package
"""

import os
import math
import time
import sys
import random
import numpy as np
import re
import csv
import scipy.stats


def func_name():
    """
    Returns the name of function its called from
    """
    return sys._getframe(1).f_code.co_name

def printstatus(message,message_type,func_nm,verbose_flg=0):
    """
    Prints the Status/Warning/Error message along with the function and time of the message 
    """
    if message_type=='E':
        print '%s : Error   : %30s : %s'%(time.asctime(),func_nm,message)
    if message_type=='W' and verbose_flg:
        print '%s : Warning : %30s : %s'%(time.asctime(),func_nm,message)
    if message_type=='S':
        print '%s : Status  : %30s : %s'%(time.asctime(),func_nm,message)     
    if message_type=='F':
        print '%s : Fatal Error  : %30s : %s'%(time.asctime(),func_nm,message) 
        sys.exit()
        

def fl2str(inlist,fmt='%8.4f'):
    """
    Prints a list of floats as comma separated string with given float format
    """
    return _tfl2str(inlist,fmt)[:-1]

def _tfl2str(inlist,fmt='%8.4f'):
    """
    A recursive function used by function fl2str
    """
    outstr='['
    for el in inlist:
        if isinstance(el,list):
            outstr+=_tfl2str(el,fmt)
        else:
            outstr+=fmt%el+','
    outstr=outstr[:-1]+'],'
    return outstr

def str2fl(str,type='float'):
    """
    Converts a list string (int or float) to a list
    """
    str=str.replace(' ','')
    strlist=str[1:-1].split(',')
    if type=='float':
        list=[float(elem) for elem in strlist]
    if type=='int':
        list=[int(elem) for elem in strlist]        
    return list

def uniquify(list):
   """
   Converts a list to a list of unique elements of the list
   """
   keys = {}
   for e in list:
       keys[e] = 1
   return keys.keys()

def findheaderdict(headerline,coltype=['wt1','wt2','flow']):
    """
    Find the indices of columns in a header line
    Input
    """
    headerdict={}
    for colname in headerline:
        if ':' in colname:
            col=colname.split(':')
            if col[0] not in headerdict:
                headerdict[col[0]]=[-1 for i in range(len(coltype))]
            if col[1] in coltype:
                headerdict[col[0]][coltype.index(col[1])]=headerline.index(colname)
    for h in headerdict.keys():
        if min(headerdict[h])==-1:
            del headerdict[h] 
    return headerdict   

def median(list):
    list.sort()
    l=len(list)
    if l%2==0:
        return 1.0/2*(list[l/2]+list[l/2-1])
    else:
        return list[l/2]

def shannon_entropy(V):
    """
    Computes shannon entropy of a vector
    """
    E = 0
    for v in V:
        if v<1E-10:
            v=1E-10
        E -= v * math.log(v,2)
    return E   

def JSD(V1,V2):
    """
    Computes Jensen-Shannon Divergence between two vectors of same size
    """
    if len(V1)!= len(V2):
        printstatus('Vectors are of unequal length','E',func_name())
        return -1
    else:
        avg_V = []
        for i in range(len(V1)):
            avg_V.append(0.5*(V1[i]+V2[i]))
        JSD = shannon_entropy(avg_V) - 0.5 * (shannon_entropy(V1) + shannon_entropy(V2))
        return JSD

def  DaviesBouldinIndex(pointlist):
    '''
    The sublists are points in the cluster 
    '''
    if [] in pointlist:
        pointlist.remove([])
    if [] in pointlist:
        pointlist.remove([])        
    numclusters=len(pointlist)
    if isinstance(pointlist[0][0],list):
        dim=len(pointlist[0][0])
    else:
        dim=1
    centroidlist=[]
    for sublist in pointlist:
        if dim>1:
            centroid=[0 for i in range(dim)]
            for i in range(dim):
                for point in sublist:
                    centroid[i]+=point[i]/(1.0*len(sublist))
        else:
            centroid=1.0*sum(sublist)/len(sublist)
        centroidlist.append(centroid)
    
    scatterlist=[]
    for i in range(numclusters):
        sublist=pointlist[i]
        centroid=centroidlist[i]
        scatter=0.0
        npts=1.0*len(sublist)
        for point in sublist:
            scatter+=1/npts*lndistance(point,centroid)
        scatterlist.append(scatter)
    Ri=[0 for i in range(numclusters)]
    R=[Ri for i in range(numclusters)]
    for i in range(numclusters):
        for j in range(numclusters):
            if i!=j:
                R[i][j]=(scatterlist[i]+scatterlist[j])/(lndistance(centroidlist[i],centroidlist[j]))
    #D is average over Ri
    if numclusters>1:
        D=[1.0/(numclusters-1)*sum(R[i]) for i in range(numclusters)]
        DBIndex=1.0/numclusters*sum(D)
    else:
        return (100.0,[],[])
        
    mincentroiddist=2.0
    for i in range(numclusters-1):
        for j in range(i+1,numclusters):
            dist=lndistance(centroidlist[i],centroidlist[j])
            if dist<mincentroiddist:
                mincentroiddist=dist    
                
    return (DBIndex,mincentroiddist,centroidlist)

    def DEL_DBIndexpair(pt1,pt2):
        '''
        Assumes pt1 to be smaller list
        '''
        cqual,midfdm,centroidlist=self._DaviesBouldinIndex([pt1,pt2])
        fdmlist=[lndistance(pt1[i],pt2[i]) for i in range(len(pt1))]
        minfdm=min(fdmlist)
        medfdm=median(fdmlist)
        return cqual,minfdm,medfdm

 
def toss(pct=0.5):
    x=random.random()
    if x<pct:
        return 1
    else:
        return 0

def randgroup(list,numgroup):
    '''
    Used in two classes
    '''
    inlist=list[0:]
    if len(inlist)>0 and len(inlist)<numgroup:
        inlist=inlist*numgroup
    random.shuffle(inlist)
    newlist = []
    splitsize = 1.0*len(inlist)/numgroup
    for i in range(numgroup):
        newlist.append(inlist[int(round(i*splitsize)):int(round((i+1)*splitsize))])
    return newlist

def getlistidx(list,element): 
    for i in range(len(list)):
        if list[i] > element:
            return i-1
    return len(list)-1

def mergefilelist(filelist,outcols,outfileptr,keycollist=[0],sep='\t'):
    '''
    All files have the same keys in first column in Order
    No key matching performed
    '''
    fstrlist=[]
    for file in filelist:
        fstrlist.append(open(file))
    for linetxt in fstrlist[0]:
        if len(linetxt)==0:
            continue
        line=linetxt.rstrip('\n').split(sep)
        outline=[]
        for keycol in keycollist:
            outline.append(line[keycol])
        for col in outcols:
            if col<=len(line)-1:
                outline.append(line[col])
            else:
                outline.append('_NA_')
        #print line,outline
        for fstr in fstrlist[1:]:
            linetxt=fstr.readline()
            line=linetxt.rstrip('\n').split(sep)
            found=1
            for keycol in keycollist:
                if line[keycol]!=outline[keycol]:
                    found=0
            if found==0:
                for col in outcols:
                    outline.append('_NA_')
            else:
                for col in outcols:
                    if col<=len(line)-1:
                        outline.append(line[col])
                    else:
                        outline.append('_NA_')
        outlinetxt=sep.join(outline)
        #print outline
        outfileptr.write('%s\n'%outlinetxt)

def mergedict(a, b, path=None):
    "merges b into a"
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                mergedict(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a

def mergedictlist(dictlist):
    out={}
    for dicti in dictlist:
        mergedict(out,dicti)
    return out

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def concatprefixfiles(outfile,folder,prefix,headerflg=1):
    dirfiles = [ f for f in os.listdir(folder) if os.path.isfile('%s/%s'%(folder,f)) and f.startswith(prefix) ]
    sortedfiles=natural_sort(dirfiles)
    concatfiles(outfile,sortedfiles,headerflg)

def getgroupshufflepvalue(group1,group2,numiter=100000):
    """
        Computes how well 1's are separated from 0s
        group1, group2: List of 1s and 0s only
    """
    diff=sum(group1)-sum(group2)
    rnddifflist=[]
    concatgroup=group1+group2
    for i in range(numiter):
        randidx1=random.sample(range(len(concatgroup)),len(group1))
        randgroup1=[concatgroup[j] for j in range(len(concatgroup)) if j in randidx1]
        randgroup2=[concatgroup[j] for j in range(len(concatgroup)) if j not in randidx1]
        rnddiff=sum(randgroup1)-sum(randgroup2)
        rnddifflist.append(rnddiff)
    rnddifflist.sort()
    for i in range(numiter):
        if rnddifflist[i]>=diff:
            break
    if rnddifflist[i]<diff:
        i+=1
    return 1-i*1.0/numiter

def combinepvalues(pvaluelist,method):
    '''
    method='fisher','liptak'
    '''
    eps=0.0001
    plist=[min(1-eps,max(eps,x)) for x in pvaluelist]
    if method=='fisher':
        fisherchi2=min(114,-2.0*sum([math.log(p) for p in plist]))
        fisherpval=1-scipy.stats.chi2.cdf(fisherchi2, 2*len(plist))
        return fisherpval
    if method=='liptak':
        w=[1.0/len(plist) for x in plist]
        Zilist=[scipy.stats.norm.ppf(1-p) for p in plist]
        Z=min(8.2,sum([x*y for x,y in zip(Zilist,w)])/math.sqrt(sum([wi*wi for wi in w])))
        liptakpval=1-scipy.stats.norm.cdf(Z)
        return liptakpval     

def concatfiles(outfile,filelist,headerflg=1):
    fout=open(outfile,'w')
    filecnt=0
    for infile in filelist:
        filecnt+=1
        fin=open(infile)
        if filecnt!=1:
            fin.readline()
        for lntxt in fin:
            fout.write(lntxt)
        fin.close()

def ifthenelse(ifclause,thenclause,elseclause):
    if ifclause:
        return thenclause
    else:
        return elseclause

def csv2html(csvfile,htmlfile,delim='\t'):
    reader = csv.reader(csvfile,delimiter=delim)
    rownum = 0
    htmlfile.write('<table border="1">\n')
    for row in reader:
        if rownum == 0:
            htmlfile.write('<tr>') 
            for column in row:
                htmlfile.write('<th>' + column + '</th>')
            htmlfile.write('</tr>\n')
        else:
            htmlfile.write('<tr>')    
            for column in row:
                htmlfile.write('<td>' + column + '</td>')
            htmlfile.write('</tr>\n')
        rownum += 1
    htmlfile.write('</table>\n')
    print "Created " + str(rownum) + " row table." 

def normalize_vector(V):
    Z=1E-20
    Vector = V[0:]
    if sum(Vector)==0:
        for i in range(len(Vector)):
            Vector[i]=Z
            
    sumv = sum(Vector)
    for i in range(len(Vector)):
        Vector[i]=round(Vector[i]*1.0/sumv,4)
    return Vector


def lndistance(point1,point2,n=1):
    if isinstance(point1,list):
        if n==1:
            distance=1.0*sum([abs(point1[i]-point2[i]) for i in range(len(point1))])
        elif n==2:
            distance=math.sqrt(sum([(point1[i]-point2[i])**2 for i in range(len(point1))]))
    else:
        distance=abs(point1-point2)
    return max(distance,0.000001)
             

def islandlist2chr1Mdict(islandlist):
    chr1Mdict={}
    for island in islandlist:
        chrnm,startb,endb, gene,strand = island 
        #chrnm='chr%s'%chrnm
        if chrnm not in chr1Mdict:
            chr1Mdict[chrnm]={}
        start1M=startb/1000000
        end1M=endb/1000000
        for i in range(start1M,end1M+1):
            if i not in chr1Mdict[chrnm]:
                chr1Mdict[chrnm][i]=[]
            chr1Mdict[chrnm][i].append(gene)
    return chr1Mdict

def intv2gene(chrpos,chr1Mdict,islanddict):
    chrnm,startb,endb=chrpos
    chrnode1=(chrnm,startb)
    genelist1=node2gene(chrnode1,chr1Mdict,islanddict)
    chrnode2=(chrnm,endb)
    genelist2=node2gene(chrnode2,chr1Mdict,islanddict)    
    pos1andpos2=[x for x in genelist1 if x in genelist2]
    pos1minuspos2=[x for x in genelist1 if x not in genelist2]
    pos2minuspos1=[x for x in genelist2 if x not in genelist1]
    return [genelist1,genelist2,pos1andpos2,pos1minuspos2,pos2minuspos1]

def node2gene(chrnode,chr1Mdict,islanddict):
    chrnm,pos=chrnode
    pos1M=pos/1000000
    genelist=[]; tgenelist=[]
    if chrnm in chr1Mdict:
        if pos1M in chr1Mdict[chrnm]:
            tgenelist=chr1Mdict[chrnm][pos1M] 
    outgenelist=[]
    for gene in tgenelist:
        if islanddict[gene][1]<pos and islanddict[gene][2]>pos:
            outgenelist.append(gene)
    return outgenelist

def cigar2juncexon(cigarstr):
    '''
    Alert: Junction Start should have been -1; two downstream user defs correct it
    '''
    covlist = [int(x) for x in cigarstr.rstrip('M').replace('M', 'N').split('N')]
    jtuplist=[]
    elist=[]
    offset=0   
    if len(covlist)==1:
        elist.append((offset,offset+covlist[0]))
    else:
        while len(covlist)!=1:
            elist.append((offset,offset+covlist[0]))
            jtuplist.append((offset+covlist[0],offset+covlist[0]+covlist[1]))
            offset+=covlist[0]+covlist[1]
            covlist=covlist[2:]
        elist.append((offset,offset+covlist[0]))
    return (jtuplist,elist)    

def exonlist2cigar(exonlist):
    cigarstr='%dM'%(exonlist[0][1]-exonlist[0][0]+1)
    exon2gap=[(x[1][0]-x[0][1]-1,x[1][1]-x[1][0]+1) for x in zip(exonlist[:-1],exonlist[1:])]
    txsize=exonlist[0][1]-exonlist[0][0]+1
    for gap in exon2gap:
        cigarstr+='%dN%dM'%gap
        txsize+=gap[1]
    return (exonlist[0][0],cigarstr,txsize)
    
                  
def exonbaseoverlap(exonlist1,exonlist2,exonwtlist2):
    '''Output is bases covered for each exon1
    '''
    i=0;j=0;size=0
    covlist=[]
    while i<len(exonlist1) and j<len(exonlist2):
        exon1=exonlist1[i]
        exon2=exonlist2[j]
        currsize = max(0,min(exon1[1],exon2[1])-max(exon1[0],exon2[0])+1)
        size+=currsize*exonwtlist2[j]         
        if exon1[1]<exon2[1]:
            i+=1
            covlist.append(size)
            size=0
        elif exon1[1]>exon2[1]:
            j+=1
        else:
            i+=1;j+=1
            covlist.append(size)
            size=0            
    while i<len(exonlist1):
        covlist.append(size)
        size=0;i+=1
    return covlist

def convgenericcigar(cigar):
    ocig = cigar[0:]
    cigval = ocig.replace('N','M').replace('I','M').replace('D','M').replace('S','M').replace('S','M').replace('P','M').split('M')
    cigval = cigval[:-1]
    cigsep = cigar[0:]
    for val in cigval:
        cigsep=cigsep.replace(val,'',1)
    outcigar=''
    currsep=''; sum=0
    for i in range(len(cigsep)):
        if cigsep[i] in ['M','D']:
            currsep='M'
            sum+=int(cigval[i])
        elif cigsep[i]=='N':
            outcigar+='%dM%sN'%(sum,cigval[i])
            currsep=''; sum=0
    outcigar+='%dM'%sum
    return outcigar

def cigarsubstr(cigarstr,substrstart,substrsize):
    cigartuple=cigarstr.replace('N','M').split('M')[:-1]
    exonic=[int(cigartuple[i]) for i in range(len(cigartuple)) if i%2==0]
    intronic=[int(cigartuple[i]) for i in range(len(cigartuple)) if i%2==1]
    substrend=substrstart+substrsize
    subcigarstr=''
    cutleft=substrstart
    foundexon=0
    for i in range(len(exonic)):
        if foundexon==1:
            if cutleft>exonic[i]:
                subcigarstr+='%dM'%(exonic[i])
                if i<len(intronic):
                    subcigarstr+='%dN'%(intronic[i])
            elif cutleft>0:
                subcigarstr+='%dM'%cutleft
            else:
                break
            cutleft-=exonic[i]           
        if foundexon==0:
            if cutleft>=exonic[i]:
                cutleft-=exonic[i]
            else:
                foundexon=1
                subcigarstr='%dM'%min(exonic[i]-cutleft,substrsize)
                cutleft=substrsize-min(exonic[i]-cutleft,substrsize)
                if cutleft>0:
                    if i<len(intronic):
                        subcigarstr+='%dN'%intronic[i]
        #print subcigarstr
    return subcigarstr

def aligntuple2coverage(aligntuple,exonlist):
    refbase=int(aligntuple[1])
    gencigar=convgenericcigar(aligntuple[2])
    juncexon=cigar2juncexon(gencigar)
    exoncovlist=[]
    for exon in juncexon[1]:
        exoncovlist.append([refbase+exon[0], refbase+exon[1]])
    junccovlist=[]
    for junc in juncexon[0]:
        junccovlist.append([refbase+junc[0]-1, refbase+junc[1]])   
    exoncovwtlist=[1 for elem in exoncovlist]
    exoncoveragelist= exonbaseoverlap(exonlist,exoncovlist,exoncovwtlist)
    #print 'junc',junccovlist
    return (exoncoveragelist,junccovlist)
                 
def aligntuples2coverage(aligntuples,exonlist,splicelist):
    '''aligntuples = (chr,refbase,cigar)
       genomeintervals=(exonic,intronic): exonic first exon gets priority of base
    '''
    exoncoveragelist=[0]*len(exonlist)
    splicecoveragelist=[0]*len(splicelist)
    for aligntuple in aligntuples:
        coveragelist=aligntuple2coverage(aligntuple,exonlist)
        texoncoveragelist=coveragelist[0]
        exoncoveragelist=[i+j for i,j in zip(texoncoveragelist, exoncoveragelist)]
        tsplicecoveragelist=coveragelist[1]
        for splice in tsplicecoveragelist:
            try:
                idx=splicelist.index(splice)
            except:
                idx=-1
            if idx!=-1:
                splicecoveragelist[idx]+=1
    for i in range(len(exoncoveragelist)):
        exoncoveragelist[i]=1.0*exoncoveragelist[i]/(exonlist[i][1]-exonlist[i][0]+1)
    return (exoncoveragelist,splicecoveragelist)

        

def splitexonwithinsplice(graphstruct):
    exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist=graphstruct
    for  splice in splicelist:
        coverexon=[]
        if [splice[0],splice[1]] in exonlist+intronlist:
            coverexon=[splice[0],splice[1]]
        elif [splice[0]+1,splice[1]] in exonlist+intronlist:
            coverexon=[splice[0]+1,splice[1]]
        elif [splice[0]+1,splice[1]-1] in exonlist+intronlist:
            coverexon=[splice[0]+1,splice[1]-1]
        elif [splice[0],splice[1]-1] in exonlist+intronlist:
            coverexon=[splice[0],splice[1]-1]        
        if coverexon!=[]:
            if coverexon in exonlist:
                exonlist.remove(coverexon)
                if coverexon[1]-coverexon[0]>=2:
                    dummyleftnode=(coverexon[0]+coverexon[1])/2
                    dummyrightnode=dummyleftnode+1
                    exonlist.append([coverexon[0],dummyleftnode])
                    exonlist.append([dummyrightnode,coverexon[1]])
                    novelnodelist.append(dummyleftnode)
            if coverexon in intronlist:
                intronlist.remove(coverexon)
                if coverexon[1]-coverexon[0]>=2:
                    dummyleftnode=(coverexon[0]+coverexon[1])/2
                    dummyrightnode=dummyleftnode+1
                    intronlist.append([coverexon[0],dummyleftnode])
                    intronlist.append([dummyrightnode,coverexon[1]])
                    novelnodelist.append(dummyleftnode)
    exonlist.sort()
    intronlist.sort()
    novelnodelist.sort()
    return (exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist)

def divdict2graphstruct(divdict):
    '''
    divdict=position:[[incoming/outgoing=0,1,exonstart=0=no,1=yes,2=start transcript,3=end transcript][exon, flowlist]]
    exonstart=0=no,1=yes,2=start transcript and exonstart,3=end transcript and exonstart,4=start transcript and no exonstart(insplice),
    5=end transcript and no exonstart outsplice 
    flowlist=[exon/splicelist]
    wtgraphstruct=(exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist)
    move bam file to class bam
    '''
    exonlist=[];splicelist=[]
    for pos in divdict.keys():
        outflg,exonflg=divdict[pos][0]
        exon,flowlist=divdict[pos][1]
        if exon not in exonlist:
            exonlist.append(exon)
        if exonflg in [1,2,3]:
            if flowlist[0] not in exonlist:
                exonlist.append(flowlist[0])
        else:
            if flowlist[0] not in splicelist:
                splicelist.append(flowlist[0])
        for splice in flowlist[1:]:
            if splice not in splicelist:
                splicelist.append(splice)
    
    exonlist.sort(); splicelist.sort()
    return [exonlist,[],splicelist,[],[],[]]
           
def graphstruct2divdict(graphstruct):
    '''startnode is always incoming = 2,4
    endnode is always outgoing = 3,5
    divdict=position:[incoming/outgoing=0,1,
    exonstart=0=no,1=yes,2=start transcript and exonstart,3=end transcript and exonstart,4=start transcript and no exonstart(insplice),
    5=end transcript and no exonstart outsplice 
    flow exon/splicelist]
    rewrite as outflow and then inflow
    warning only if both inflow and outflow at same point
    chromosome inclusion
    '''
    #print '******',func_name, __name__
    exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist=graphstruct
    divdict={}
    #print len(exonlist)
    #Find outgoing  
    for i in range(len(exonlist)-1):
        if exonlist[i][1]-exonlist[i][0]<5:
            continue
        oflowlist=[]
        if (exonlist[i+1][0]-exonlist[i][1] in [0,1]):
            oflowlist.append(exonlist[i+1])
            exonstartflag=1
        else:
            exonstartflag=0
        for splice in splicelist:
            if splice[0]==exonlist[i][1]:
                oflowlist.append(splice)
        exontxflag=exonstartflag
        if exonlist[i][1] in endnodelist or exonlist[i+1][0] in endnodelist:
            if exonstartflag==1:
                exontxflag=3
            else:
                exontxflag=5
        if len(oflowlist)>1 or exontxflag==3:
            divdict[exonlist[i][1]]=[[1,exontxflag],[exonlist[i],oflowlist]]    
    for i in range(len(exonlist)-1):
        if exonlist[i+1][1]-exonlist[i+1][0]<5:
            continue
        message=''
        iflowlist=[]
        if (exonlist[i+1][0]-exonlist[i][1] in [0,1]):
            iflowlist.append(exonlist[i])
            exonstartflag=1
        else:
            exonstartflag=0            
        for splice in splicelist:
            if splice[1]==exonlist[i+1][0]:
                iflowlist.append(splice)
        exontxflag=exonstartflag
        if exonlist[i][1] in startnodelist or exonlist[i+1][0] in startnodelist:
            if exonstartflag==1:
                exontxflag=2
            else:
                exontxflag=4
        if (len(iflowlist)>1 or exontxflag==2):
            if exonlist[i+1][0] in divdict:
                message='Divergence Node exists at %d: %s \n'%(exonlist[i+1][0],str(divdict[exonlist[i+1][0]]))
                message+='Trying to insert new divergence: %s'%str([[0,exontxflag],[exonlist[i+1],iflowlist]])
                printstatus(message,'W',func_name())  
            else:
                divdict[exonlist[i+1][0]]=[[0,exontxflag],[exonlist[i+1],iflowlist]]            
    return divdict
   
def divdictnalignlines2flow(divdict,alignlines):
    '''
    This needs to change similar to wtgraphstruct2flow
    get 
    '''
    divexonlist=[]; divsplicelist=[]
    for position in divdict.keys():
        outgoingflg,exonstartflg=divdict[position][0]
        exon,flowlist=divdict[position][1]
        if exon not in divexonlist:
            divexonlist.append(exon)
        if exonstartflg in [1,2,3]:
            if flowlist[0]  not in divexonlist:
                divexonlist.append(flowlist[0])
            elif flowlist[0]  not in divsplicelist:
                divsplicelist.append(flowlist[0])
        for splice in flowlist[1:]:
            if splice not in divsplicelist:
                divsplicelist.append(splice)
    divexonlist.sort()
    divsplicelist.sort()
    coveragelist=aligntuples2coverage(aligntuples,divexonlist,divsplicelist)
    
    flowdict={}
    for position in divdict.keys():
        outgoingflg,exonstartflg=divdict[position][0]
        exon,flowlist=divdict[position][1]
        wt1=coveragelist[0][divexonlist.index[exon]]
        flowvec=[]
        if exonstartflg in [2,3]:
            wtoth=coveragelist[0][divexonlist.index[flowlist[0]]]
            if exonstartflg==2:
                flowvec=[wt1,max(wtoth-wt1,0)]
                if wtoth-wt1<0:
                    message='Flow decreases at transcript start %s:%s; Weight Start=%10.4f, Before=%10.4f'%(str(exon),fl2str(flowlist[0]),wt1,wtoth)
                    printstatus(message,'W',func_name())     
            else:
                flowvec=[wt1,max(wt1-wtoth,0)]
                if wt1-wtoth<0:
                    message='Flow increases at transcript end %s:%s; Weight End=%10.4f, After=%10.4f'%(str(exon),fl2str(flowlist[0]),wt1,wtoth)
                    printstatus(message,'W',func_name())  
        if exonstartflg==1:
            flowvec.append(coveragelist[0][divexonlist.index[flowlist[0]]])
        elif exonstartflg==0:
            flowvec.append(coveragelist[1][divsplicelist.index[flowlist[0]]])
        for flowedge in flowlist[1:]:
            flowvec.append(coveragelist[1][divsplicelist.index[flowlist[0]]])
        wt2=sum(flowvec) 
        nflowvec=normalize_vector(flowvec)
        flowdict[position]=[[outgoingflg,exonstartflg],[[wt1,wt2],nflowvec]]    
    return flowdict

def make_ucsctrack(tslist,run_prefix,rootdir,outdir,httpdir):
    if not os.path.isdir(outdir):
        os.system('mkdir -p %s'%outdir)
    outfilename='%s/%s_tracks.txt'%(outdir,run_prefix)
    fout=open(outfilename,'w')
    fout.write('browser position chr20:24986874-25039616\n')
    fout.write('browser hide all\n')
    fout.write('browser pack refGene encodeRegions mrna\n')
    for ts in tslist:
        cmd='ln -s %s/%s.bb %s/%s.bb'%(cfg.datadir_dict[ts],ts,outdir,ts)
        os.system(cmd)
        cmd='ln -s %s/%s.bw %s/%s.bw'%(cfg.datadir_dict[ts],ts,outdir,ts)
        os.system(cmd)
        bbtrack='track type=bigBed name="%s splices" visibility=pack color=200,100,0 altColor=0,100,200 useScore=1 bigDataUrl=%s/%s.bb'%(ts,httpdir,ts)
        fout.write('%s\n'%bbtrack)
        bwtrack='track type=bigWig name="%s Coverage" visibility=full color=200,100,0 altColor=0,100,200 bigDataUrl=%s/%s.bw'%(ts,httpdir,ts)
        fout.write('%s\n'%bwtrack)
    
    fdmresultdir='%s/result/fdm/%s'%(rootdir,run_prefix)
    fdmfilelist=[filename for filename in os.listdir(fdmresultdir) if filename.startswith(run_prefix) and filename.endswith('.bw')]
    fdmlist.sort()
    for fdmfile in fdmfilelist:
        cmd='ln -s %s/%s %s/%s'%(fdmresultdir,fdmfile,outdir,fdmfile)
        os.system(cmd)
        fdmset1=fdmfile.split('__')[1]
        fdmset2=fdmfile.split('__')[2]
        if fdmfile.split('__')[3]=='__fdm.bdg':
            ftrack='track type=bigWig name="FDM %s: %s-%s" visibility=full color=0,200,100 altColor=0,100,200 bigDataUrl=%s/%s'%(run_prefix,fdmset1,fdmset2,httpdir,fdmfile)
            fout.write('%s\n'%ftrack)
        elif fdmfile.split('__')[3]=='__pvl.bdg':
            ptrack='track type=bigWig name="p-value %s: %s-%s" visibility=full color=0,200,100 altColor=0,100,200 bigDataUrl=%s/%s'%(run_prefix,fdmset1,fdmset2,httpdir,fdmfile)
            fout.write('%s\n'%ptrack)

def fdmstat2(pvaldict,groups,nperm=1000,pval1=0.05,pval2=0.05):
    '''
    pval= {(1,2):pval}
    group [[],[]]
    '''
    tslist=groups[0]+groups[1]
    sdict={}
    for i in range(len(tslist)-1):
        for j in range(i+1,len(tslist)):
            #print i,j
            if (tslist[i],tslist[j]) in pvaldict:
                if pvaldict[(tslist[i],tslist[j])]<pval1:
                    sdict[(i,j)]=1
                else:
                    sdict[(i,j)]=0
            elif (tslist[j],tslist[i]) in pvaldict:
                if pvaldict[(tslist[i],tslist[j])]<pval1:
                    sdict[(i,j)]=1
                else:
                    sdict[(i,j)]=0
            else:
                message='FDM not run for %s:%s'%(tslist[i],tslist[j])
                printstatus(message,'W',func_name())
    group1=[tslist.index(ts) for ts in groups[0]]
    group2=[tslist.index(ts) for ts in groups[1]]
    ow1=0
    for i in group1:
        for j in group1:
            ow1+=sdict.get((i,j),0)
    ow2=0
    for i in group2:
        for j in group2:
            ow2+=sdict.get((i,j),0)    
    ow1w2=0
    for i in group1:
        for j in group2:
            ow1w2+=sdict.get((i,j),0)
            ow1w2+=sdict.get((j,i),0)
    odiff=-1*(ow1w2-ow1-ow2)
    #print w1,w2,w1w2
    stat2null=[]
    cntts=len(tslist)
    for k in range(nperm):
        group1=random.sample(range(cntts),cntts/2)
        group2=[i for i in range(cntts) if i not in group1]
        w1=0
        for i in group1:
            for j in group1:
                w1+=sdict.get((i,j),0)
        w2=0
        for i in group2:
            for j in group2:
                w2+=sdict.get((i,j),0)    
        w1w2=0
        for i in group1:
            for j in group2:
                w1w2+=sdict.get((i,j),0)
                w1w2+=sdict.get((j,i),0)
        diff=-1*(w1w2-w1-w2)
        stat2null.append(diff)
    stat2null.sort()
    pval=(getlistidx(stat2null,odiff)+1.0)/len(stat2null)
    return (pval,odiff,ow1,ow2,ow1w2)   
         