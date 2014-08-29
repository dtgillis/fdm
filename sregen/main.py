#!/usr/bin/python
#sregen

import cPickle
import os
import math
import time
import sys
import random
import numpy as np

import ConfigParser
from optparse import OptionParser

sys.path.insert(0, '..')
import common
import gtffile
import plot

def RNAmetasource2source(parameterlist):
    """
    Output:
        Gene Dict: Transcript Dict: chromosome, strand, startbase, cigar, txsize
    """
    # empirical parameter
    genewithnovelskippedexonpct=50
    metasourcefilename=parameterlist[0][0]
    filetype=parameterlist[1][0]
    if parameterlist[1][1]=='no_novel_transcript':
        novel_skipped_exon_flg=0
    else:
        novel_skipped_exon_flg=1
    numgenes=int(parameterlist[1][2])
    genelist=parameterlist[1][3]
    
    if filetype=='gtf':
        metasource=gtffile.gtfFile(metasourcefilename)
        genetranscriptdict=metasource.getgenetranscriptdict()
    elif filetype=='pck':
        genetranscriptdict=cPickle.load(open(metasourcefilename))
    
    allgenes=genetranscriptdict.keys()
    if len(genelist)==0:
        if numgenes==0:
            numgenes=len(genetranscriptdict)
        else:
            numgenes=min(numgenes,len(genetranscriptdict))
        choosegenes=random.sample(allgenes,numgenes)
    else:
        choosegenes=genelist[0:]
    
    #print len(genelist),len(choosegenes)
    
    outtranscriptdict={}
    for gene in choosegenes:
        txdict=genetranscriptdict[gene]
        txexonslist=[txdict[tx][2] for tx in txdict]
        if novel_skipped_exon_flg==1 and common.toss(genewithnovelskippedexonpct/100.0):
            newtxflg=0
            trys=0
            while not newtxflg and trys<10:
                trys+=1
                tx=random.choice(txdict.keys())
                txexons=txdict[tx][2]
                if len(txexons)>2:
                    skipped=random.randint(1,len(txexons)-2)
                    newtx=txexons[:skipped-1]+txexons[skipped:]
                if newtx not in txexonslist:
                    newtxflg=1
                    txdict['%s_skip'%tx]=[txdict[tx][0],txdict[tx][1],newtx]
        outtranscriptdict[gene]={}
        for tx in txdict:
            cigartup=common.exonlist2cigar(txdict[tx][2])
            outtranscriptdict[gene][tx]=[txdict[tx][0],txdict[tx][1],cigartup[0],cigartup[1],cigartup[2]]
    return outtranscriptdict                    

def getsimulatedgenereadcounts(numgenes,numreads):
    """
    Computes number of simulated reads for all genes
    input:
        numgenes=total number of genes for which reads are generated
        numreads==total number of reads generated
    output:
        a list of size numgenes containing the number of reads for each gene
    """
    # empirical parameters 
    lognormmu=3.4
    lognormsigma=0.95
    lornormrange=[0.1,5.1]
    
    logreadstemp1=[random.gauss(lognormmu,lognormsigma) for i in range(2*numgenes)]
    logreadstemp2=[x for x in logreadstemp1 if x>lornormrange[0] and x<lornormrange[1]]
    logreads=random.sample(logreadstemp2,numgenes)
    reads=[pow(10,x) for x in logreads]
    genereads=[int(x*numreads/sum(reads)) for x in reads]
    
    return genereads
    
def getsimulatedtxexpressionpct(numtx):
    """
    Computes the simulated relative expressions of transcripts in a gene
    input:
        numtx: number of transcripts
    output:
        a vector of size numtx containing percent share of each transcript
    """
    #empirical parameter
    paramlist=[(0.0,0.4,0.0,0.1),(0.4,07,0.00,0.33),(0.7,0.9,0.0,1.0),(0.9,1.0,0.0,2.0)]
    chooser=random.random()
    for param in paramlist:
        if param[0]< chooser <=param[1]:
            out=[pow(10,random.uniform(param[2],param[3])) for i in range(numtx)]
    
    return [x/sum(out) for x in out]

def getsimulatedoneexpression(genetxcountdict,numreads):
    numgenes=len(genetxcountdict)
    genereadslist=getsimulatedgenereadcounts(numgenes,numreads)
    
    outdict={}
    i=0
    for gene in genetxcountdict:
        v1=getsimulatedtxexpressionpct(genetxcountdict[gene])
        outdict[gene]=[genereadslist[i],v1]
        i+=1
    return outdict
    

def getsimulatedtwoexpressions(genetxcountdict,covdistparam,JSDdistparam,numreads):
    """
    Computes gene and transcript coverages for two datasets
    input:
        genetxcountdict: count of transcripts for each gene, make sure than 
        covdistparam: the parameters describing truncated normal distribution - mu,sigma,left,right
                      of ratio of rpkm of gene between two datasets
        JSDdistparam: the parameters describing uniform plus triangular distribution - left, right, mode
                      of JSD of trascript expression of a gene between two datasets
    output:
        genepairexpdict=gene:genereads1,genereads2,txexp1,txexp2
    """
    numgenes=len(genetxcountdict)
    genereadslist1=getsimulatedgenereadcounts(numgenes,numreads)
    #print covdistparam
    covratiolisttemp1=[random.gauss(covdistparam[0],covdistparam[1]) for i in range(10*numgenes)]
    covratiolisttemp2=[x for x in covratiolisttemp1 if covdistparam[2]<=x<=covdistparam[3]]
    covratiolist=random.sample(covratiolisttemp2,numgenes)
    # 1: gene read count <= 5% of total is acceptable 
    # 2: gene read count >= 60% of total is not acceptable
    # 3: if between 1 and 2, gene read count < 10 times fair share is acceptable
    genefairshare=numreads*1.0/numgenes
    genereadslist2temp1=[genereadslist1[i]*pow(10,covratiolist[i]) for i in range(numgenes)]
    genereadslist2temp2=[min(x,numreads*0.60) for x in genereadslist2temp1]
    genereadslist2temp3=[min(x,max(10*genefairshare,numreads*0.05)) for x in genereadslist2temp2]
    genereadslist2=[max(1,int(x*numreads/sum(genereadslist2temp3))) for x in genereadslist2temp3]
    
    JSDlist=[]
    for gene in genetxcountdict:
        if genetxcountdict[gene]>1:
            if common.toss():
                JSDlist.append(random.uniform(JSDdistparam[0],JSDdistparam[1]))
            else:
                JSDlist.append(random.triangular(JSDdistparam[0],JSDdistparam[1],JSDdistparam[2]))
                
    numbins=50
    JSDbins=np.linspace(0,1,numbins+1)
    JSDbinscount=[]
    for i in range(numbins):
        JSDbinscount.append(len([x for x in JSDlist if JSDbins[i]<x<=JSDbins[i+1]]))
    
    #print JSDbinscount
    
    outdict={}
    genelist=genetxcountdict.keys()
    for i in  range(numgenes):
        gene=genelist[i]
        if genetxcountdict[gene]==1:
            outdict[gene]=[genereadslist1[i],genereadslist2[i],[1.0],[1.0],0.0]
        else:
            binned=0; numattempts=0
            while binned==0:
                numattempts+=1
                v1=getsimulatedtxexpressionpct(genetxcountdict[gene])
                v2=getsimulatedtxexpressionpct(genetxcountdict[gene])
                genejsd=common.JSD(v1,v2)
                binidx=int(genejsd*50)
                if JSDbinscount[binidx]>0:
                    binned=1
                    JSDbinscount[binidx]-=1
                if numattempts>10000 and binidx>1:
                    binned=1
                    JSDbinscount[binidx]-=1
                
            outdict[gene]=[genereadslist1[i],genereadslist2[i],v1,v2,genejsd]     
    return outdict

def getfabases(faptr,startbase,endbase,rowsz=50):
    '''
    rowsz = fa file row size
    '''
    startidx=startbase+int((startbase-1)/rowsz)            #second part is number of \n
    endidx=endbase+int((endbase-1)/rowsz)
    strsize=endidx-startidx+1
    faptr.seek(0)
    faptr.readline()
    faptr.seek(startidx-1, os.SEEK_CUR)
    basestr=faptr.read(strsize)
    basestr=basestr.replace('\n','')
    return basestr 
    
def getRNAtranscriptstring(chrfafileptr,startbase,cigarstr):
    cigarlist=[int(x) for x in cigarstr.replace('M','|').replace('N','|').split('|')[:-1]]
    cigarlist.insert(0,0)
    unsplicedstr=getfabases(chrfafileptr,startbase,startbase+sum(cigarlist))
    cigarcumsum=list(np.cumsum(cigarlist))
    outstr=''.join([unsplicedstr[cigarcumsum[2*i]:cigarcumsum[2*i+1]] for i in range(len(cigarcumsum)/2)])
    #print 'here',cigarstr,cigarlist,len(unsplicedstr),cigarcumsum,len(outstr)
    return outstr.upper()

def getfragments(txstr,numfragments,strand,cutscoredict,fragmentsizerange,txstartbias,numcutmu,numcutsig):
    """
    input:
        txstr: transcript string
        numfragments: number of fragments output
        strand: strand of transcript
        cutscoredict: dict(xmer:score) where score is cuttability at the center of the xmer, default 1
        fragmentsizerange: size of fragments output
        txstartbias: bias of fragments extracted from the start of the transcript
    output:
        fragments: start,end
    """
    outfragments=[]; fragmentscount=0
    if numfragments==0:
        return outfragments
    sttime=time.time()
    #compute cutpoints and probabilities
    scorelist=[1 for i in range(len(txstr))]
    minxmer=min([len(xmer) for xmer in cutscoredict])
    for i in range(len(txstr)-minxmer+1):
        for xmer in cutscoredict:
            if txstr[i:i+len(xmer)]==xmer:
                scorelist[i]=cutscoredict[xmer]
                break
    scorecumlist=[x*1.0/sum(scorelist) for x in np.cumsum(scorelist)]
    scorecumlist.insert(0,0)
    
    numcuts=max(numcutmu,int(random.gauss(numcutmu,numcutsig)))
    
    tries=0
    while fragmentscount<numfragments:
        tries+=1
        txdegradepoint=np.random.exponential(0.005*txstartbias)
        if strand=='+':
            cutrange=[0,int(round(len(txstr)*(1-txdegradepoint)))]
        else:
            cutrange=[int(round(len(txstr)*txdegradepoint)),len(txstr)-1]
        cutpoints=list(np.random.uniform(0,1,numcuts))
        cutpositions=[] 
        #print scorecumlist
        #print cutpoints
        for cutpoint in cutpoints:
            for i in range(cutrange[0],cutrange[1]):
                if cutpoint>scorecumlist[i] and cutpoint<scorecumlist[i+1]:
                    cutpositions.append(i)
        cutpositions.sort()
        cutpositions.insert(0,cutrange[0])
        cutpositions.append(cutrange[1]+1)
        #print 'here',cutpositions
        goodfragments=[[x[0],x[1]-1] for x in zip(cutpositions[:-1],cutpositions[1:]) if fragmentsizerange[0]<=x[1]-x[0]+1<= fragmentsizerange[1]]
        #print 'here',len(goodfragments)
        outfragments+=goodfragments
        fragmentscount+=len(goodfragments)
    if debug_flg==1:
        eltime=time.time()-sttime
        message='Start Time=%s'%sttime
        common.printstatus(message,'S',common.func_name())
        message='Average fragment: %6.2f, Tries per read: %6.2f, Reads per second: %10.4f, Num reads: %d'%(len(txstr)/(numcuts+0.1),tries/(numfragments+0.1),numfragments/eltime,numfragments)
        common.printstatus(message,'S',common.func_name())
    outfragments=random.sample(outfragments,numfragments)
    outfragments.sort()
    return outfragments

def fragment2read(allparameterlist):
    """
    input:
        transcript string, transcript cigar, fragmentlist
    output:
        list of (readstr,readcigar)
    """
    #read_extract_parameters=SE,100,300,400       ; single:size and fragment range
    #read_extract_parameters=PE,100,300,400       ; paired:size and fragment range
    #read_extract_parameters=PB,500,1000          ; pacbio fragment range
    txstr,txcigar,fragmentlist,parameterlist=allparameterlist
    readlist=[]          
    if parameterlist[0]=='SE':
        substrsize=int(parameterlist[1])
        for fragment in fragmentlist:
            if common.toss():
                substrstart=fragment[0]
            else:
                substrstart=fragment[1]-substrsize
            readcigar=common.cigarsubstr(txcigar,substrstart,substrsize)
            readstr=txstr[substrstart:substrstart+substrsize]
            readlist.append([substrstart,readstr,readcigar])
    if parameterlist[0]=='PE':
        substrsize=int(parameterlist[1])
        for fragment in fragmentlist:
            substrstart1=fragment[0]
            readcigar1=common.cigarsubstr(txcigar,substrstart1,substrsize)
            readstr1=txstr[substrstart1:substrstart1+substrsize]
            #readlist.append([substrstart,readstr,readcigar])
            substrstart2=fragment[1]-substrsize
            readcigar2=common.cigarsubstr(txcigar,substrstart2,substrsize)
            readstr2=txstr[substrstart2:substrstart2+substrsize]
            readlist.append([[substrstart1,readstr1,readcigar1],[substrstart2,readstr2,readcigar2]])       
    if parameterlist[0]=='PB':
        for fragment in fragmentlist:
            substrstart=fragment[0]
            substrsize=fragment[1]-fragment[0]
            readcigar=common.cigarsubstr(txcigar,substrstart,substrsize)
            readstr=txstr[substrstart:substrstart+substrsize]
            readlist.append([substrstart,readstr,readcigar])       
    readlist.sort() 
    if parameterlist[0]=='PE':
        pe_readlist=readlist[0:]  
        readlist=[] 
        for pe_read in pe_readlist:
            read1,read2=pe_read
            readlist.append(read1)
            readlist.append(read2)
    return readlist
            
        
def adderrortoread(readstr,readcigar,readqualitydeteriorationrate,indelrate):
    switchdict={'A':'t','T':'a','C':'g','G':'c'}
    outstr=''
    errorrate=0
    indellist=[0]
    indelflaglist=[-1]
    for i in range(len(readstr)):
        if common.toss(indelrate) and i!=0 and indellist[-1]!=(i-1):
            indellist.append(i)
            if common.toss():
                outstr+=outstr[-1]
                outstr+=readstr[i]
                indelflaglist.append(1)
            else:
                indelflaglist.append(0)
        else:
            if common.toss(errorrate):
                outstr+=switchdict[readstr[i]]
            else:
                outstr+=readstr[i]
        errorrate+=readqualitydeteriorationrate
    indellist.append(len(readstr))
    indelflaglist.append(-1)
    outcigarstr=''
    for i in range(1,len(indellist)):
        subcig=common.cigarsubstr(readcigar,indellist[i-1],indellist[i]-indellist[i-1])
        #print 'here',readcigar,indellist[i-1],indellist[i],subcig
        outcigarstr+=subcig
        if indelflaglist[i]==1:
            outcigarstr+='1I'
        elif indelflaglist[i]==0:
            outcigarstr+='1D'
    
    return [outstr,outcigarstr]

def runfunc(funcname,parameterlist):
    return function_map[funcname](parameterlist)

def createfoldersetup(outputdir,runid,numtypes,numdatasets,configfilename,biasfilename):
    if os.path.exists('%s/%s'%(outputdir,runid)):
        message='%s/%s already exists'%(outputdir,runid)
        common.printstatus(message,'F',common.func_name()) 

    cmd='mkdir -p %s/%s'%(outputdir,runid)
    os.system(cmd)
    cmd='mkdir -p %s/%s/config'%(outputdir,runid)
    os.system(cmd)
    cmd='mkdir -p %s/%s/data'%(outputdir,runid)
    os.system(cmd)
    cmd='mkdir -p %s/%s/metadata'%(outputdir,runid)
    os.system(cmd)
    cmd='cp %s %s/%s/config'%(configfilename,outputdir,runid)
    os.system(cmd)
    cmd='cp %s %s/%s/config'%(biasfilename,outputdir,runid)
    os.system(cmd)
    outdict={}
    outdict['metadata']='%s/%s/metadata'%(outputdir,runid)
    if numtypes==0:
        for ds in range(1,numdatasets+1):
            foldername='%s/%s/data/T%02dS01'%(outputdir,runid,ds)
            cmd='mkdir -p %s'%foldername
            os.system(cmd)
    else:
        for type in range(1,numtypes+1):
            for ds in range(1,numdatasets+1):
                foldername='%s/%s/data/T%02dS%02d'%(outputdir,runid,type,ds)
                cmd='mkdir -p %s'%foldername
                os.system(cmd)
    outdict['data']=[numtypes,numdatasets,'%s/%s/data'%(outputdir,runid)]
    return outdict

def genetxreads2txgenereads(genetxreaddict,genetxcigardict):
    '''
    input:
        genetxreaddict:: gene:genereads,txvector
        genetxcigardict:: gene: tx: chr,strand,start,cigar,size
    '''
    stranddict={0:'-',1:'+'}
    outdict={}
    for gene in genetxreaddict:
        outdict[gene]={}
        genereads,txvect=genetxreaddict[gene]
        txlist=genetxcigardict[gene].keys()
        txsizelist=[genetxcigardict[gene][tx][4] for tx in txlist]
        denom=sum([txvect[i]*txsizelist[i] for i in range(len(txlist))])
        txreadsvec=[genereads*txvect[i]*txsizelist[i]/denom for i in range(len(txlist))]
        for i in range(len(txlist)):
            tx=txlist[i]
            outdict[gene][tx]=[gene,genetxcigardict[gene][tx][0],stranddict[genetxcigardict[gene][tx][1]],tx,
                               txsizelist[i],txvect[i],int(round(txreadsvec[i])),genetxcigardict[gene][tx][2],genetxcigardict[gene][tx][3]]
    return outdict

def genetxdict2chrgenedict(genetxcigardict):
    outdict={}
    for gene in genetxcigardict:
        txlist=genetxcigardict[gene].keys()
        chrname=genetxcigardict[gene][txlist[0]][0]
        outdict[chrname]=outdict.get(chrname,[])+[gene]
    return outdict


def gensimulatedreadsdata(sample_id,replicate_id,folderdict,genetxreaddict,genetxcigardict,config_object):
    if debug_flg==2:
        funcstarttime=time.time()
    fragmentsizerange=[int(config_object.read_extract_parameters[2]),int(config_object.read_extract_parameters[3])]
    txstartbias=config_object.txstartbias
    cutscoredict=dict([(ln.split('\t')[0],float(ln.rstrip('\n').split('\t')[1])) for ln in open(config_object.cutpreferencefile)])
    
    alignsam='%s/T%02dS%02d/alignments.sam'%(folderdict['data'][2],sample_id,replicate_id)
    alignbam='%s/T%02dS%02d/alignments.bam'%(folderdict['data'][2],sample_id,replicate_id)
    alignfout=open(alignsam,'w')
    alignerrsam='%s/T%02dS%02d/alignments_with_errors.sam'%(folderdict['data'][2],sample_id,replicate_id)
    alignerrbam='%s/T%02dS%02d/alignments_with_errors.bam'%(folderdict['data'][2],sample_id,replicate_id)
    alignerrfout=open(alignerrsam,'w')
    txgenereaddict=genetxreads2txgenereads(genetxreaddict,genetxcigardict)
    metadatafile='%s/expression_T%02dS%02d.txt'%(folderdict['metadata'],sample_id,replicate_id)
    fout=open(metadatafile,'w')
    for gene in txgenereaddict:
        for tx in txgenereaddict[gene]:
            fout.write('%s:%s:%s\t%s\tTX\t%d\t%6.4f\t%d\t%d\t%s\n'%tuple(txgenereaddict[gene][tx]))
    fout.close()
    chrgenedict=genetxdict2chrgenedict(genetxcigardict)
    chrlist=chrgenedict.keys()
    chrlist.sort()
    if debug_flg==2:
        tottxstrgettime=0
        totfraglistgettime=0
        totreadlistgettime=0
        toterrorreadgettime=0       
    for chrname in chrlist:
        if not(os.path.isfile('%s/%s.fa'%(config_object.reference_genome_dir,chrname))):
            continue
        chrfafileptr=open('%s/%s.fa'%(config_object.reference_genome_dir,chrname))
        for gene in chrgenedict[chrname]:
            for tx in genetxcigardict[gene]:
                if debug_flg==2:
                    loopstarttime=time.time()
                numreads=txgenereaddict[gene][tx][6]
                startbase,cigarstr=genetxcigardict[gene][tx][2],genetxcigardict[gene][tx][3]
                txstr=getRNAtranscriptstring(chrfafileptr,startbase,cigarstr)
                if debug_flg==2:
                    txstrgettime=time.time()
                    tottxstrgettime+=txstrgettime-loopstarttime
                    message='Gene:%s, Transcript:%s; Transcript Size:%d, Num reads=%d'%(gene,tx,len(txstr),numreads)
                    common.printstatus(message,'S',common.func_name())
                if len(txstr)<fragmentsizerange[0]:
                    continue
                #print 'here',tx,startbase,cigarstr,len(txstr)
                #empirical number of cuts
                numcutmu=max(1,int(len(txstr)*4.0/(fragmentsizerange[0]+fragmentsizerange[1])))
                numcutsig=numcutmu/2.0
                #print gene,tx,len(txstr),fragmentsizerange[0],fragmentsizerange[1],numcutmu,numreads
                fragmentlist=getfragments(txstr,numreads,txgenereaddict[gene][tx][2],cutscoredict,fragmentsizerange,txstartbias,numcutmu,numcutsig)
                if debug_flg==2:
                    fraglistgettime=time.time()
                    totfraglistgettime+=fraglistgettime-txstrgettime
                #print 'here',len(txstr),numcutmu,numreads,len(fragmentlist)
                allparameterlist=[txstr,cigarstr,fragmentlist,config_object.read_extract_parameters]
                readlist=runfunc(config_object.read_extract_method,allparameterlist)
                if debug_flg==2:
                    readlistgettime=time.time()
                    totreadlistgettime+=readlistgettime-fraglistgettime
                inum=0
                for read in readlist:
                    rstartbase,readstr,readcigar=read
                    readstartlocation=sum([int(x) for x in common.cigarsubstr(cigarstr,0,rstartbase+1).replace('M','N').split('N')[:-1]])-1
                    inum+=1
                    alignfout.write('%s:%d:%d\t1\t%s\t%d\t0\t%s\t*\t0\t0\t%s\t*\n'%
                                    (tx,inum,rstartbase,chrname,startbase+readstartlocation,readcigar,readstr))
                    if max(config_object.readqualitydeteriorationrate,config_object.indelrate)>0.0:
                        errreadstr,errcigarstr=adderrortoread(readstr,readcigar,config_object.readqualitydeteriorationrate,config_object.indelrate)
                        alignerrfout.write('%s:%d:%d\t1\t%s\t%d\t0\t%s\t*\t0\t0\t%s\t*\n'%
                                    (tx,inum,rstartbase,chrname,startbase+readstartlocation,errcigarstr,errreadstr))
                if debug_flg==2:
                    errorreadgettime=time.time()
                    toterrorreadgettime+=errorreadgettime-readlistgettime
    alignfout.close()
    alignerrfout.close() 
    if debug_flg==2:
        functime=time.time()-funcstarttime
        othertime=functime-(tottxstrgettime+totfraglistgettime+totreadlistgettime+toterrorreadgettime)
        message='\n\nFunc Time=%s,Txstr=%4.2f,Fragtime=%4.2f,Readtime=%4.2f,Errtime=%4.2f,Oth=%4.2f\n\n'% \
                 (functime,tottxstrgettime*100/functime,totfraglistgettime*100/functime,
                  totreadlistgettime*100/functime,toterrorreadgettime*100/functime,othertime*100/functime)     
        common.printstatus(message,'S',common.func_name())    
    cmd='%s/samtools view -bt %s %s >%s'%(config_object.pathsamtools,config_object.chrfaifile,alignsam,alignbam)
    message='Running %s'%cmd
    common.printstatus(message,'S',common.func_name())
    os.system(cmd)
    if max(config_object.readqualitydeteriorationrate,config_object.indelrate)>0.0:
        cmd='%s/samtools view -bt %s %s >%s'%(config_object.pathsamtools,config_object.chrfaifile,alignerrsam,alignerrbam)
        message='Running %s'%cmd
        common.printstatus(message,'S',common.func_name())
        os.system(cmd)
    return 1

def deletesamfiles(folderdict):
    numtypes,replicates,datadir=folderdict['data']
    samples=max(1,numtypes)
    for i in range(1,samples+1):
        for j in range(1,replicates+1):
            cmd='rm %s/T%02dS%02d/*.sam'%(datadir,i,j)
            os.system(cmd)
    return 1
            

class runparam:
    def __init__(self, config_file):
        if os.path.isfile(config_file):
            self.config_file=config_file
        else:
            message='Config file %s does not exist'%config_file
            common.printstatus(message,'F',common.func_name())     
            
    def _checkdir(self,dirname,dirpath,function_name):
        if not(os.path.exists(dirpath)):
            message='%s does not exist'%dirname
            common.printstatus(message,'F',function_name) 
        else:
            return dirpath
        
    def _checkfile(self,filename,filepath,function_name):
        if not(os.path.isfile(filepath)):
            message='%s does not exist'%filename
            common.printstatus(message,'F',function_name) 
        else:
            return filepath       
        
    def _checknum(self,x,type,function_name):
        if type=='int':
            try:
                int(x)
                return int(x)
            except:
                message='%s is not integer'%x
                common.printstatus(message,'F',function_name) 
        if type=='float':
            try:
                float(x)
                return float(x)
            except:
                message='%s is not float'%x
                common.printstatus(message,'F',function_name) 
    
    def parse(self):
        config = ConfigParser.ConfigParser()
        config.read(self.config_file)
        
        self.pathsamtools=self._checkdir('samtools path',config.get('tools','pathsamtools'),common.func_name())
        
        self.reference_genome_dir=self._checkdir('reference genome dir',config.get('reference','reference_genome_dir'),common.func_name())
        self.chrfaifile =self._checkfile('chromosome index file',config.get('reference','chrfaifile'),common.func_name())
        #self.chromosome_size_file =self._checkfile('chromosome size file',config.get('reference','chromosome_size_file'),common.func_name())
        
        metasourcedict=dict(config.items('metasource'))
        for filename in metasourcedict:
            self._checkfile(filename,metasourcedict[filename],common.func_name())
        self.metasourcedict=metasourcedict
        
        random_seed=self._checknum(config.get('source','random_seed'),'int',common.func_name())
        print random_seed
        if random_seed<0:
            seed=random.randint(0,10000)
        else:
            seed=random_seed
        random.seed(seed)
        
        self.generate_method=config.get('source','generate_method')
        self.generate_method_files=config.get('source','generate_method_files').split(',')
        self.generate_method_parameters=config.get('source','generate_method_parameters').split(',')
        self.generate_method_parameters[2]=int(self.generate_method_parameters[2])
        if len(self.generate_method_parameters)==3:
            self.generate_method_parameters.append([])
        else:
            genefile=self._checkfile('gene list file',self.generate_method_parameters[3],common.func_name())
            ingenefilelist=list(set([line.strip().split('\t')[0].split(':')[0] for line in open(genefile)]))
            if len(ingenefilelist)>self.generate_method_parameters[2]:
                self.generate_method_parameters[3]=random.sample(ingenefilelist,self.generate_method_parameters[2])
            else:
                self.generate_method_parameters[3]=ingenefilelist[0:]
                
            
        self.output_dir=config.get('RNAseqdatasets','output_dir')
        self.run_name=config.get('RNAseqdatasets','run_name')
        self.numtypes=self._checknum(config.get('RNAseqdatasets','numtypes'),'int',common.func_name())
        self.numdatasets=self._checknum(config.get('RNAseqdatasets','numdatasets'),'int',common.func_name())
        self.coverageratiodistribution=[self._checknum(x,'float',common.func_name()) for x in config.get('RNAseqdatasets','coverageratiodistribution').split(',')]
        self.jsddistribution=[self._checknum(x,'float',common.func_name()) for x in config.get('RNAseqdatasets','jsddistribution').split(',')]
        
        self.cutpreferencefile=self._checkfile('cut preference file',config.get('bias','cutpreferencefile'),common.func_name())
        #self.pctcutpoints=self._checknum(config.get('bias','pctcutpoints'),'float',common.func_name())
        self.txstartbias=self._checknum(config.get('bias','txstartbias'),'float',common.func_name())
        
        self.read_extract_method=config.get('datareads','read_extract_method')
        self.read_extract_parameters=config.get('datareads','read_extract_parameters').split(',')
        if self.read_extract_parameters[0]=='PE':
            self.readcount=self._checknum(config.get('datareads','readcount'),'int',common.func_name())/2
        else:
            self.readcount=self._checknum(config.get('datareads','readcount'),'int',common.func_name())
        
        self.readqualitydeteriorationrate=self._checknum(config.get('quality','readqualitydeteriorationrate'),'float',common.func_name())
        self.indelrate=self._checknum(config.get('quality','indelrate'),'float',common.func_name())

def makemetadataplot(folderdict,totreads,readsize):
    '''
    In PE readsize is twice
    '''
    numtypes,replicates,datadir=folderdict['data']
    samples=max(1,numtypes)
    metadatadir=folderdict['metadata']
    genealldatadict={}
    for i in range(1,samples+1):
        exp_file=open('%s/expression_T%02dS01.txt'%(metadatadir,i))
        genedict={}
        for lntxt in exp_file:
            ln=lntxt.rstrip('\n').split('\t')
            gene=ln[0].split(':')[0]
            txsize=int(ln[3])
            txexp=float(ln[4])
            numreads=int(ln[5])
            if gene in genedict:
                genedict[gene][0].append(txsize)
                genedict[gene][1].append(numreads)
                genedict[gene][2].append(txexp)
            else:
                genedict[gene]=[[txsize],[numreads],[txexp]]
        histbins=len(genedict)/5
        genedatadict={}
        for gene in genedict:
            numtx=len(genedict[gene][0])
            entropy=common.shannon_entropy(genedict[gene][2])
            rpkm=1000000.0/totreads*sum([genedict[gene][1][k]*1.0/genedict[gene][0][k] for k in range(numtx)])
            coverage=sum([genedict[gene][1][k]*1.0/genedict[gene][0][k] for k in range(numtx)])*readsize
            genedatadict[gene]=[numtx,entropy,rpkm,coverage]
            genealldatadict[gene]=genealldatadict.get(gene,[])+[coverage]
        plot.plot_histogram([genedatadict[gene][0] for gene in genedatadict],histbins,
                            'transcripts count','# genes','Number of Transcripts','%s/plt_T%02d_transcript_count.jpg'%(metadatadir,i))
        plot.plot_histogram([genedatadict[gene][1] for gene in genedatadict if genedatadict[gene][1]>1],histbins,
                            'entropy','# genes','Entropy Distribution','%s/plt_T%02d_entropy_dist.jpg'%(metadatadir,i))
        plot.plot_histogram([genedatadict[gene][2] for gene in genedatadict],histbins,
                            'RPKM','# genes','RPKM Distribution','%s/plt_T%02d_rpkm_dist.jpg'%(metadatadir,i),1)            
        plot.plot_histogram([genedatadict[gene][3] for gene in genedatadict],histbins,
                            'coverage','# genes','Coverage distribution','%s/plt_T%02d_coverage_dist.jpg'%(metadatadir,i),1)         
    if samples==2:
        jsd_file=open('%s/expression_jsd.txt'%metadatadir)
        for lntxt in jsd_file:
            ln=lntxt.rstrip('\n').split('\t')
            gene=ln[0]; jsd=float(ln[5])
            #print gene, jsd
            genealldatadict[gene].append(jsd)
        
        #print genealldatadict
        todellist=[]
        covlist1=[];covlist2=[];jsdlist=[]; mincovlist=[]
        for gene in genealldatadict:
            #print genealldatadict[gene]
            if min(genealldatadict[gene][0],genealldatadict[gene][1])>0.5 and genealldatadict[gene][2]>0.02:
                covlist1.append(genealldatadict[gene][0])
                covlist2.append(genealldatadict[gene][1])
                mincovlist.append(min(genealldatadict[gene][0],genealldatadict[gene][1]))
                jsdlist.append(genealldatadict[gene][2])
            else:
                todellist.append(gene)
        numofgenestr='%d of %d'%(len(genedatadict)-len(todellist),len(genedatadict))
        
        imagefilename='%s/plt_jsd_coverage1.pdf'%(metadatadir)
        plot.plotscatterwithhistogram(jsdlist,covlist1,'jsd','coverage1','cov1-jsd\n%s'%numofgenestr,imagefilename,markertup=('.',5),logscaleyflg=1,logscalexflg=0)
        imagefilename='%s/plt_jsd_coverage2.pdf'%(metadatadir)
        plot.plotscatterwithhistogram(jsdlist,covlist2,'jsd','coverage2','cov2-jsd\n%s'%numofgenestr,imagefilename,markertup=('.',5),logscaleyflg=1,logscalexflg=0)
        imagefilename='%s/plt_coverage1_coverage2.pdf'%(metadatadir)
        plot.plotscatterwithhistogram(covlist1,covlist2,'coverage1','coverage2','cov1-cov2\n%s'%numofgenestr,imagefilename,markertup=('.',5),logscaleyflg=1,logscalexflg=1)
        imagefilename='%s/plt_jsd_mincoverage.pdf'%(metadatadir)
        plot.plotscatterwithhistogram(jsdlist,mincovlist,'jsd','mincoverage','mincov-jsd\n%s'%numofgenestr,imagefilename,markertup=('.',5),logscaleyflg=1,logscalexflg=0)

if __name__ == "__main__":
    debug_flg=2
    function_map={'RNAmetasource2source':RNAmetasource2source,'fragment2read':fragment2read}
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-c', '--config', dest='cfgfile',default='config/sregen.cfg')
    (options, args) = parser.parse_args()
    config_file=options.cfgfile         
    
    cfg=runparam(config_file)
    cfg.parse()
        
    message='Program Started'
    common.printstatus(message,'S',common.func_name())
    
    
    metasourcefilelist=[cfg.metasourcedict[x] for x in cfg.generate_method_files]
    parameters=cfg.generate_method_parameters
    parameterlist=[metasourcefilelist,parameters]
    genetxcigardict=runfunc(cfg.generate_method,parameterlist) #chr,strand,start,cigar,size
    
    folderdict=createfoldersetup(cfg.output_dir,cfg.run_name,cfg.numtypes,cfg.numdatasets,config_file,cfg.cutpreferencefile)
    
    genetxcountdict=dict([(gene,len(genetxcigardict[gene])) for gene in genetxcigardict])
    
    if cfg.numtypes==0:
        for sample_id in range(1,cfg.numdatasets+1):
            genetxreaddict=getsimulatedoneexpression(genetxcountdict,cfg.readcount)
            replicate_id=1
            gensimulatedreadsdata(sample_id,replicate_id,folderdict,genetxreaddict,genetxcigardict,cfg)
    if cfg.numtypes==1:
        genetxreaddict=getsimulatedoneexpression(genetxcountdict,cfg.readcount)
        for replicate_id in range(1,cfg.numdatasets+1):
            sample_id=1
            gensimulatedreadsdata(sample_id,replicate_id,folderdict,genetxreaddict,genetxcigardict,cfg)
    if cfg.numtypes==2:
        #print cfg.coverageratiodistribution,cfg.jsddistribution
        genetwoexpressions=getsimulatedtwoexpressions(genetxcountdict,cfg.coverageratiodistribution,cfg.jsddistribution,cfg.readcount)
        genetxreaddict1=dict([(gene,[genetwoexpressions[gene][0],genetwoexpressions[gene][2]]) for gene in genetwoexpressions])
        genetxreaddict2=dict([(gene,[genetwoexpressions[gene][1],genetwoexpressions[gene][3]]) for gene in genetwoexpressions])
        jsdmetadatafile='%s/expression_jsd.txt'%(folderdict['metadata'])
        fout=open(jsdmetadatafile,'w')
        for gene in genetwoexpressions:
            fout.write('%s\t% 7d\t% 7d\t%s\t%s\t%6.4f\n'%(gene,genetwoexpressions[gene][0],genetwoexpressions[gene][1],
                                                          common.fl2str(genetwoexpressions[gene][2]),common.fl2str(genetwoexpressions[gene][3]),
                                                          genetwoexpressions[gene][4]))
        fout.close()
        #Todo Plot
        sample_id=1
        for replicate_id in range(1,cfg.numdatasets+1):
            gensimulatedreadsdata(sample_id,replicate_id,folderdict,genetxreaddict1,genetxcigardict,cfg)            
        sample_id=2
        for replicate_id in range(1,cfg.numdatasets+1):
            gensimulatedreadsdata(sample_id,replicate_id,folderdict,genetxreaddict2,genetxcigardict,cfg)                
    
#    if debug_flg==0:
#        deletesamfiles(folderdict)
    #print folderdict
    deletesamfiles(folderdict)
    
    if cfg.read_extract_parameters[0]=='PE':
        makemetadataplot(folderdict,cfg.readcount,int(cfg.read_extract_parameters[1])*2)
    else:
        makemetadataplot(folderdict,cfg.readcount,int(cfg.read_extract_parameters[1]))
    
    message='Done'
    common.printstatus(message,'S',common.func_name())
       
