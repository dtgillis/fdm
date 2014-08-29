import os
import cPickle
import math
import ast
import itertools
import sys
import scipy.cluster.vq
import scipy.stats
import numpy
import random

import common
import bamfile
import bedgraph
import junctionfile
import actgraph
import plot

#############################
# Useful functions

def extractNdominantsplices(flowfile,dominanttslist,N=2):
    """
    outdict[genepos]=[outflowedges,tsdict]
    genepos=gene,pos
    outflowedges=most dominant edges
    tsdict= 
    """
    outdict={}
    lncnt=0
    for lntxt in open(flowfile):
        lncnt+=1
        ln=lntxt.strip('\n').split('\t')
        if lncnt==1:
            hdict=common.findheaderdict(ln,['wt1','wt2','flow'])
            tslist=hdict.keys()
        else:
            genepos=(ln[2],int(ln[1]))
            flowedges=ast.literal_eval(ln[6])
            tsdict={}
            for ts in tslist:
                flowvalues=ast.literal_eval(ln[hdict[ts][2]])
                wt=min(float(ln[hdict[ts][0]]),float(ln[hdict[ts][1]]))
                tsdict[ts]=[flowvalues,wt]
            domflowvaluelist=[tsdict[ts][0] for ts in dominanttslist]
            #print flowvaluelist,flowvaluelist[0]
            edgevaluelist=[[j[i] for j in domflowvaluelist] for i in range(len(domflowvaluelist[0]))]
            edgesumlist=[-1*sum(x) for x in edgevaluelist]
            rankedge=sorted(range(len(edgesumlist)), key=edgesumlist.__getitem__)
            outedges=rankedge[0:N]
            for ts in tsdict:
                flowvalues=tsdict[ts][0]
                tsdict[ts][0]=[flowvalues[i] for i in outedges]
            if len(flowedges)>=len(outedges):
                outflowedges=[flowedges[i] for i in outedges]
            else:
                outflowedges=flowedges
            outdict[genepos]=[outflowedges,tsdict]
    return outdict
                

def distr2pvalue(distrlist,x):
    eps=0.00000001
    distrlist.sort()
    numperm=len(distrlist)
    for i in range(numperm):
        if x<=distrlist[i]+eps:
            break
    return 1-i*1.0/numperm

#############################

class project:
    '''
    Performs all the tasks of the project
    '''
    def __init__(self,config_object,gtfobj):
        '''
        groups: list of list of ts 
        '''
        self.cfg=config_object
        self.tslist=[]
        for list in self.cfg.project_groups:
            self.tslist+=list
        self.tslist.sort()
        self.dir='%s/project/%s'%(self.cfg.root_dir,self.cfg.project_name)
        if not os.path.isdir(self.dir):
            os.system('mkdir -p %s'%self.dir)
        else:
            message='Project directory %s already exists'%self.dir
            common.printstatus(message,'W',common.func_name())
        
        fdmfastdir='%s/project/%s/ffast'%(self.cfg.root_dir,self.cfg.project_name)
        if not os.path.exists(fdmfastdir):
            cmd='mkdir -p %s'%fdmfastdir  
            os.system(cmd)  
            
        ffulldir='%s/project/%s/ffull'%(self.cfg.root_dir,self.cfg.project_name)
        if not os.path.exists(ffulldir):
            cmd='mkdir -p %s'%ffulldir  
            os.system(cmd)   
            
        reportdir='%s/project/%s/report'%(self.cfg.root_dir,self.cfg.project_name)
        if not os.path.exists(reportdir):
            cmd='mkdir -p %s'%reportdir
            os.system(cmd)
            
        imagedir='%s/project/%s/report/image'%(self.cfg.root_dir,self.cfg.project_name)
        if not os.path.exists(imagedir):
            cmd='mkdir -p %s'%imagedir
            os.system(cmd)     
            
        
        
        listparamdir='%s/project/%s/list'%(self.cfg.root_dir,self.cfg.project_name)
        if not os.path.exists(listparamdir):
            os.system('mkdir -p %s'%listparamdir)
            
            statusdir='%s/project/%s/status'%(self.cfg.root_dir,self.cfg.project_name)
            os.system('mkdir -p %s'%statusdir)
            
            chrlist=gtfobj.getchrlist()
            chrlist=[('chr%s'%x).replace('chrchr','chr') for x in chrlist]
            
            if self.cfg.run_type in [1,2]:
                tspairlist=['%s::%s'%(x[0],x[1]) for x in itertools.combinations(self.tslist,2)]
            elif self.cfg.run_type in [3]:
                tspairlist=['%s::%s'%('+'.join(self.cfg.project_groups[0]),'+'.join(self.cfg.project_groups[1]))]
                tspairlist+=['%s::%s'%(x,y) for x,y in zip(self.cfg.project_groups[0],self.cfg.project_groups[1])]
            elif self.cfg.run_type in [4,5]:
                tspairlist=['BLANK']
                            
            f1=open('%s/alllist.txt'%listparamdir,'w')
            for ts in self.tslist:
                f1.write('%s\n'%ts)
            f1.close()
    
            f1=open('%s/chrlist.txt'%listparamdir,'w')
            for chr in chrlist:
                f1.write('%s\n'%chr)
            f1.close()     
            
            f1=open('%s/allpair.txt'%listparamdir,'w')
            for tspair in tspairlist:
                f1.write('%s\n'%tspair)
            f1.close()
                                    
        
        self.geneprojgrsdict={}
        self.geneprojdivdict={}
        self.gtfobj=gtfobj
        
        
            
    def _createunionjuncfile(self):
        bedfilelist=['%s/dataout/%s/%s.bed'%(self.cfg.root_dir,ts,ts) for ts in self.tslist]
        prfdict={'N':0,'A':1,'S':2,'P':3}    
        juncdict={}
        for bedfile in bedfilelist:
            for lntxt in open(bedfile):
                ln=lntxt.rstrip('\n').split('\t')
                junckey='%s:%011d:%011d'%(ln[0].ljust(5,' '),int(ln[1]),int(ln[2]))
                if junckey not in juncdict:
                    juncdict[junckey]=prfdict[ln[3][0:1]]
        fout=open('%s/%s.jun'%(self.dir,self.cfg.project_name),'w')
        junclist=juncdict.keys()
        junclist.sort()
        for junckey in junclist:
            ln=junckey.split(':')
            ln[0]=ln[0].strip(' ')
            ln[1]=int(ln[1])
            ln[2]=int(ln[2])
            fout.write('%s\t%d\t%d\t%d\n'%(ln[0],ln[1],ln[2],juncdict[junckey]))
    
    def _addsplice2graphstruct(self,graphstruct,addsplicelist):
        '''Add a splicelist to graph structure to get new graph structure
        '''
        if len(addsplicelist)>1000:
            message='Novel splices more than 1000'
            common.printstatus(message,'W',common.func_name())
            return graphstruct
        dummyexonsz=5
        exonlist=graphstruct[0][0:] 
        intronlist=graphstruct[1][0:]
        splicelist= graphstruct[2][0:]
        startnodelist=graphstruct[3][0:]
        endnodelist=graphstruct[4][0:]
        novelnodelist=graphstruct[5][0:]
        exonlist.sort()
        intronlist.sort()
        splicelist.sort()
        startnodelist.sort()
        endnodelist.sort()
 
        for splice in addsplicelist:
            if splice[0]>=splice[1]:
                continue
            if splice not in splicelist:
                splicelist.append(splice)
                newexonlist=[]
                if splice[0]<exonlist[0][0]:
                    newexonlist.append([splice[0]-dummyexonsz,splice[0]])
                    intronlist.insert(0,[splice[0]+1,max(splice[0]+1,exonlist[0][0]-1)])      
                    novelnodelist.append(splice[0])          
                for exon in exonlist:
                    if exon[0]<splice[0] and exon[1]>splice[0]:
                        newexonlist.append([exon[0],splice[0]])
                        newexonlist.append([splice[0]+1,exon[1]])
                        novelnodelist.append(splice[0])
                    else:
                        newexonlist.append(exon)
                exonlist=newexonlist[0:]
                newexonlist=[]
                for exon in exonlist:
                    if exon[0]<splice[1] and exon[1]>splice[1]:
                        newexonlist.append([exon[0],splice[1]-1])
                        newexonlist.append([splice[1],exon[1]])
                        novelnodelist.append(splice[1])
                    else:
                        newexonlist.append(exon)
                exonlist=newexonlist[0:]
                if splice[1]>exonlist[-1][1]:
                    exonlist.append([splice[1],splice[1]+dummyexonsz])
                    intronlist.append([exonlist[-2][1]+1,min(splice[1]-1,exonlist[-1][1]+1)])  
                    novelnodelist.append(splice[1])
        for splice in addsplicelist:
            newintronlist=[]
            for intron in intronlist:
                if intron[0]<splice[0] and intron[1]>splice[0]:
                    newintronlist.append([intron[0],splice[0]])
                    newintronlist.append([splice[0]+1,intron[1]])
                    novelnodelist.append(splice[0])
                else:
                    newintronlist.append(intron)
            intronlist=newintronlist[0:]
            newintronlist=[]
            for intron in intronlist:
                if intron[0]<splice[1] and intron[1]>splice[1]:
                    newintronlist.append([intron[0],splice[1]-1])
                    newintronlist.append([splice[1],intron[1]])
                    novelnodelist.append(splice[1])
                else:
                    newintronlist.append(intron)
            intronlist=newintronlist[0:]                             
        splicelist.sort()
        for intron in intronlist:
            if intron[0]>intron[1]:
                intronlist[:]=(value for value in intronlist if value!=intron)
        for exon in exonlist:
            if exon[0]>exon[1]:
                exonlist[:]=(value for value in exonlist if value!=exon)
        return (exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist)
    
    def splitexonwithinsplice(self,graphstruct):
        exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist=graphstruct
        for [splice[0],splice[1]] in splicelist:
            coverexon=[]
            if splice in exonlist+intronlist:
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
                
                
    
    def setprojmetadict(self,collate_flg=2):
        '''
        create proj annotated graphstruct and div dicts
            include non-annotated skipped exons
        create proj all graphstruct and div dicts
        1: create junction union + graphs/divdict (if different)
        2: create junction union + graphs/divdict (always)
        3: graphs/divdict (if different)
        4: graphs/divdict - always
        Note: 3 and 4 assume union juntion exists
        
        '''
        gtfo=self.gtfobj
        if collate_flg in [1,2]:
            self._createunionjuncfile()
        unionjuncfile='%s/%s.jun'%(self.dir,self.cfg.project_name)
        ujunc=junctionfile.junFile(unionjuncfile,[self.cfg.pathucsctools,self.cfg.chromszfile])
        if not ujunc.checkindex():
            ujunc.buildindex()
        geneprojannogrsdict={}
        geneprojallgrsdict={}
        geneprojannodivdict={}
        geneprojalldivdict={}
        islandlist=gtfo.getislandlist()
        for island in islandlist:
            gene=island[3]
            #debug
            #if gene not in ['B3GALT6','FAM132A','BC033949']:
            #    continue
            #debug
            annographstruct=gtfo.getgene2annographstruct(gene)
            key='%s:%d-%d'%(island[0],island[1],island[2])
            allsplices=ujunc.retrieve(key)
            exonstrlist=[exon[0] for exon in annographstruct[0]]
            exonendlist=[exon[1] for exon in annographstruct[0]]
            projannosplices=[]
            projallsplices=[]
            for splice in allsplices:
                if splice[3] in [1,2]:
                    if (splice[2] in exonstrlist) and (splice[1] in exonendlist):
                        projannosplices.append([splice[1],splice[2]])
                projallsplices.append([splice[1],splice[2]])
            
            tprojannographstruct=self._addsplice2graphstruct(annographstruct,projannosplices)
            projannographstruct=common.splitexonwithinsplice(tprojannographstruct)
            tprojallgraphstruct=self._addsplice2graphstruct(annographstruct,projallsplices)
            projallgraphstruct=common.splitexonwithinsplice(tprojallgraphstruct)
            #debug
            #print gene
            #print 'anno'
            #print annographstruct
            #print projannosplices
            #print projannographstruct
            #print 'all'
            #print annographstruct
            #print projallsplices
            #print projallgraphstruct
            #debug
            
            if annographstruct!=projannographstruct or (collate_flg in [2,4]):
                projannodivdict=common.graphstruct2divdict(projannographstruct)
            else:
                projannodivdict=gtfo.getgene2annodivdict(gene)

            if annographstruct!=projallgraphstruct or (collate_flg in [2,4]):
                projalldivdict=common.graphstruct2divdict(projallgraphstruct)
            else:
                projalldivdict=gtfo.getgene2annodivdict(gene)             
            
            geneprojannogrsdict[gene]=projannographstruct
            geneprojallgrsdict[gene]=projallgraphstruct
            geneprojannodivdict[gene]=projannodivdict
            geneprojalldivdict[gene]=projalldivdict
            
        cPickle.dump(geneprojannogrsdict,open('%s/geneprojannogrsdict.pck'%self.dir,'w'))  
        cPickle.dump(geneprojallgrsdict,open('%s/geneprojallgrsdict.pck'%self.dir,'w'))  
        cPickle.dump(geneprojannodivdict,open('%s/geneprojannodivdict.pck'%self.dir,'w'))  
        cPickle.dump(geneprojalldivdict,open('%s/geneprojalldivdict.pck'%self.dir,'w'))  
    
    def getgeneprojgrs(self,gene):
        self._setgeneprojdict()
        return self.geneprojgrsdict[gene]
    
    def getgeneprojdiv(self,gene):
        self._setgeneprojdict()
        return self.geneprojdivdict[gene]    
            
    def _setgeneprojdict(self):
        '''
        type    = 1: annotated + skipped exons
                = 2: All novels included
        '''
        #self._setprojmetadict()
        if self.geneprojgrsdict=={}:
            self.geneprojgrsdict=cPickle.load(open('%s/geneprojallgrsdict.pck'%self.dir))
            #graph always complete
        if self.geneprojdivdict=={}:
            if self.cfg.project_type==1:
                self.geneprojdivdict=cPickle.load(open('%s/geneprojannodivdict.pck'%self.dir))
            else:
                self.geneprojdivdict=cPickle.load(open('%s/geneprojalldivdict.pck'%self.dir))            
         
    def _bdglinesexonlist2wt(self,bdglines,exonlist,meanflg=1):
        '''
        mean_flg
        =1 for mean
        =0 for median
        median not yet implemented
        '''
        invlist=[]; invwtlist=[]
        for bdgline in bdglines:
            invlist.append([bdgline[1],bdgline[2]-1])
            invwtlist.append(bdgline[3])
        exoncovlist=common.exonbaseoverlap(exonlist,invlist,invwtlist)
        exonwtlist=[exoncovlist[i]/(1.0*(max(exonlist[i][1]-exonlist[i][0],0)+1)) for i in range(len(exonlist))]
        return exonwtlist

    def _junlinessplicelist2wt(self,junlines,splicelist):
        '''
        '''
        junsplicelist=[]; junwtlist=[]
        for junline in junlines:
            junsplicelist.append([junline[1],junline[2]])
            junwtlist.append(junline[3])
        splicewtlist=[]
        for splice in splicelist:
            if splice in junsplicelist:
                splicewtlist.append(junwtlist[junsplicelist.index(splice)])
            else:
                splicewtlist.append(0)
        return splicewtlist
            
    def makeAct(self, ts, islandlist=['ALL']):
        '''type= 1 : Annotated and skipped exons
               = 2 : Include Novel
        '''
        maxexons=2000
        if ts not in self.tslist:
            message='Sample %s is not part of the project %s'%(ts,self.cfg.project_name)
            common.printstatus(message,'F',common.func_name())
            return

        actfilename='%s/dataout/%s/%s_%s.act'%(self.cfg.root_dir,ts,self.cfg.project_name,ts)
        wgsfilename='%s/dataout/%s/%s_%s.wgs'%(self.cfg.root_dir,ts,self.cfg.project_name,ts)
        rpkfilename='%s/dataout/%s/%s_%s.rpk'%(self.cfg.root_dir,ts,self.cfg.project_name,ts)
       
        bdgfilename='%s/dataout/%s/%s.bdg'%(self.cfg.root_dir,ts,ts)
        junfilename='%s/dataout/%s/%s.jun'%(self.cfg.root_dir,ts,ts)
        
        bam0=bamfile.bamFile(self.cfg.data_dict[ts],[self.cfg.pathsamtools,self.cfg.pathbedtools])
        nummapped,readsize=bam0.getstats()
        
        if islandlist==['ALL']:
            islandlist=self.gtfo.getislandlist()
        
        wgsdict={} 
        fout=open(actfilename,'a')
        frpk=open(rpkfilename,'w')
        bdgfile=bedgraph.bdgFile(bdgfilename,[self.cfg.pathucsctools,self.cfg.chromszfile])
        junfile=junctionfile.junFile(junfilename,[self.cfg.pathucsctools,self.cfg.chromszfile])
        #islandlist.sort()
        for island in islandlist:
            gene=island[3]
            
            key='%s:%d-%d'%(tuple(island[0:3]))
            grs=self.getgeneprojgrs(gene)       
            #print gene, grs 
            exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist=grs
            for intron in intronlist:
                if intron[1]<intron[0]:
                    intronlist[:]=(value for value in intronlist if value!=intron)
            if len(exonlist)>maxexons or len(intronlist)>maxexons:
                message='Gene %s is rejected. Number of exons(introns) is too large %d(%d)'%(gene,len(exonlist),len(intronlist))
                common.printstatus(message,'E',common.func_name())
                continue
            bdglines=bdgfile.retrieve(key)     
            junlines=junfile.retrieve(key)     
            exonwtlist=self._bdglinesexonlist2wt(bdglines,exonlist)      
            intronwtlist=self._bdglinesexonlist2wt(bdglines,intronlist)      
            splicewtlist=self._junlinessplicelist2wt(junlines,splicelist)        
            wgrstuple=(exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist)
            wgsdict[gene]=wgrstuple
            wtgrs=actgraph.wtgraphstruct(wgrstuple)      
            actfilelines=wtgrs.Toactfilelines(island)
            for actfileline in actfilelines:
                fout.write('%s\n'%actfileline)    
            totalcov=0.0; totalexon=0
            for i in range(len(exonlist)):
                exonsz=exonlist[i][1]-exonlist[i][0]+1
                totalexon+=exonsz
                totalcov+=exonwtlist[i]*exonsz
            avgcov=totalcov/totalexon
            rpkm=avgcov*1000*1000000/(nummapped*readsize)
            frpk.write('%s\tchr%s\t%d\t%d\t%10.4f\t%10.4f\n'%(island[3],island[0],island[1],island[2],avgcov,rpkm))
            #message='MakeAct gene done %s'%gene
            #common.printstatus(message,'W',common.func_name())  
        cPickle.dump(wgsdict,open(wgsfilename,'w'))
        return actfilename

                    
    def extractflows(self,flowprefix,islandlist):
        '''
         rootdir,chrlist,islanddict,genedivdict,outfilename)
        geneflowdict[gene]=position:[[outgoingflg,exonstartflg],[[wt1,wt2],nflowvec]]
        genedivdict[gene]=position:[[incoming/outgoing=0,1;exonstart=0=no,1=yes,2=start transcript,3=end transcript][exon, flowlist coordinates]]
            exonstart=0=no,1=yes,2=start transcript and exonstart,3=end transcript and exonstart,4=start transcript and no exonstart(insplice),
        5=end transcript and no exonstart outsplice 
        '''
        tempfiledir='%s/flows/temp.%s'%(self.dir,flowprefix)
        os.system('mkdir -p %s'%tempfiledir)
        outfilename='%s/flows/%s_flows.txt'%(self.dir,flowprefix)
        fout=open(outfilename,'w')
        fout.write('chr\tpos\tgene\toutflag\texonstartflag\tdivexon\tflowedges')
        for ts in self.tslist:
            fout.write('\t%s:wt1\t%s:wt2\t%s:flow'%(ts,ts,ts))
        fout.write('\n') 
        chrlist=self.gtfobj.getchrlist()
        islanddict=self.gtfobj.getislanddict()
        #print chrlist
        #print islandlist[0]
        #genelist=genedivdict.keys()
        #tempfix
        chrlist=[('chr%s'%x).replace('chrchr','chr') for x in chrlist]
        for chr in chrlist:
            chrislanddict={}
            chrgenedivdict={}     
            for island in islandlist:
                #print island
                gene=island[3]
                if island[0]==chr:
                    chrgenedivdict[gene]=self.getgeneprojdiv(gene)
                    chrislanddict[gene]=islanddict[gene]
            if len(chrislanddict.keys())==0:
                continue
            chrgenelist=chrgenedivdict.keys()
            chrgenelist.sort()
            geneflowdictlist=[]
            ftemplist=[]
            for ts in self.tslist:
                actfilename='%s/%s_%s.act'%(self.cfg.datadir_dict[ts],self.cfg.project_name,ts)
                message='Extracting Flows for %s: %s'%(chr,actfilename)
                common.printstatus(message,'S',common.func_name())
                ftempname='%s.temp.%s'%(outfilename,ts)
                ftempout=open(ftempname,'w')
                ftemplist.append(ftempname)
                act0=actgraph.actFile(actfilename)
                geneflowdict=act0.genedivdictTogeneflowdict(chrgenedivdict)
                for gene in geneflowdict:
                    divdict=chrgenedivdict[gene]
                    poslist=divdict.keys()
                    poslist.sort()
                    for pos in poslist:
                        if pos in geneflowdict[gene]:
                            ftempout.write('%s\t%d\t%s\t%d\t%d\t%s\t%s'%   \
                                       (chr,pos,gene,divdict[pos][0][0],divdict[pos][0][1],str(divdict[pos][1][0]),str(divdict[pos][1][1])))            
                            flowposlist=geneflowdict[gene][pos]
                            ftempout.write('\t%10.4f\t%10.4f\t%s\n'%(flowposlist[1][0][0],flowposlist[1][0][1],common.fl2str(flowposlist[1][1]))) 
                ftempout.close()
            common.mergefilelist(ftemplist,range(7,10),fout,range(7))
        
        for tfilename in ftemplist:
            os.system('mv %s %s'%(tfilename,tempfiledir))
    
    def concatchrflows(self):
        chrlist=self.gtfobj.getchrlist()
        chrlist=[('chr%s'%x).replace('chrchr','chr') for x in chrlist]
        wildcardfilelist=['%s/flows/%s_%s_flows.txt'%(self.dir,self.cfg.flow_prefix,chr.upper()) for chr in chrlist]
        flowfilelist=[f for f in wildcardfilelist if os.path.isfile(f)]
        outfilename='%s/flows/%s_ALL_flows.txt'%(self.dir,self.cfg.flow_prefix)
        common.concatfiles(outfilename,flowfilelist,headerflg=1)
        tempfiledir='%s/flows/temp'%(self.dir)
        os.system('mkdir -p %s'%tempfiledir)
        for tfilename in flowfilelist:
            os.system('mv %s %s'%(tfilename,tempfiledir))
        
    def deletetempflows(self):
        folder='%s/flows'%(self.dir)
        tempdirs = [ d for d in os.listdir(folder) if not os.path.isfile('%s/%s'%(folder,d)) and d.startswith('temp') ]   
        for dir in tempdirs:
            os.system('rm -r %s/%s'%(folder,dir))
    
        
    def _flowveclist2flowvec(self,flowveclist):
        #print flowveclist
        flowvec=[]
        for i in range(len(flowveclist[0])):
            flowsum=0
            for j in range(len(flowveclist)):
                flowsum+=flowveclist[j][i]
            flowvec.append(flowsum/len(flowveclist))
        return flowvec
    
    def _flowveclistlist2flowveclist(self,flowveclistlist):
        '''
        ith vector of all lists is averaged to give ith vector 
        [[[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]], [[0.5, 0.5], [0.2, 0.8], [0.0, 1.0]]]
        result:
        [[0.75, 0.25], [0.6, 0.4], [0.5, 0.5]]
        '''
        flowvectlist=[]
        for i in range(len(flowveclistlist[0])):
            tempflowveclist=[]
            for j in range(len(flowveclistlist)):
                tempflowveclist.append(flowveclistlist[j][i])
            flowvectlist.append(self._flowveclist2flowvec(tempflowveclist))
        return flowvectlist
        
    
    def _flowpair2fdm(self,flow1,flow2):
        #print flow1, flow2
        fdm=0
        if len(flow1)!=len(flow2):
            message='Flow sizes are different %s:%s'%(common.fl2str(flow1),common.fl2str(flow2))
            common.printstatus(message,'E',common.func_name())   
        else:
            for i in range(len(flow1)):
                fdm+=1.0/2*math.fabs(flow1[i]-flow2[i]) 
        return fdm
        
    def _flowveclists2fdmpval(self,infdm,flowveclist1,flowveclist2,numpart,numperm):
        #print 'xln1820',flowveclist1,flowveclist2
        if numpart==1:
            pvalue=-1.0
            return (infdm,pvalue)
        else:
            epsilon=0.000000001
            fdmnulllist=[]
            for i in range(numperm):
                allflowvec=flowveclist1+flowveclist2
                flowgrp=common.randgroup(allflowvec,2)
                flowgrp1=self._flowveclist2flowvec(flowgrp[0])
                flowgrp2=self._flowveclist2flowvec(flowgrp[1])
                fdmnulllist.append(self._flowpair2fdm(flowgrp1,flowgrp2))
                
            pvalue=distr2pvalue(fdmnulllist,infdm) 
            return (infdm,pvalue)      
         
    def compute_fast_fdm(self,flowprefix,ffastprefix,tspair):
        #ffast: no p-values   
        flowfile='%s/project/%s/flows/%s_ALL_flows.txt'%(self.cfg.root_dir,self.cfg.project_name,flowprefix)
        self.ffast_dict={}
        lncnt=0
        for lntxt in open(flowfile):
            lncnt+=1
            ln=lntxt.rstrip('\n').split('\t')
            if lncnt==1:
                hdict=common.findheaderdict(ln)
                tslist=hdict.keys()
                tslist.sort()
                ffastfile='%s/project/%s/ffast/%s_ALL_fdm.txt'%(self.cfg.root_dir,self.cfg.project_name,ffastprefix)
                ffastpck='%s/project/%s/ffast/%s_ALL_fdm.pck'%(self.cfg.root_dir,self.cfg.project_name,ffastprefix)
                fout=open(ffastfile,'w')
                fout.write('\t'.join(ln[0:7]))
                tsi='_'.join(tspair[0])
                tsj='_'.join(tspair[1])
                if tsi<tsj:
                    fout.write('\t%s__%s:minwt\t%s__%s:fdm'%(tsi,tsj,tsi,tsj))
                else:
                    fout.write('\t%s__%s:minwt\t%s__%s:fdm'%(tsj,tsi,tsj,tsi))
                fout.write('\n')
            else:
                fout.write('\t'.join(ln[0:7]))
                gene=ln[2]; pos=int(ln[1])
                if gene not in self.ffast_dict:
                    self.ffast_dict[gene]={}
                if pos not in self.ffast_dict[gene]:
                    self.ffast_dict[gene][pos]={}
                wtlist=[]
                flowlist1=[]; flowlist2=[]
                for ts in tspair[0]:
                    wtlist+=[float(ln[hdict[ts][0]]),float(ln[hdict[ts][1]])]
                    flow=common.str2fl(ln[hdict[ts][2]])   
                    weightedflow=[x*wtlist[-1] for x in flow]
                    flowlist1.append(weightedflow)
                for ts in tspair[1]:
                    wtlist+=[float(ln[hdict[ts][0]]),float(ln[hdict[ts][1]])]
                    flow=common.str2fl(ln[hdict[ts][2]])   
                    weightedflow=[x*wtlist[-1] for x in flow]
                    flowlist2.append(weightedflow)     
                tsi='_'.join(tspair[0])
                tsj='_'.join(tspair[1])
                if tsi<tsj:
                    tskey='%s__%s'%(tsi,tsj)
                else:
                    tskey='%s__%s'%(tsj,tsi)                         

                minwt=min(wtlist)   
                flow1=common.normalize_vector([sum(x) for x in zip(*flowlist1)])
                flow2=common.normalize_vector([sum(x) for x in zip(*flowlist2)])
                fdm=self._flowpair2fdm(flow1,flow2) 
                self.ffast_dict[gene][pos][tskey]=[minwt,fdm]
                fout.write('\t%8.2f\t%6.4f'%(minwt,fdm))
                fout.write('\n') 
        cPickle.dump(self.ffast_dict,open(ffastpck,'w'))  
        
    def merge_fast_fdm(self):
        folder='%s/project/%s/ffast'%(self.cfg.root_dir,self.cfg.project_name)
        prefix=self.cfg.project_name
        dirfiles = ['%s/%s'%(folder,f) for f in os.listdir(folder) if os.path.isfile('%s/%s'%(folder,f)) and f.startswith(prefix) and f.split('.')[1]=='txt' ]
        dirfiles.sort()
        fout=open('%s/project/%s/ffast/%s_ALL_fdm.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name),'w')
        common.mergefilelist(dirfiles,range(7,9),fout,range(7))
        fout.close()
        
        dictfiles = ['%s/%s'%(folder,f) for f in os.listdir(folder) if os.path.isfile('%s/%s'%(folder,f)) and f.startswith(prefix) and f.split('.')[1]=='pck' ]
        dictlist=[cPickle.load(open(f)) for f in dictfiles]
        fout=open('%s/project/%s/ffast/%s_ALL_fdm.pck'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name),'w')
        mergeddict=common.mergedictlist(dictlist)
        cPickle.dump(mergeddict,fout)
        fout.close()
        

    def filter_fast_fdm(self):
        ffastall='%s/project/%s/ffast/%s_ALL_fdm.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        fout=open('%s/project/%s/ffast/%s_FIL_fdm.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name),'w')
        lncnt=0
        for lntxt in open(ffastall):
            lncnt+=1
            ln=lntxt.rstrip('\n').split('\t')
            if lncnt==1:
                fout.write(lntxt)
                hdict=common.findheaderdict(ln,['minwt','fdm'])
                tskeylist=hdict.keys()
                if self.cfg.run_type==3:
                    project_group_keys=['_'.join(y) for y in self.cfg.project_groups]
                    project_group_keys.sort()
                    fdmtskeylist=['__'.join(x) for x in itertools.combinations(project_group_keys,2)]
                else:
                    fdmtskeylist=hdict.keys()
                #print fdmtskeylist
            else:
                minwt =min([float(ln[hdict[tskey][0]]) for tskey in tskeylist])
                maxfdm=max([float(ln[hdict[tskey][1]]) for tskey in fdmtskeylist])
                if minwt>self.cfg.ffast_min_cov and maxfdm>self.cfg.ffast_min_fdm:
                    fout.write(lntxt)       
                    
    def splitgenelist(self):
        fdmfil='%s/project/%s/ffast/%s_FIL_fdm.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        if os.path.exists(fdmfil):
            genelist=list(set([lntxt.split('\t')[2] for lntxt in open(fdmfil) if lntxt.split('\t')[2]!='gene']))
        else:
            genelist=[]
        splitgenedir='%s/project/%s/list/splitgene'%(self.cfg.root_dir,self.cfg.project_name)
        os.system('mkdir -p %s'%splitgenedir)
        numfiles=int(len(genelist)/self.cfg.ffull_genesplit_size)+1
        for i in range(numfiles):
            gfile=open('%s/G%04d'%(splitgenedir,i+1),'w')
            sgenelist=genelist[i*self.cfg.ffull_genesplit_size:(i+1)*self.cfg.ffull_genesplit_size]
            for gene in sgenelist:
                gfile.write('%s\n'%gene)
            gfile.close()

    def extractFEF(self):
        flowfile='%s/project/%s/flows/%s_ALL_flows.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        FEFfile1='%s/project/%s/flows/%s_ALL_FEFs.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        FEFfile2='%s/project/%s/flows/%s_ALL_FEF_processed.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        fout1=open(FEFfile1,'w')
        fout2=open(FEFfile2,'w')
        lncnt=0
        for lntxt in open(flowfile):
            lncnt+=1
            ln=lntxt.rstrip('\n').split('\t')
            if lncnt==1:
                hdict=common.findheaderdict(ln)
                tslist=hdict.keys()
                tslist.sort()
                fout1.write('\t'.join(ln[0:7]))
                fout1.write('\tedge used')
                for ts in tslist:
                    fout1.write('\t%s:FEF\t%s:wt'%(ts,ts))
                fout1.write('\n')
                fout2.write('\t'.join(ln[0:7]))
                fout2.write('\tedge used')
                for ts in tslist:
                    fout2.write('\t%s:FEF'%(ts))
                fout2.write('\n')                
            else:
                edges=ast.literal_eval(ln[6])
                numrows=max(1,len(edges)-1)
                for i in range(numrows):
                    fout1.write('\t'.join(ln[0:7]))
                    fout1.write('\t%s'%str(edges[i]))
                    fout2.write('\t'.join(ln[0:7]))
                    fout2.write('\t%s'%str(edges[i]))                    
                    for ts in tslist:
                        wt1=float(ln[hdict[ts][0]])
                        wt2=float(ln[hdict[ts][1]])
                        flow=ast.literal_eval(ln[hdict[ts][2]])
                        FEF=flow[i]
                        fout1.write('\t%6.4f\t%6.4f'%(FEF,wt2))
                        if wt2<self.cfg.ffast_min_cov:
                            FEF=-1
                        fout2.write('\t%6.4f'%(FEF))
                    fout1.write('\n')
                    fout2.write('\n')         
        fout1.close()
        fout2.close()
           
    def getFEFlinecluster(self,headerdict,FEFline,numclusters):
        tslist=headerdict.keys()
        tslist.sort()
        FEFlist=[]
        goodtslist=[]
        missinglist=[]
        #print FEFline
        #print headerdict
        for ts in tslist:
            if float(FEFline[headerdict[ts][1]])>self.cfg.ffast_min_cov:
                FEFlist.append(float(FEFline[headerdict[ts][0]]))
                goodtslist.append(ts)
            else:
                missinglist.append(ts)
        FEFlist=numpy.array(FEFlist)
        #print 'F',FEFlist
        centroids,distortion=scipy.cluster.vq.kmeans(FEFlist,numclusters)
        indexes,distortionlist=scipy.cluster.vq.vq(FEFlist,centroids)
        clusterlist=[]; mediantslist=[]; FEFclusterlist=[]
        #print indexes
        for i in range(numclusters):
            cluster=[]
            for j in range(len(indexes)):
                if indexes[j]==i:
                    cluster.append(goodtslist[j])
            clusterlist.append(cluster)
            FEFcluster=FEFlist[indexes==i]
            FEFclusterlist.append(FEFcluster.tolist())
            try:
                min_d,idx=min((val, id) for (id, val) in enumerate(distortionlist) if indexes[id]==i)
            except:
                idx=goodtslist.index(goodtslist[-1])
            mediantslist.append(goodtslist[idx])
        DBidx,_,mclist=common.DaviesBouldinIndex(FEFclusterlist)
        #print centroids,mclist
        return [DBidx,FEFline[0:8],clusterlist,missinglist,FEFclusterlist,centroids,mediantslist]
            
    def getallFEFclusterreport(self):
        '''
        Add 3 cluster
        '''
        reportdir='%s/project/%s/report'%(self.cfg.root_dir,self.cfg.project_name)
        clusterdir='%s/project/%s/report/cluster'%(self.cfg.root_dir,self.cfg.project_name)
        cmd='mkdir -p %s'%clusterdir
        os.system(cmd)
        minsampvar=0.1
        mingooddatapct=0.6
        minclusterszpct=0.2
        FEFfile1='%s/project/%s/flows/%s_ALL_FEFs.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        worklist2=[]
        worklist3=[]
        lncnt=0
        for lntxt in open(FEFfile1):
            lncnt+=1
            ln=lntxt.rstrip('\n').split('\t')
            if lncnt==1:
                hdict=common.findheaderdict(ln,coltype=['FEF','wt'])
                tslist=hdict.keys()
                tslist.sort()
                csvhlist2=ln[0:8]+['% Good Data','Not Enough Data','Sample STD','Cluster2 DBI','Cluster2 image','Cluster2.1','Cluster2.2',
                                  'Cluster2.1 ACT','Cluster2.2 ACT']
                csvhline2='\t'.join(csvhlist2)
                csvhlist3=ln[0:8]+['% Good Data','Not Enough Data','Sample STD','Cluster3 DBI','Cluster3 image','Cluster3.1','Cluster3.2','Cluster3.3',
                                  'Cluster3.1 ACT','Cluster3.2 ACT','Cluster3.3 ACT']
                csvhline3='\t'.join(csvhlist3)                            
            else:
                allFEFdata=[float(ln[hdict[ts][0]]) for ts in tslist if float(ln[hdict[ts][1]])>self.cfg.ffast_min_cov]
                percentagegooddata=1.0*len(allFEFdata)/len(tslist)
                samplevariance=scipy.stats.nanstd(allFEFdata)
                if samplevariance>minsampvar and percentagegooddata>mingooddatapct:
                    #print ln
                    #print 'AF',allFEFdata
                    FEFdatadict=dict([(ts,[float(ln[hdict[ts][0]]),float(ln[hdict[ts][1]])]) for ts in tslist])
                    FEFclustertuple=self.getFEFlinecluster(hdict,ln,2)
                    clusterlist=FEFclustertuple[2]
                    minclussize=min([len(x) for x in clusterlist])
                    if 1.0*minclussize/len(FEFdatadict)>minclusterszpct:
                        worklist2.append(FEFclustertuple+[FEFdatadict,samplevariance])
                    FEFclustertuple=self.getFEFlinecluster(hdict,ln,3)
                    clusterlist=FEFclustertuple[2]
                    minclussize=min([len(x) for x in clusterlist])
                    if 1.0*minclussize/len(FEFdatadict)>minclusterszpct:                    
                        worklist3.append(FEFclustertuple+[FEFdatadict,samplevariance])
        worklist2.sort()
        toplist=worklist2[0:self.cfg.report_top_x]
        csvfile='%s/top_cluster2_report.txt'%reportdir
        csvptr=open(csvfile,'w')
        csvptr.write('#\t%s\n'%csvhline2)
        lncnt=0
        for line in toplist:
            lncnt+=1
            filenameprefix='%s_%s__%010d_%010d'%(line[1][2],line[1][0],ast.literal_eval(line[1][7])[0],ast.literal_eval(line[1][7])[0])
            fclust21=open('%s/%s_2_1.txt'%(clusterdir,filenameprefix),'w')
            for ts in line[2][0]:
                fclust21.write('%s\n'%ts)
            fclust21.close()
            fclust22=open('%s/%s_2_2.txt'%(clusterdir,filenameprefix),'w')
            for ts in line[2][1]:
                fclust22.write('%s\n'%ts) 
            fclust22.close()
            fmiss=open('%s/%s_miss.txt'%(clusterdir,filenameprefix),'w')
            for ts in line[3]:
                fmiss.write('%s\n'%ts)     
            fmiss.close()
            fimage2='%s/%s_cluster2.pdf'%(clusterdir,filenameprefix)  
            FEFdatadict=line[7]  
            genechredge=[line[1][2],line[1][0],ast.literal_eval(line[1][7])]
            if lncnt<=self.cfg.graph_top_x:
                plot.plotFEFclusters(FEFdatadict,line[2],line[6],genechredge,fimage2)   
            l2=ast.literal_eval(line[1][6])
            l2.append(ast.literal_eval(line[1][5]))
            minx=min([min(x) for x in l2])
            maxx=max([max(x) for x in l2])
            urlimg='<a href="cluster/%s_cluster2.pdf">cluster</a>'%(filenameprefix)
            urlactlist=['<a href="image/%s_%s__%s__%d-%d.pdf">actg</a>'%(self.cfg.project_name,ts,line[1][2],minx,maxx) for ts in line[6]]
            urlmissing='<a href="cluster/%s_miss.txt">missing</a>'%(filenameprefix)
            urlclust1='<a href="cluster/%s_2_1.txt">cluster 1</a>'%(filenameprefix)
            urlclust2='<a href="cluster/%s_2_2.txt">cluster 2</a>'%(filenameprefix)
            csvline=line[1][0:8]+['%6.2f'%(100*(1-1.0*len(line[3])/len(FEFdatadict))),urlmissing,'%6.4f'%line[8],'%6.4f'%line[0]] 
            csvline.append(urlimg)
            csvline.append(urlclust1)
            csvline.append(urlclust2)
            csvline+=urlactlist
            csvptr.write('%d\t%s\n'%(lncnt,'\t'.join(csvline)))
            genepos=[line[1][2],[minx,maxx],int(line[1][1])] 
            if lncnt<=self.cfg.graph_top_x:
                self.createactgraphs(line[6],[genepos])
        csvptr.close()
        htmlfile=open('%s/top_cluster2_report.html'%reportdir,'w')    
        common.csv2html(open(csvfile),htmlfile)
        htmlfile.close()
        
        worklist3.sort()
        toplist=worklist3[0:self.cfg.report_top_x]
        csvfile='%s/top_cluster3_report.txt'%reportdir
        csvptr=open(csvfile,'w')
        csvptr.write('#\t%s\n'%csvhline3)
        lncnt=0
        for line in toplist:
            lncnt+=1
            filenameprefix='%s_%s__%010d_%010d'%(line[1][2],line[1][0],ast.literal_eval(line[1][7])[0],ast.literal_eval(line[1][7])[0])
            fclust31=open('%s/%s_3_1.txt'%(clusterdir,filenameprefix),'w')
            for ts in line[2][0]:
                fclust31.write('%s\n'%ts)
            fclust31.close()
            fclust32=open('%s/%s_3_2.txt'%(clusterdir,filenameprefix),'w')
            for ts in line[2][1]:
                fclust32.write('%s\n'%ts) 
            fclust32.close()
            fclust33=open('%s/%s_3_3.txt'%(clusterdir,filenameprefix),'w')
            for ts in line[2][2]:
                fclust33.write('%s\n'%ts)                 
            fclust33.close()
            fmiss=open('%s/%s_miss.txt'%(clusterdir,filenameprefix),'w')
            for ts in line[3]:
                fmiss.write('%s\n'%ts)     
            fmiss.close()
            fimage3='%s/%s_cluster3.pdf'%(clusterdir,filenameprefix)  
            FEFdatadict=line[7]  
            genechredge=[line[1][2],line[1][0],ast.literal_eval(line[1][7])]
            if lncnt<=self.cfg.graph_top_x:
                plot.plotFEFclusters(FEFdatadict,line[2],line[6],genechredge,fimage3)   
            l2=ast.literal_eval(line[1][6])
            l2.append(ast.literal_eval(line[1][5]))
            minx=min([min(x) for x in l2])
            maxx=max([max(x) for x in l2])
            urlimg='<a href="cluster/%s_cluster3.pdf">cluster</a>'%(filenameprefix)
            urlactlist=['<a href="image/%s_%s__%s__%d-%d.pdf">actg %s</a>'%(self.cfg.project_name,ts,line[1][2],minx,maxx,ts) for ts in line[6]]
            urlmissing='<a href="cluster/%s_miss.txt">missing</a>'%(filenameprefix)
            urlclust1='<a href="cluster/%s_3_1.txt">cluster 1</a>'%(filenameprefix)
            urlclust2='<a href="cluster/%s_3_2.txt">cluster 2</a>'%(filenameprefix)
            urlclust3='<a href="cluster/%s_3_3.txt">cluster 3</a>'%(filenameprefix)            
            csvline=line[1][0:8]+['%6.2f'%(100*(1-1.0*len(line[3])/len(FEFdatadict))),urlmissing,'%6.4f'%line[8],'%6.4f'%line[0]] 
            csvline.append(urlimg)
            csvline.append(urlclust1)
            csvline.append(urlclust2)
            csvline.append(urlclust3)            
            csvline+=urlactlist
            csvptr.write('%d\t%s\n'%(lncnt,'\t'.join(csvline)))
            genepos=[line[1][2],[minx,maxx],int(line[1][1])] 
            if lncnt<=self.cfg.graph_top_x:
                self.createactgraphs(line[6],[genepos])
        csvptr.close()
        htmlfile=open('%s/top_cluster3_report.html'%reportdir,'w')    
        common.csv2html(open(csvfile),htmlfile)
        htmlfile.close()



        
    def gettwolargegroupreport(self,twogroupdifflist):
        """
        Make boxplot, html
        """
        reportdir='%s/project/%s/report'%(self.cfg.root_dir,self.cfg.project_name)
        csvhlist=['chr','pos','gene','outflag','exonstartflag','divexon','flowedges','edge used']
        csvhlist+=['Size Group1','Size Group2','Group1 Median','Group2 Median','Random Median1','Random Median2','Sample STD']
        csvhlist+=['Sample Range','Kruskal pvalue','Median Random pvalue','Group1 Median ACT','Group2 Median ACT','Boxplot']
        csvhline='\t'.join(csvhlist)
        worklist=[]
        for line in twogroupdifflist:
            worklist.append([[float(line[16]),1-float(line[15])],line])
        worklist.sort()
        toplist=worklist[0:self.cfg.report_top_x]
        csvfile='%s/top_twolargegroup_report.txt'%reportdir
        csvptr=open(csvfile,'w')
        csvptr.write('#\t%s\n'%csvhline)
        lncnt=0
        for line in toplist:
            lncnt+=1
            l2=ast.literal_eval(line[1][6])
            l2.append(ast.literal_eval(line[1][5]))
            minx=min([min(x) for x in l2])
            maxx=max([max(x) for x in l2])
            #print line
            #print line[1][18],line[1][0],line[1][2],minx,maxx
            urlact1='<a href="image/%s_%s__%s__%d-%d.pdf">actg %s</a>'%(self.cfg.project_name,line[1][18],line[1][2],minx,maxx,line[1][18])
            urlact2='<a href="image/%s_%s__%s__%d-%d.pdf">actg %s</a>'%(self.cfg.project_name,line[1][19],line[1][2],minx,maxx,line[1][19])
            imagefilename='%s/image/%s_jitter__%s__%d-%d.pdf'%(reportdir,self.cfg.project_name,line[1][2],minx,maxx)
            ylabel='Fractional Edge Weight'
            xlabel='Group1 vs Group2'
            title='%s__%d-%d\np-value=%s'%(line[1][2],minx,maxx,line[1][16])
            if lncnt<=self.cfg.graph_top_x:
                plot.makejitter(line[1][20],line[1][21],[],ylabel,xlabel,title,imagefilename)
            urlimg='<a href="image/%s_jitter__%s__%d-%d.pdf">boxplot</a>'%(self.cfg.project_name,line[1][2],minx,maxx)
            csvline=line[1][0:18] 
            csvline.append(urlact1)
            csvline.append(urlact2)
            csvline.append(urlimg)
            csvptr.write('%d\t%s\n'%(lncnt,'\t'.join(csvline)))
            genepos=[line[1][2],[minx,maxx],int(line[1][1])]
            mediantslist=[line[1][18],line[1][19]]
            #print mediantslist
            if lncnt<=self.cfg.graph_top_x:
                self.createactgraphs(mediantslist,[genepos])
        csvptr.close()
        htmlfile=open('%s/top_twolargegroup_report.html'%reportdir,'w')    
        common.csv2html(open(csvfile),htmlfile)
        htmlfile.close()
                

    def twolargegroupdifference(self):
        outlinelist=[]
        mingroupsz=10
        numiter=101
        minstd=0.3
        FEFfile2='%s/project/%s/flows/%s_ALL_FEF_processed.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        twogroupdiff='%s/project/%s/report/%s_Two_Group_diff.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)        
        fout=open(twogroupdiff,'w')
        lncnt=0
        for lntxt in open(FEFfile2):
            lncnt+=1
            ln=lntxt.rstrip('\n').split('\t')
            if lncnt==1:
                hdict=common.findheaderdict(ln,coltype=['FEF'])
                tslist=hdict.keys()
                tslist.sort()
                fout.write('\t'.join(ln[0:8]))
                fout.write('Size Group1\tSize Group2\tGroup1 Median\tGroup2 Median\tRandom Median1\tRandom Media2\tSample Variance')
                fout.write('\tSample Range\tKruskal pvalue\tMedian Random pvalue\tGroup1 Median ACT\tGroup2 Median ACT\n')
            else:
                allgroup1=[float(ln[hdict[ts][0]]) for ts in tslist if ts in self.cfg.project_groups[0]]
                allgroup2=[float(ln[hdict[ts][0]]) for ts in tslist if ts in self.cfg.project_groups[1]]
                group1=[x for x in allgroup1 if x>=0.0]
                group2=[x for x in allgroup2 if x>=0.0]
                group1.sort()
                group2.sort()
                if len(group1)<mingroupsz or len(group2)<mingroupsz:
                    continue 
                medgroup1=common.median(group1)
                medgroup2=common.median(group2)
                medts1=[ts for ts in tslist if ts in self.cfg.project_groups[0] and float(ln[hdict[ts][0]])-medgroup1<0.001 ][0]
                medts2=[ts for ts in tslist if ts in self.cfg.project_groups[1] and float(ln[hdict[ts][0]])-medgroup2<0.001 ][0]
                bothgroup=group1+group2
                samplestd=scipy.stats.nanstd(bothgroup)
                if samplestd>minstd:
                    #print group1, group2
                    _,grouppval=scipy.stats.kruskal(numpy.array(group1),numpy.array(group2))
                    rndpvaluelist=[]
                    for j in range(numiter):
                        #print j
                        random.shuffle(bothgroup)
                        rndgrp1=bothgroup[0:len(group1)]
                        rndgrp2=bothgroup[len(group1):]
                        _,rndpval=scipy.stats.kruskal(numpy.array(rndgrp1),numpy.array(rndgrp2))
                        rndpvaluelist.append(rndpval)
                    
                    outline=ln[0:8]+['%d'%len(group1),'%d'%len(group2),'%5.4f'%medgroup1,'%5.4f'%medgroup2,'%5.4f'%common.median(rndgrp1),'%5.4f'%common.median(rndgrp2),
                             '%5.4f'%samplestd,'%5.4f'%(max(bothgroup)-min(bothgroup)),'%5.4f'%grouppval,'%5.4f'%common.median(rndpvaluelist),medts1,medts2,group1,group2]
                    outlinelist.append(outline)
                    fout.write('%s\n'%'\t'.join(outline[0:-2]))        
        fout.close()     
        return outlinelist
        
        
    def compute_full_fdm(self,tspair,genelist):
        '''
        _flowveclist2flowvec
        genefdmdict = gene:pos: [(fdm,pvalue)] for pool,[(fdm,pvalue)1st,(fdm,pvalue)2nd,(fdm,pvalue)12],[posflgs],[wtpairlist1, 2],[[flowvec allparts] list1, 2]
        '''
        self._setgeneprojdict() 
        
        #print self.cfg.flow_prefix
        #print self.cfg.ffast_prefix
        #print self.cfg.ffull_prefix
        tspair[0].sort()
        tspair[1].sort()
        fdmset1='_'.join(tspair[0])
        fdmset2='_'.join(tspair[1])
        if fdmset2<fdmset1:
            fdmset1,fdmset2=fdmset2,fdmset1
            
        ffulldir='%s/project/%s/ffull'%(self.cfg.root_dir,self.cfg.project_name)
        
        fdmrawfilename='%s/%s__%s__%s__raw.txt'%(ffulldir,self.cfg.ffull_prefix,fdmset1,fdmset2)
        fdmsmalloutfilename='%s/%s__%s__%s__summary.txt'%(ffulldir,self.cfg.ffull_prefix,fdmset1,fdmset2)
        if os.path.exists(fdmrawfilename):
            os.system('rm %s'%fdmrawfilename)
        
        ffastpck='%s/project/%s/ffast/%s_ALL_fdm.pck'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
        num_partition=self.cfg.ffull_partition
        num_permutation=self.cfg.ffull_permutation
        
        self.ffast_dict=cPickle.load(open(ffastpck))
        
        islanddict=self.gtfobj.getislanddict()
        useislanddict={}
        usegenedivdict={}
        for gene in genelist:
            useislanddict[gene]=islanddict[gene]
            usegenedivdict[gene]=self.geneprojdivdict[gene]
        
        #genelist=self.genedivdict.keys()
        geneposflgdictdict={}
        genewtpairdict1={}
        genewtpairdict2={}
        geneflowvecdict1={}
        geneflowvecdict2={}
        genefdmdict={}
        for gene in genelist:
            geneposflgdictdict[gene]={}
            genewtpairdict1[gene]={}
            geneflowvecdict1[gene]={}
            genewtpairdict2[gene]={}
            geneflowvecdict2[gene]={}  
            genefdmdict[gene]={}          
        

        #ffull: Working with p-values
        for i in range(len(tspair[0])):
            bamfilename=self.cfg.data_dict[tspair[0][i]]
            #inconsistent tooldirs
            bamfile1=bamfile.bamFile(bamfilename,[self.cfg.pathsamtools,self.cfg.pathbedtools])
            #print len(usegenedivdict), len(useislanddict), num_permutation
            geneflowlistdict=bamfile1.Togeneflowlistdict(usegenedivdict, useislanddict,num_partition)
            #print geneflowlistdict.keys()
            for gene in genelist:
                flowdictlist=geneflowlistdict[gene]
                #All flowdict are identifical in keys
                poslist=flowdictlist[0].keys()
                for pos in poslist:
                    wtpairlist=[]
                    flowveclist=[]
                    for flowdict in flowdictlist:
                        if i==0:
                            geneposflgdictdict[gene][pos]=flowdict[pos][0]
                            genewtpairdict1[gene][pos]=[]
                            geneflowvecdict1[gene][pos]=[]
                        wtpairlist.append(flowdict[pos][1][0])
                        flowveclist.append(flowdict[pos][1][1])
                    genewtpairdict1[gene][pos].append(wtpairlist)
                    geneflowvecdict1[gene][pos].append(flowveclist)
        
        for i in range(len(tspair[1])):
            bamfilename=self.cfg.data_dict[tspair[1][i]]
            bamfile2=bamfile.bamFile(bamfilename,[self.cfg.pathsamtools,self.cfg.pathbedtools])
            geneflowlistdict=bamfile2.Togeneflowlistdict(usegenedivdict, useislanddict,num_partition)
            for gene in genelist:
                flowdictlist=geneflowlistdict[gene]
                #All flowdict are identifical in keys
                poslist=flowdictlist[0].keys()
                for pos in poslist:
                    wtpairlist=[]
                    flowveclist=[]
                    for flowdict in flowdictlist:
                        if i==0:
                            geneposflgdictdict[gene][pos]=flowdict[pos][0]
                            genewtpairdict2[gene][pos]=[]
                            geneflowvecdict2[gene][pos]=[]
                        wtpairlist.append(flowdict[pos][1][0])
                        flowveclist.append(flowdict[pos][1][1])
                    genewtpairdict2[gene][pos].append(wtpairlist)
                    geneflowvecdict2[gene][pos].append(flowveclist)
                             

        for gene in genelist:
            for pos in geneflowvecdict1[gene].keys():
                flowveclist1=self._flowveclistlist2flowveclist(geneflowvecdict1[gene][pos])
                flowveclist2=self._flowveclistlist2flowveclist(geneflowvecdict2[gene][pos])
                tskey='%s__%s'%(fdmset1,fdmset2)
                infdm=self.ffast_dict[gene][pos][tskey][1]
                fdmpval=self._flowveclists2fdmpval(infdm,flowveclist1,flowveclist2,num_partition,num_permutation)
                genefdmdict[gene][pos]=[fdmpval,[],geneposflgdictdict[gene][pos], \
                                    [genewtpairdict1[gene][pos],genewtpairdict2[gene][pos]], \
                                    [flowveclist1,flowveclist2]] 

        self.genefdmdict=genefdmdict
        
        self.printgenefdm(fdmrawfilename)
        
        fsmallout=open(fdmsmalloutfilename,'w')
        fdmpair='%s__%s'%(fdmset1,fdmset2)
        fsmallout.write('gene\tpos\t%s:fdm\t%s:pvalue\t%s:significance\n'%(fdmpair,fdmpair,fdmpair))
        lncnt=0
        for lntxt in open(fdmrawfilename):
            lncnt+=1
            if lncnt==1:
                continue
            ln=lntxt.rstrip('\n').split('\t')
            if float(ln[3])<self.cfg.ffull_pvalue:
                sig_flg='SIG'
            else:
                sig_flg='NOT'
            fsmallout.write('%s\t%s\t%s\t%s\t%s\n'%(ln[0],ln[1],ln[2],ln[3],sig_flg))         
 
        
    def printgenefdm(self,filename):
        fout=open(filename,'w')
        fout.write('Gene\tPosition\tFDM\tp-value\toutflag\texonflg\twtlist1\twtlist2\t\tflowlist1\tflowlist2\n')
        genelist=self.genefdmdict.keys()
        genelist.sort()   
        for gene in genelist:
            fdmdict=self.genefdmdict[gene]
            poslist=fdmdict.keys()
            poslist.sort()
            for pos in poslist:
                fdm,pvalue=fdmdict[pos][0]
                fdmpairlist=fdmdict[pos][1]
                outflg,exonflg=fdmdict[pos][2]
                wtlist1,wtlist2=fdmdict[pos][3]
                flowlist1,flowlist2=fdmdict[pos][4]
                fout.write('%s\t%10d\t%6.4f\t%10.6f\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n'% 
                           (gene,pos,fdm,pvalue,outflg,exonflg,common.fl2str(wtlist1),common.fl2str(wtlist2),common.fl2str(flowlist1),
                            common.fl2str(flowlist2),common.fl2str(fdmpairlist)))
        fout.close()
                
    def ToFDMbedgraph(self,filename,islanddict):
        # to deprecate
        foutfdm=open(filename[:-4]+'__fdm'+filename[-4:],'a')
        foutpval=open(filename[:-4]+'__pvl'+filename[-4:],'a')
        genelist=self.genefdmdict.keys()
        genelist.sort()
        for gene in genelist:
            chrnm='%s'%islanddict[gene][0]
            fdmdict=self.genefdmdict[gene]
            poslist=fdmdict.keys()
            poslist.sort()
            for pos in poslist:
                fdm,pvalue=fdmdict[pos][0]
                outflg,exonflg=fdmdict[pos][2]
                if outflg==0:
                    foutfdm.write('%s\t%d\t%d\t%6.4f\n'%(chrnm,pos,pos+1,fdm))
                    foutpval.write('%s\t%d\t%d\t%6.4f\n'%(chrnm,pos,pos+1,pvalue))
                else:
                    foutfdm.write('%s\t%d\t%d\t%6.4f\n'%(chrnm,pos,pos+1,fdm))
                    foutpval.write('%s\t%d\t%d\t%6.4f\n'%(chrnm,pos,pos+1,pvalue))

    def createactgraphs(self,tslist,geneposlist):
        actfilelist=['%s/dataout/%s/%s_%s.act'%(self.cfg.root_dir,ts,self.cfg.project_name,ts) for ts in tslist]
        imagedir='%s/project/%s/report/image'%(self.cfg.root_dir,self.cfg.project_name)
        for actfilename in actfilelist:
            act0=actgraph.actFile(actfilename)
            for genepos in geneposlist:
                act0.Toimage(genepos[0],imagedir,genepos[1],highlightnodelist=[genepos[2]])     
    
    def mergeffullfiles(self):
        #merge all fdm full common to all run types
        ffulldir='%s/project/%s/ffull'%(self.cfg.root_dir,self.cfg.project_name)
        ffullfiles = [f for f in os.listdir(ffulldir) if f.endswith('summary.txt')]
        ffullkeylist=list(set(['%s__%s'%(f.split('__')[1],f.split('__')[2]) for f in os.listdir(ffulldir) if f.endswith('summary.txt')]))
        ffullkeylist.sort()
        ffullfiledict={}  
        for fkey in ffullkeylist:
            for f in ffullfiles:
                if fkey=='%s__%s'%(f.split('__')[1],f.split('__')[2]):
                    ffullfiledict[fkey]=ffullfiledict.get(fkey,[])+['%s/%s'%(ffulldir,f)]
        allgenefilelist=[]
        for fkey in ffullkeylist:
            outfile='%s/%s_ALL__%s__ffull.txt'%(ffulldir,self.cfg.project_name,fkey)
            allgenefilelist.append(outfile)
            #print ffullfiledict[fkey]
            ffullfiledict[fkey].sort()
            #print ffullfiledict[fkey]
            common.concatfiles(outfile,ffullfiledict[fkey])
        ffullmergefile='%s/%s_MERGED_ffull.txt'%(ffulldir,self.cfg.project_name)
        #print allgenefilelist
        allgenefilelist.sort()
        if not os.path.isfile(ffullmergefile):
            fout=open(ffullmergefile,'w')
            common.mergefilelist(allgenefilelist,range(2,5),fout,range(2))
            fout.close()
        #all additional columns
        ffullallfile='%s/%s_ALL_ffull.txt'%(ffulldir,self.cfg.project_name)
        if not os.path.isfile(ffullallfile):
            flowfile='%s/project/%s/flows/%s_ALL_flows.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
            flowdict={}
            lncnt=0
            for lntxt in open(flowfile):
                lncnt+=1
                ln=lntxt.rstrip('\n').split('\t')
                if lncnt==1:
                    hdict=common.findheaderdict(ln)
                    tslist=hdict.keys()
                    tslist.sort()
                    allhline='%s'%'\t'.join(ln[0:7])
                    for ts in tslist:
                        allhline+='\t%s:wt1\t%s:wt2'%(ts,ts)
                else:
                    genepos=(ln[2],int(ln[1]))
                    outline=ln[0:7]
                    for ts in tslist:
                        outline+=[ln[hdict[ts][0]],ln[hdict[ts][1]]]
                    if genepos not in flowdict:
                        flowdict[genepos]=outline
            ffullout=open(ffullallfile,'w')
            lncnt=0
            for lntxt in open(ffullmergefile):
                lncnt+=1
                ln=lntxt.rstrip('\n').split('\t')
                if lncnt==1:
                    hdict=common.findheaderdict(ln,['fdm','pvalue','significance'])
                    tstslist=hdict.keys()
                    tstslist.sort()
                    for tsts in tstslist:
                        allhline+='\t%s:fdm\t%s:pvalue\t%s:significance'%(tsts,tsts,tsts)                
                    ffullout.write('%s\n'%allhline)
                else:
                    genepos=(ln[0],int(ln[1]))
                    outline=ln[0:7]
                    fdmlist=[]
                    for tsts in tstslist:
                        fdmlist+=[ln[hdict[tsts][0]],ln[hdict[tsts][1]],ln[hdict[tsts][2]]]
                    outline=flowdict[genepos]+fdmlist
                    ffullout.write('%s\n'%'\t'.join(outline))
            ffullout.close()
            # remove unmerged files
            
            os.system('mkdir %s/temp'%(ffulldir))
            for f in [x for x in os.listdir(ffulldir) if x!='%s_ALL_ffull.txt'%self.cfg.project_name]:
                #print f
                os.system('mv %s/%s %s/temp'%(ffulldir,f,ffulldir))

            del_flg=0
            if del_flg:
                os.system('rm -r %s/temp'%(ffulldir))
   
                 
    def run_report(self):
        
        ffulldir='%s/project/%s/ffull'%(self.cfg.root_dir,self.cfg.project_name)
        ffullallfile='%s/%s_ALL_ffull.txt'%(ffulldir,self.cfg.project_name)
        if not os.path.isfile(ffullallfile):
            self.mergeffullfiles()
           
        # run_type==1: N bam files
        if self.cfg.run_type==1:
            #concat gene files    
            reportdir='%s/project/%s/report'%(self.cfg.root_dir,self.cfg.project_name)
            ffulltopfile='%s/%s_TOP_ffull.txt'%(reportdir,self.cfg.project_name)
            worklist=[]
            lncnt=0
            for lntxt in open(ffullallfile):
                lncnt+=1
                ln=lntxt.strip('\n').split('\t')
                if lncnt==1:
                    hdict=common.findheaderdict(ln,['fdm','significance'])
                    hkeys=hdict.keys()
                    csvhlist=ln+['%s:actg'%ts for ts in self.tslist]
                    csvhline='\t'.join(csvhlist)
                else:
                    l2=ast.literal_eval(ln[6])
                    l2.append(ast.literal_eval(ln[5]))
                    minx=min([min(x) for x in l2])
                    maxx=max([max(x) for x in l2])
                    fdmvalues=[float(ln[hdict[h][0]]) for h in hkeys]
                    sigvalues=[ln[hdict[h][1]] for h in hkeys]   
                    if 'SIG' in sigvalues: 
                        urllist=['<a href="image/%s_%s__%s__%d-%d.pdf">actg</a>'%(self.cfg.project_name,ts,ln[2],minx,maxx) for ts in self.tslist]
                        urllist.sort()
                        worklist.append([max(fdmvalues),[ln[2],[minx,maxx],int(ln[1])],ln+urllist])
            worklist.sort()
            worklist.reverse()
            toplist=worklist[0:self.cfg.report_top_x]
            #print worklist
            #print toplist
            csvfile='%s/top_report.txt'%reportdir
            csvptr=open(csvfile,'w')
            csvptr.write('#\t%s\n'%csvhline)
            lncnt=0
            for line in toplist:
                lncnt+=1
                csvptr.write('%d\t%s\n'%(lncnt,'\t'.join(line[-1])))
            csvptr.close()
            htmlfile=open('%s/top_report.html'%reportdir,'w')    
            common.csv2html(open(csvfile),htmlfile)
            htmlfile.close()
            
            #make act image files
            geneposlist=[x[-2] for x in toplist]
            self.createactgraphs(self.tslist,geneposlist[0:self.cfg.graph_top_x])
            
        # run_type==2: Two set of bam files
        if self.cfg.run_type==2:
            reportdir='%s/project/%s/report'%(self.cfg.root_dir,self.cfg.project_name)
            ffulltopfile='%s/%s_TOP_ffull.txt'%(reportdir,self.cfg.project_name)
            sigdict={'SIG':1,'NOT':0}
            worklist=[]
            betweenpair=[]
            for x,y in itertools.product(self.cfg.project_groups[0],self.cfg.project_groups[1]):
                if x<y:
                    betweenpair.append('%s__%s'%(x,y))
                else:
                    betweenpair.append('%s__%s'%(y,x))
            withinpair=[]
            for x,y in itertools.combinations(self.cfg.project_groups[0],2):
                if x<y:
                    withinpair.append('%s__%s'%(x,y))
                else:
                    withinpair.append('%s__%s'%(y,x))    
            for x,y in itertools.combinations(self.cfg.project_groups[1],2):
                if x<y:
                    withinpair.append('%s__%s'%(x,y))
                else:
                    withinpair.append('%s__%s'%(y,x))                        
            lncnt=0
            for lntxt in open(ffullallfile):
                lncnt+=1
                ln=lntxt.strip('\n').split('\t')
                if lncnt==1:
                    hdict=common.findheaderdict(ln,['fdm','significance'])
                    hkeys=hdict.keys()
                    csvhlist=ln+['%s:actg'%ts for ts in self.tslist]+['Between','Within','Between-Within','Group Diff p-value']
                    csvhline='\t'.join(csvhlist)
                else:
                    l2=ast.literal_eval(ln[6])
                    l2.append(ast.literal_eval(ln[5]))
                    minx=min([min(x) for x in l2])
                    maxx=max([max(x) for x in l2])
                    betsigvalues=[sigdict[ln[hdict[h][1]]] for h in hkeys if h in betweenpair]  
                    withinsigvalues=[sigdict[ln[hdict[h][1]]] for h in hkeys if h in withinpair]  
                    ranker=sum(betsigvalues)-sum(withinsigvalues)
                    groupdiffpvalue=common.getgroupshufflepvalue(betsigvalues,withinsigvalues,10000)
                    if (groupdiffpvalue < 0.05) or (0 not in betsigvalues and 1 not in withinsigvalues): 
                        urllist=['<a href="image/%s_%s__%s__%d-%d.pdf">actg</a>'%(self.cfg.project_name,ts,ln[2],minx,maxx) for ts in self.tslist]
                        urllist.sort()
                        betweengroupcnt=sum(betsigvalues)
                        withingroupcnt=sum(withinsigvalues)
                        worklist.append([ranker,[ln[2],[minx,maxx],int(ln[1])],ln+urllist+['%2d'%betweengroupcnt,'%2d'%withingroupcnt,'%2d'%ranker,'%5.4f'%groupdiffpvalue]])
            worklist.sort()
            worklist.reverse()
            toplist=worklist[0:self.cfg.report_top_x]
            #print worklist
            #print toplist
            csvfile='%s/top_report.txt'%reportdir
            csvptr=open(csvfile,'w')
            csvptr.write('#\t%s\n'%csvhline)
            lncnt=0
            for line in toplist:
                lncnt+=1
                csvptr.write('%d\t%s\n'%(lncnt,'\t'.join(line[-1])))
            csvptr.close()
            htmlfile=open('%s/top_report.html'%reportdir,'w')    
            common.csv2html(open(csvfile),htmlfile)
            htmlfile.close()
            
            #make act image files
            geneposlist=[x[-2] for x in toplist]
            self.createactgraphs(self.tslist,geneposlist[0:self.cfg.graph_top_x])
            
        # run_type==3: paired groups
        if self.cfg.run_type==3:
            #concat gene files
            pvaluethr=0.05
            reportdir='%s/project/%s/report'%(self.cfg.root_dir,self.cfg.project_name)
            ffulltopfile='%s/%s_TOP_ffull.txt'%(reportdir,self.cfg.project_name)
            numpairs=len(self.cfg.project_groups[0])
            sigtranslatedict={'SIG':1,'NOT':0}
            flowfile='%s/project/%s/flows/%s_ALL_flows.txt'%(self.cfg.root_dir,self.cfg.project_name,self.cfg.project_name)
            splicetupdict=extractNdominantsplices(flowfile,self.cfg.project_groups[0],2)
            worklist=[]
            lncnt=0
            for lntxt in open(ffullallfile):
                lncnt+=1
                ln=lntxt.strip('\n').split('\t')
                if lncnt==1:
                    hdict=common.findheaderdict(ln,['fdm','significance','pvalue'])
                    hkeys=hdict.keys()
                    hkeys.sort()
                    csvhlist=ln+['%s:actg'%ts for ts in self.tslist]+['arrow plot1','arrow plot2','Dominant Edge1','Dominant edge2']+['Fisher combined p-value','Liptak combined p-value','Fisher Bonferroni Significant','Liptak Bonferroni Significant']
                    csvhline='\t'.join(csvhlist)
                    poolkey='%s__%s'%('_'.join(self.cfg.project_groups[0]),'_'.join(self.cfg.project_groups[1]))
                else:
                    l2=ast.literal_eval(ln[6])
                    l2.append(ast.literal_eval(ln[5]))
                    minx=min([min(x) for x in l2])
                    maxx=max([max(x) for x in l2])
                    fdmvalues=[float(ln[hdict[h][0]]) for h in hkeys if h!=poolkey]
                    fdmpool=float(ln[hdict[poolkey][0]])
                    sigvalues=[sigtranslatedict[ln[hdict[h][1]]] for h in hkeys if h!=poolkey]
                    sigpool=ln[hdict[poolkey][1]]
                    pvaluelist=[float(ln[hdict[h][2]]) for h in hkeys if h!=poolkey]   
                    #print fdmvalues,fdmpool, sigvalues, sigpool,pvaluelist
                    if sigpool=='SIG' and fdmpool>self.cfg.ffast_min_fdm: 
                        fcombp=common.combinepvalues(pvaluelist,method='fisher')
                        lcombp=common.combinepvalues(pvaluelist,method='liptak')
                        urllist=['<a href="image/%s_%s__%s__%d-%d.pdf">actg</a>'%(self.cfg.project_name,ts,ln[2],minx,maxx) for ts in self.tslist]
                        urllist.append('<a href="image/%s_arrow__%s__%s_1.pdf">arrowplot1</a>'%(self.cfg.project_name,ln[2],ln[1]))
                        urllist.append('<a href="image/%s_arrow__%s__%s_2.pdf">arrowplot2</a>'%(self.cfg.project_name,ln[2],ln[1]))
                        arrowfilename1='%s/image/%s_arrow__%s__%s_1.pdf'%(reportdir,self.cfg.project_name,ln[2],ln[1])
                        arrowfilename2='%s/image/%s_arrow__%s__%s_2.pdf'%(reportdir,self.cfg.project_name,ln[2],ln[1])
                        
                        genepos=(ln[2],int(ln[1]))
                        tsinhkeys=self.cfg.project_groups[0]+self.cfg.project_groups[1]
                        outflowedges,tsflowdict=splicetupdict[genepos]
                        edge1=tuple(outflowedges[0])
                        edge2=tuple(outflowedges[-1])
                        edge1data=[tsflowdict[ts][0][0] for ts in tsinhkeys]
                        edge2data=[tsflowdict[ts][0][1] for ts in tsinhkeys]
                        wtlist=[tsflowdict[ts][1] for ts in tsinhkeys]
                        wtflglist=[]
                        for i in range(numpairs):
                            if min(wtlist[i],wtlist[i+numpairs])<self.cfg.ffast_min_cov:
                                wtflglist.append(0)
                            else:
                                wtflglist.append(1)
                        edge1outdata={edge1:zip(edge1data[:numpairs],edge1data[numpairs:])}
                        edge2outdata={edge2:zip(edge2data[:numpairs],edge2data[numpairs:])}
                        totalmovement1=0.0;totalmovement2=0.0
                        for i in range(numpairs):
                            if wtflglist[i]==1:
                                totalmovement1+=(edge1data[numpairs+i]-edge1data[i])
                                totalmovement2+=(edge2data[numpairs+i]-edge2data[i])
                        totalmovement=abs(totalmovement1)+abs(totalmovement2)
                        if totalmovement>=4.0:
                            movement=0
                        elif totalmovement>=2.0 and totalmovement<4.0:
                            movement=1
                        elif totalmovement>=1.0 and totalmovement<2.0:
                            movement=2
                        else:
                            movement=3
                        worklist.append([[movement,math.log10(fcombp)+math.log10(lcombp),fcombp,lcombp,fdmpool],[arrowfilename1,arrowfilename2,sigvalues],[ln[2],[minx,maxx],int(ln[1])],ln+urllist+['%.3e'%fcombp,'%.3e'%lcombp]])
            benfthr=pvaluethr/len(worklist)
            
            #outdict[genepos]=[outflowedges,tsdict]
            
            for line in worklist:
                sigwt=0; sigtup=[]
                if line[0][1]<benfthr:
                    sigwt+=1;sigtup.append('F BSIG')
                else:
                    sigtup.append('F BNOT')
                if line[0][2]<benfthr:
                    sigwt+=1;sigtup.append('L BSIG')
                else:
                    sigtup.append('L BNOT')
                line[0].insert(0,-1*sigwt)
                line[-1]+=sigtup                
            worklist.sort()
            #worklist.reverse()
            toplist=worklist[0:self.cfg.report_top_x]
            csvfile='%s/top_report.txt'%reportdir
            csvptr=open(csvfile,'w')
            csvptr.write('#\t%s\n'%csvhline)
            lncnt=0
            for line in toplist:
                lncnt+=1
                genepos=(line[2][0],line[2][2])
                outline=line[-1]
                tsinhkeys=self.cfg.project_groups[0]+self.cfg.project_groups[1]
                outflowedges,tsflowdict=splicetupdict[genepos]
                edge1=tuple(outflowedges[0])
                edge2=tuple(outflowedges[-1])
                edge1data=[tsflowdict[ts][0][0] for ts in tsinhkeys]
                edge2data=[tsflowdict[ts][0][1] for ts in tsinhkeys]
                wtlist=[tsflowdict[ts][1] for ts in tsinhkeys]
                wtflglist=[]
                for i in range(numpairs):
                    if min(wtlist[i],wtlist[i+numpairs])<self.cfg.ffast_min_cov:
                        wtflglist.append(0)
                    else:
                        wtflglist.append(1)
                edge1outdata={edge1:zip(edge1data[:numpairs],edge1data[numpairs:])}
                edge2outdata={edge2:zip(edge2data[:numpairs],edge2data[numpairs:])}
                outline.insert(-4,str(edge1outdata))
                outline.insert(-4,str(edge2outdata))
                csvptr.write('%d\t%s\n'%(lncnt,'\t'.join(outline)))
                sigflglist=line[1][2]
                xlabellist=['P%s'%i for i in range(1,len(self.cfg.project_groups[0])+1)]
                ylabel='Fraction Splice weight'
                title='%s:%d-%d'%(line[2][0],edge1[0],edge1[1])
                imagefilename=line[1][0]
                y1y2data=zip(edge1data[:numpairs],edge1data[numpairs:])
                if lncnt<=self.cfg.graph_top_x:
                    plot.arrowplot(title,xlabellist,ylabel,y1y2data,sigflglist,wtflglist,imagefilename)
                title='%s:%d-%d'%(line[2][0],edge2[0],edge2[1])
                imagefilename=line[1][1]
                y1y2data=zip(edge2data[:numpairs],edge2data[numpairs:])
                if lncnt<=self.cfg.graph_top_x:
                    plot.arrowplot(title,xlabellist,ylabel,y1y2data,sigflglist,wtflglist,imagefilename)
            csvptr.close()
            htmlfile=open('%s/top_report.html'%reportdir,'w')    
            common.csv2html(open(csvfile),htmlfile)
            htmlfile.close()
            
            #make act image files
            geneposlist=[x[-2] for x in toplist]
            self.createactgraphs(self.tslist,geneposlist[0:self.cfg.graph_top_x])
        if self.cfg.run_type==4:
            self.extractFEF()
            self.getallFEFclusterreport()
        if self.cfg.run_type==5:
            twogroupdifflist=self.twolargegroupdifference()
            self.gettwolargegroupreport(twogroupdifflist)
            
            
            
