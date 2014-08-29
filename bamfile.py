import os
import common
import actgraph

class bamFile:
    '''Class to index and retrieve data from large files
    input is file and chr:start-end and gene for act file
    1) bam : uses samtools to build index and retrieve data
    2) bdg : bedgraph file chr start end wt
    3) jun : junction file chr start end wt
    4) act : actg file : gff format: only indexed by genes
    samtools and bedtools paths are known
    '''
    def __init__(self, filename,tooldirs):
        self.name=filename
        self.ext=filename[-3:]
        self.idxname=filename+'.'+ self.ext[:-1]+'i'
        self.tooldirs=tooldirs
    
    def checkindex(self):
        if os.path.isfile(self.idxname):
            return 1
        return 0
    
    def buildindex(self):
        os.system('mv %s %s.unsorted'%(self.name,self.name))
        os.system('%s/samtools sort %s.unsorted %s'%(self.tooldirs[0],self.name,self.name[:-4]))
        os.system('%s/samtools index %s'%(self.tooldirs[0],self.name))

    def retrieve(self,key):
        outlines=[]
        cmd='%s/samtools view %s %s'%(self.tooldirs[0],self.name,key)
        stdo=os.popen(cmd)
        for linetxt in stdo:
            line=linetxt.rstrip('\n').split('\t')
            outlines.append([line[2],line[3],line[5]])
        return outlines
      
    def getchrlist(self):
        cmd='%s/samtools idxstats %s'%(self.tooldirs[0],self.name)
        chrlist=[]
        stdo=os.popen(cmd)
        for linetxt in stdo:
            line=linetxt.rstrip('\n').split('\t')
            chrlist.append(line[0])
        chrlist.sort()
        if '*' in chrlist:
            chrlist.remove('*')
        return chrlist 
    
    def getstats(self):
        cmd='%s/samtools idxstats %s '%(self.tooldirs[0],self.name)
        stdo=os.popen(cmd)
        nummapped=0; minmappedcnt=10000000000
        for linetxt in stdo:
            line=linetxt.rstrip('\n').split('\t')
            nummappedchr=int(line[2])
            nummapped+=nummappedchr
            if minmappedcnt>nummappedchr and nummappedchr>0:
                minmappedcnt=nummappedchr
                minmappedchr=line[0]
        cmd='%s/samtools view %s %s'%(self.tooldirs[0],self.name,minmappedchr)
        stdo=os.popen(cmd)
        for linetxt in stdo:
            line=linetxt.rstrip('\n').split('\t')
            break
        readsize=len(line[9])
        return (nummapped,readsize)   
    
    def Tobdg(self,outname):
        bdgfile='/'.join(self.name.split('/')[:-1])+'/'+outname+'.bdg'
        if os.path.exists(bdgfile):
            os.system('rm %s'%bdgfile)
        cmd='%s/bedtools genomecov -split -bg -ibam %s > %s'%(self.tooldirs[1],self.name,bdgfile)
        os.system(cmd)

    def Tojun(self,outname):
        junfile='/'.join(self.name.split('/')[:-1])+'/'+outname+'.jun'
        if os.path.exists(junfile):
            os.system('rm %s'%junfile)
        fjun=open(junfile,'w')
        chrlist=self.getchrlist()
        for chrnm in chrlist:
            cmd='%s/samtools view %s %s'%(self.tooldirs[0],self.name,chrnm)
            stdo=os.popen(cmd)
            alljuncwtlist=[]
            for linetxt in stdo:
                line=linetxt.rstrip('\n').split('\t')
                cigar=line[5]
                cigarstr=common.convgenericcigar(cigar)
                juncexon=common.cigar2juncexon(cigarstr)
                if len(juncexon[0])==0:
                    continue
                else:
                    refbase=int(line[3])
                    junclist=[(junc[0]+refbase-1,junc[1]+refbase) for junc in juncexon[0]]
                    for junc in junclist:
                        if len(alljuncwtlist)==0:
                            alljuncwtlist.append([junc,1])
                        else:    
                            if junc==alljuncwtlist[-1][0]:
                                alljuncwtlist[-1][1]+=1                             
                            else:
                                alljuncwtlist.append([junc,1])
            alljuncwtlist.sort()
            j=len(alljuncwtlist)-1
            while j>0:
                if alljuncwtlist[j][0]==alljuncwtlist[j-1][0]:
                    alljuncwtlist[j-1][1]+=alljuncwtlist[j][1]
                    del alljuncwtlist[j]
                j-=1
            for juncwt in alljuncwtlist:
                fjun.write('%s\t%d\t%d\t%d\n'%(chrnm,juncwt[0][0],juncwt[0][1],juncwt[1]))
        fjun.close()
        
    def _alignlinesngraphstruct2wtgraphstruct(self,graphstruct,alignlines):
        exonlist=graphstruct[0]; splicelist=graphstruct[2]
        exonwtlist,splicewtlist=common.aligntuples2coverage(alignlines,exonlist,splicelist)
        return [exonlist,[],splicelist,[],[],[],exonwtlist,[],splicewtlist]
    
    def _alignlinesndivdict2flowdict(self,alignlines,divdict):
        '''
        Alert: graph struct is created only for nodes with flows
        May need improvement with introns - using flow file to find divergences
        FIXME
        '''
        graphstruct=common.divdict2graphstruct(divdict)
        wtgs=self._alignlinesngraphstruct2wtgraphstruct(graphstruct, alignlines)
        wgs=actgraph.wtgraphstruct(wtgs)
        return wgs.divdictToflow(divdict)
    
    def Togeneflowlistdict(self,genedivdict,islanddict,numpart):
        '''
        Also splits it into numparts
        flowdict[position]=[[outgoingflg,exonstartflg],[[wt1,wt2],nflowvec]]
        geneflowdict=gene:pos:[flowdictlist]
        '''
        geneflowlistdict={}
        for gene in genedivdict.keys():
            divdict=genedivdict[gene]
            key='%s:%d-%d'%(tuple(islanddict[gene][0:3]))
            alignlines=self.retrieve(key) 
            alignlineparts=common.randgroup(alignlines,numpart)
            flowdictlist=[]
            flowlistdict={}
            for alignlinepart in alignlineparts:
                flowdict=self._alignlinesndivdict2flowdict(alignlinepart, divdict)
                flowdictlist.append(flowdict)
            geneflowlistdict[gene]=flowdictlist             
        return geneflowlistdict
