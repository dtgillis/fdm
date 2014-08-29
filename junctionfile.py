import os
import cPickle

class junFile:
    '''Class to index and retrieve data from large files
    input is file and chr:start-end and gene for act file
    1) bam : uses samtools to build index and retrieve data
    2) bdg : bedgraph file chr start end wt
    3) jun : junction file chr start end wt
    4) act : actg file : gff format: only indexed by genes
    samtools and bedtools paths are known
    '''
    def __init__(self, filename,tooldirs):
        self.tooldirs=tooldirs
        self.name=filename
        self.ext=filename[-3:]
        self.idxname=filename+'.'+ self.ext[:-1]+'i'
        if not self.checkindex():
            self.buildindex()
    
    def checkindex(self):
        if os.path.isfile(self.idxname):
            if os.path.getmtime(self.idxname) > os.path.getmtime(self.name) and os.path.getsize(self.name)>0:
                return 1
        return 0
    
    def buildindex(self):
        indexdict={}
        total=0
        for line in open(self.name):
            ln=line.rstrip('\n').split('\t')
            chrnm=ln[0]; chstr=int(ln[1]); chrstr1M = chstr/1000000
            if chrnm not in indexdict.keys():
                indexdict[chrnm]={}
                indexdict[chrnm][chrstr1M]=[total,1]
            else:
                if chrstr1M not in indexdict[chrnm].keys():
                    indexdict[chrnm][chrstr1M]=[total,1]
                else:
                    indexdict[chrnm][chrstr1M][1]+=1
            total += len(line)
        cPickle.dump(indexdict,open(self.idxname,'w'))
        
    def retrieve(self,key):
        outlines=[]
        indexdict=cPickle.load(open(self.idxname))
        chrnm=key.split(':')[0]
        chrstr=int(key.split(':')[1].split('-')[0])
        chrend=int(key.split(':')[1].split('-')[1]) 
        chrstr1M=chrstr/1000000
        chrend1M=chrend/1000000
        if chrnm not in indexdict.keys():
            return []
        offset=-1; numlines=0
        for chrstrkey in range(chrstr1M,chrend1M+1):
            if chrstrkey in indexdict[chrnm]:
                if offset==-1:
                    offset=indexdict[chrnm][chrstrkey][0]
                    numlines=indexdict[chrnm][chrstrkey][1]
                else:
                    numlines+=indexdict[chrnm][chrstrkey][1]
        if offset==-1:
            return []
        filein=open(self.name)
        filein.seek(offset)
        for i in range(numlines):
            line = filein.readline().rstrip('\n').split('\t')
            line[1]=int(line[1]); line[2]=int(line[2]); line[3]=int(line[3])
            if min(line[2],chrend)>=max(line[1],chrstr):
                outlines.append(line)
            
        return outlines    
    
    def _chrjuncdict2chrsplicenodedict(self,chrjuncdict):
        '''
        Input: {chr:strand:splicelist}
        Output:{chr:[node1list,node2list]}
        '''
        chrsplicenodedict={}
        for chr in chrjuncdict:
            splicelist=chrjuncdict[chr][0]+chrjuncdict[chr][1]
            node1list=[splice[0] for splice in splicelist]
            node2list=[splice[1] for splice in splicelist]
            chrsplicenodedict[chr]=[node1list,node2list]
        return chrsplicenodedict
            
    def Tobed(self,annochrjuncdict={}):
        '''
        for every gene, put a different score
        annochrjuncdict has list of annotated junctions for correct naming
        '''

        outfile=self.name[:-3]+'bed'
        fout=open(outfile,'w')
        #annoflsz=5; novelflsz=2
        #annoprf='A_';novelprf='N_'
        #replaced by dict
        flankszdict={-1:5,0:2,1:5,2:5,3:5}
        prfdict={-1:'A',0:'N',1:'A',2:'S',3:'P'}        
        annochrsplicenodedict=self._chrjuncdict2chrsplicenodedict(annochrjuncdict)
        for juntxt in open(self.name):
            jun=juntxt.rstrip('\n').split('\t')
            chr=jun[0]; jun1=int(jun[1]); jun2=int(jun[2])
            #print chr,jun1,jun2,score
            score=int(jun[3])
            if chr not in annochrjuncdict.keys():
                annoflg=0
            else:
                if [jun1,jun2] in  annochrjuncdict[chr][0]:
                    annoflg=-1
                elif [jun1,jun2] in  annochrjuncdict[chr][1]:
                    annoflg=1
                elif jun1 in annochrsplicenodedict[chr][0] and jun2 in annochrsplicenodedict[chr][1]:
                    annoflg=2
                elif jun1 in annochrsplicenodedict[chr][0] or jun2 in annochrsplicenodedict[chr][1]:
                    annoflg=3
                else:
                    annoflg=0
            strand=1
            prf=prfdict[annoflg]
            flanksz=flankszdict[annoflg]
            if annoflg==1:
                strand=1
            block2sz=jun2-jun1-flanksz
            if score !=0:
                fout.write('%s\t%d\t%d\t%s_%05d\t%d\t%s\t0\t0\t0\t2\t%d,%d\t0,%d\n'% \
                           (chr,jun1,jun2,prf,score,score,strand,flanksz,flanksz,block2sz))
        fout.close()
    
    def Tobigbed(self,annochrjuncdict={}):
        bedfile=self.name[:-3]+'bed'
        self.Tobed(annochrjuncdict)
        bbfile=bedfile[:-3]+'bb'
        cmd='%s/bedToBigBed %s %s %s'%(self.tooldirs[0],bedfile,self.tooldirs[1],bbfile)
        os.system(cmd)
