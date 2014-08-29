import os
import cPickle

class bdgFile:
    def __init__(self, filename,tooldirs):
        self.tooldirs=tooldirs
        self.name=filename
        self.ext=filename[-3:]
        self.idxname=filename+'.'+ self.ext[:-1]+'i'
        if not self.checkindex():
            self.buildindex()
    
    def checkindex(self):
        if os.path.isfile(self.idxname):
            if os.path.getmtime(self.idxname) > os.path.getmtime(self.name):
                return 1
        return 0
    
    def buildindex(self):
        indexdict={}
        total=0
        for line in open(self.name):
            ln=line.rstrip('\n').split('\t')
            chrnm=ln[0]; chstr=int(float(ln[1])); chrstr1M = chstr/1000000
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
        found=0
        for i in range(numlines):
            line = filein.readline().rstrip('\n').split('\t')
            line[1]=int(float(line[1])); line[2]=int(float(line[2])); line[3]=int(float(line[3]))
            if min(line[2],chrend)>=max(line[1],chrstr):
                outlines.append(line)
                found=1
            elif found==1:
                break
        return outlines   
       
    def _bdgremovedup(self):
        os.system('mv %s %s.tmp'%(self.name,self.name))
        bdgfilename='%s.tmp'%(self.name)
        fout=open(self.name,'w')
        prevpos1=1
        for linetxt in open(bdgfilename):
            line=linetxt.split('\t')
            pos1=int(line[1])
            if pos1!=prevpos1:
                fout.write(linetxt)
            prevpos1=pos1
        os.system('rm %s.tmp'%(self.name))
            
    def bdgsort(self):
        self._bdgremovedup()
        os.system('mv %s %s.tmp2'%(self.name,self.name))
        cmd="sort -t$'\t' -k 1,1 -k 2n,2n %s.tmp2 > %s"%(self.name,self.name)
        os.system(cmd)
        os.system('rm %s.tmp2'%(self.name))

    
    def Tobigwig(self):
        bwfile=self.name[:-3]+'bw'
        #cmd='cp %s %s.orig'%(self.name,self.name)
        #os.system(cmd)
        cmd='%s/bedGraphToBigWig %s %s %s'%(self.tooldirs[0],self.name,self.tooldirs[1],bwfile)
        os.system(cmd)
