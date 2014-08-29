import os
import common
import cPickle

class gtfFile:
    '''
    This processes the GTF files
    Two types of files are processed:
        Ensembl type
        Provided
    '''
    def __init__(self, filename,redo_flg=0,type='F'):
        '''
        type='F' for provided gtf file
            ='E' for ENSMBL
            ## Give examples
        '''
        self.name=filename
        self.type=type
        self.dir='%s/fdm_%s'%('/'.join(filename.split('/')[:-1]),filename.split('/')[-1])
        
        genetranscriptdict_file='%s/genetranscriptdict.pck'%self.dir
        islandlist_file='%s/islandlist.pck'%self.dir
        if not os.path.isfile(genetranscriptdict_file) or not os.path.isfile(islandlist_file) or redo_flg==1:
            do_flg=1
        else:
            do_flg=0
        if do_flg:
            if not os.path.isdir(self.dir):
                os.system('mkdir %s'%self.dir)
            genetranscriptdict=self._gtf2genetranscriptdict()
            cPickle.dump(genetranscriptdict,open(genetranscriptdict_file,'w'))  
            if type=='F':
                islandlist=self._gtf2islandlist()
            elif type=='E':
                islandlist=self._genetranscriptdict2islandlist(genetranscriptdict)
            cPickle.dump(islandlist,open(islandlist_file,'w'))   
        self._num_gtd_accessed=0
        self.genetranscriptdict={}   
        self._num_gagd_accessed=0
        self.geneannogrsdict={}
        self._num_gadd_accessed=0
        self.geneannodivdict={}
        
    def getchrlist(self):
        chrlist=[]
        for linetxt in open(self.name):
            line=linetxt.rstrip('\n').split('\t')
            if line[0] not in chrlist:
                chrlist.append(line[0])
        chrlist.sort()
        return chrlist
    
    def _attrtxt2attrdict(self,attrtxt):
        attrlist=attrtxt.split(';')
        attrtuplist=[attr.strip(' ').split(' ',1) for attr in attrlist]
        attrtuplist=[(attrtup[0],attrtup[1][1:-1]) for attrtup in attrtuplist if len(attrtup)==2]
        return dict(attrtuplist)
    
    def _gtf2islandlist(self):
        '''
        For repeated genes, if regions do not match ignore
        '''
        islanddict={}
        rej_genelist=[]
        for linetxt in open(self.name):
            line=linetxt.rstrip('\n').split('\t')
            if line[2]!='gene':
                continue
            attrdict=self._attrtxt2attrdict(line[8])
            gene=attrdict['gene_id']
            if line[6]=='+':
                strand=1
            else:
                strand=0
            gstart=int(line[3]); gend=int(line[4])
            chrid=line[0]
            if gene not in islanddict.keys() and gene not in rej_genelist:
                islanddict[gene]=[chrid,gstart,gend,gene,strand]
            else:
                if gene in rej_genelist:
                    continue
                intv1=[chrid,gstart,gend,gene,strand]
                intv2=islanddict[gene]
                intv=self._intervalintersect(intv1,intv2)
                if intv==0:
                    rej_genelist.append(gene)
                else:
                    islanddict[gene]=intv
        message='Number of repeated genes rejected: %d'%len(rej_genelist)
        common.printstatus(message,'W',common.func_name())
        message='Number of genes loaded: %d'%(len(islanddict.keys()))
        common.printstatus(message,'S',common.func_name())
        for gene in rej_genelist:
            del islanddict[gene]
        islandlist=[]
        for gene in islanddict.keys():
            islandlist.append(islanddict[gene])
        islandlist.sort()
        return islandlist
    
    def getislandlist(self):
        '''
        chrid,gstart,gend,gene,strand
        '''
        islandlist_file='%s/islandlist.pck'%self.dir
        if os.path.isfile(islandlist_file):
            islandlist=cPickle.load(open(islandlist_file))
        else:
            islandlist=self._gtf2islandlist()
            cPickle.dump(islandlist,open(islandlist_file,'w'))
        return islandlist
    
    def _intervalintersect(self,intv1,intv2):
        '''
        chrid,start,end,strand
        '''
        if intv1[0]!=intv2[0] or intv1[4]!= intv2[4]:
            return 0
        else:
            if ( intv1[1]>=intv2[1] and intv1[1]<=intv2[2] ) or ( intv2[1]>=intv1[1] and intv2[1]<=intv1[2] ) :
                return [intv1[0], min(intv1[1],intv2[1]),max(intv1[2],intv2[2]),intv1[3],intv1[4]]
            else:
                return 0    

    def getislanddict(self):
        islandlist=self.getislandlist()
        islanddict={}
        for island in islandlist:
            #islanddict[island[3]] = list(island[0:3])+[island[4]]
            islanddict[island[3]] = island
        return islanddict        

    def _gtf2genetranscriptdict(self):
        if self.type=='E':
            geneidx='gene_name'
        elif self.type=='F':
            geneidx='gene_id'
        genetranscriptdict={}
        lnum=0
        for linetxt in open(self.name):
            lnum+=1
            if lnum%200000==1:
                message='Processing GTF file for transcripts at line %d'%lnum
                common.printstatus(message,'S',common.func_name())
            line=linetxt.rstrip('\n').split('\t')
            if line[2]!='exon':
                continue
            attrdict=self._attrtxt2attrdict(line[8])
            gene=attrdict[geneidx]
            #print attrdict
            transcript=attrdict['transcript_id']
            chrnm=line[0]
            if line[6]=='+':
                strand=1
            else:
                strand=0
            #print gene, transcript, chrnm, strand
            if gene not in genetranscriptdict.keys():
                genetranscriptdict[gene]={}
            if transcript not in  genetranscriptdict[gene].keys():
                genetranscriptdict[gene][transcript]=(chrnm,strand,[])
            genetranscriptdict[gene][transcript][2].append((int(line[3]),int(line[4])))
        for gene in genetranscriptdict:
            for transcript in genetranscriptdict[gene]:
                genetranscriptdict[gene][transcript][2].sort()
        message='Completed processing GTF file for transcripts'
        common.printstatus(message,'S',common.func_name())
        
        genetranscriptdict_file='%s/genetranscriptdict.pck'%self.dir
        cPickle.dump(genetranscriptdict,open(genetranscriptdict_file,'w'))
        genetranscriptdict=self._fixgenetranscriptdict()
        return genetranscriptdict
    
    def _genetranscriptdict2islandlist(self,genetranscriptdict):
        genelist=genetranscriptdict.keys()
        islandlist=[]
        for gene in genelist:
            txlist=genetranscriptdict[gene].keys()
            startlist=[];endlist=[]
            for tx in txlist:
                chrnm=genetranscriptdict[gene][tx][0]
                strand=genetranscriptdict[gene][tx][1]
                startlist.append(genetranscriptdict[gene][tx][2][0][0])
                endlist.append(genetranscriptdict[gene][tx][2][-1][1])
            start=min(startlist)
            end=max(endlist)
            islandlist.append([chrnm,start,end,gene,strand])
        islandlist.sort()
        return islandlist
        
    def _getgenetranscriptdict(self):
        '''
        ABAT : uc002czd.3 : ('16', 1, [(8806826, 8807722), (8829556, 8829666), (8839858, 8839955), (8841965, 8841994)])
        '''
        if self.genetranscriptdict=={}:
            self._num_gtd_accessed+=1
            genetranscriptdict_file='%s/genetranscriptdict.pck'%self.dir
            if os.path.isfile(genetranscriptdict_file):
                genetranscriptdict=cPickle.load(open(genetranscriptdict_file))
            else:
                if not os.path.isdir(self.dir):
                    os.system('mkdir %s'%self.dir)
                genetranscriptdict=self._gtf2genetranscriptdict()  
            if self._num_gtd_accessed>2:
                self.genetranscriptdict=genetranscriptdict
            return genetranscriptdict
        else:
            return self.genetranscriptdict
        
    def getgenetranscriptdict(self):
        return self._getgenetranscriptdict()
    
    def _getgenetranscriptsizelist(self,gene):
        genetranscriptdict=self._getgenetranscriptdict()
        transcriptdict=genetranscriptdict[gene]
        sizelist=[]
        tkeys=transcriptdict.keys()
        for tkey in tkeys:
            exonlist=transcriptdict[tkey][2]
            size=0
            for exon in exonlist:
                size+=(exon[1]-exon[0]+1)
            sizelist.append(size)
        sizelist.sort()
        return sizelist
    
    def getmediantranscriptsize(self,gene):
        sizelist=self._getgenetranscriptsizelist(gene)
        mid=len(sizelist)/2
        return sizelist[mid]
    
    def getlargesttranscriptsize(self,gene):
        sizelist=self._getgenetranscriptsizelist(gene)
        return sizelist[-1]
    
    def getnumtranscript(self,gene):
        sizelist=self._getgenetranscriptsizelist(gene)
        return len(sizelist)
    
    def _fixgenetranscriptdict(self):
        genetranscriptdict_file='%s/genetranscriptdict.pck'%self.dir
        genetranscriptdict=cPickle.load(open(genetranscriptdict_file))
        newgenetranscriptdict={}
        for gene in genetranscriptdict:
            trascriptdict=genetranscriptdict[gene]
            transcriptkeys=trascriptdict.keys()
            transcriptkeys.sort()
            chrstrand=trascriptdict[transcriptkeys[0]][0:2]
            removelist=[]
            for key in transcriptkeys[1:]:
                if chrstrand!=trascriptdict[key][0:2]:
                    message='Gene %s in multiple chromosomes; processing only one location:%s'%(gene,str(chrstrand))
                    common.printstatus(message,'W',common.func_name())
                    removelist.append(key)
            for key in removelist:
                transcriptkeys.remove(key)
            txrangelist=[[trascriptdict[key][2][0][0],trascriptdict[key][2][-1][1]] for key in transcriptkeys]
            ztxrangekey=zip(txrangelist,transcriptkeys)
            ztxrangekey.sort()
            txrangelist=[x[0] for x in ztxrangekey]
            transcriptkeys=[x[1] for x in ztxrangekey]
            trange=txrangelist[0]
            for i in range(1,len(txrangelist)):
                txrange=txrangelist[i]
                if txrange[0]<=trange[1]:
                    trange=[trange[0],max(trange[1],txrange[1])]
                else:
                    message='Gene %s: multiple non-overlapping transcription regions; processing only one region:%s'%(gene,str(trange))
                    common.printstatus(message,'W',common.func_name())
                    transcriptkeys=transcriptkeys[0:i]
                    break
            newgenetranscriptdict[gene]={}
            for transcript in transcriptkeys:
                newgenetranscriptdict[gene][transcript]=genetranscriptdict[gene][transcript]
            
        cmd='mv %s %s.bk'%(genetranscriptdict_file,genetranscriptdict_file)
        os.system(cmd)
        cPickle.dump(newgenetranscriptdict,open(genetranscriptdict_file,'w'))
        return newgenetranscriptdict
        
    def _genetranscriptdict2graphstruct(self,gene):
        '''
        get a graph structure for a gene
        Alert: Exon Ends may need correction
        ''' 
        transcriptexonlist=[]
        genetranscriptdict=self._getgenetranscriptdict()
        trascriptdict=genetranscriptdict[gene]
        transcriptkeys=trascriptdict.keys()
        transcriptkeys.sort()
        chrstrand=trascriptdict[transcriptkeys[0]][0:2]
        removelist=[]
        for key in transcriptkeys[1:]:
            if chrstrand!=trascriptdict[key][0:2]:
                message='Gene %s in multiple chromosomes; processing only one location:%s'%(gene,str(chrstrand))
                common.printstatus(message,'W',common.func_name())
                removelist.append(key)
        for key in removelist:
            transcriptkeys.remove(key)
        txrangelist=[[trascriptdict[key][2][0][0],trascriptdict[key][2][-1][1]] for key in transcriptkeys]
        ztxrangekey=zip(txrangelist,transcriptkeys)
        ztxrangekey.sort()
        txrangelist=[x[0] for x in ztxrangekey]
        transcriptkeys=[x[1] for x in ztxrangekey]
        trange=txrangelist[0]
        for i in range(1,len(txrangelist)):
            txrange=txrangelist[i]
            if txrange[0]<=trange[1]:
                trange=[trange[0],max(trange[1],txrange[1])]
            else:
                message='Gene %s: multiple non-overlapping transcription regions; processing only one region:%s'%(gene,str(trange))
                common.printstatus(message,'W',common.func_name())
                transcriptkeys=transcriptkeys[0:i]
                break
  
        splicelist=[]
        startnodelist=[]
        endnodelist=[]
        for key in transcriptkeys:
            exonlist=trascriptdict[key][2]
            #print exonlist
            if exonlist[0][0] not in startnodelist:
                startnodelist.append(exonlist[0][0])
            if exonlist[-1][1] not in endnodelist:
                endnodelist.append(exonlist[-1][1])
            for exon in exonlist:
                if exon not in transcriptexonlist:
                    transcriptexonlist.append(exon)
            for i in range(1,len(exonlist)):
                splice = [exonlist[i-1][1],exonlist[i][0]]
                if splice not in splicelist:
                    splicelist.append(splice)
        exonlist=[]
        transcriptexonlist.sort()
        transcriptexonlist=[list(exon) for exon in transcriptexonlist]
        splicelist.sort()
        startnodelist.sort()
        endnodelist.sort()
        #print transcriptexonlist
        exonqueue=transcriptexonlist[0:1]
        for i in range(len(transcriptexonlist)-1):
            exon=transcriptexonlist[i+1]
            #print exonqueue, exon, exonlist
            tqueue=exonqueue[0:]
            for exonregion in tqueue:
                if exonregion[1]<exon[0]:
                    exonlist.append(exonregion)
                    exonqueue.remove(exonregion)
            found=0
            tqueue=exonqueue[0:]
            for exonregion in tqueue:
                if exon[0] > exonregion[0] and exon[0] < exonregion[1]:
                    if exon[1] > exonregion[0] and exon[1] < exonregion[1]:
                        exonqueue.append([exonregion[0],exon[0]-1])
                        exonqueue.append([exon[0],exon[1]])
                        exonqueue.append([exon[1]+1,exonregion[1]])
                        exonqueue.remove(exonregion)
                    else:
                        exonqueue.append([exonregion[0],exon[0]-1])
                        exonqueue.append([exon[0],exonregion[1]])
                        exonqueue.remove(exonregion)
                elif exon[1] > exonregion[0] and exon[1] < exonregion[1]:
                    exonqueue.append([exonregion[0],exon[1]])
                    exonqueue.append([exon[1]+1,exonregion[1]])
                    exonqueue.remove(exonregion)
            exonqueue.sort()     
            if len(exonqueue)>0:   
                if exon[1]>exonqueue[-1][1]:
                    exonqueue.append([exonqueue[-1][1]+1,exon[1]])
            else:
                exonqueue.append(exon)
        for exonregion in exonqueue:
            if exonregion[0]<exonregion:
                exonlist.append(exonregion)
        exonlist.sort()
        intronlist=[]
        novelnodelist=[]
        for i in range(len(exonlist)-1):
            if exonlist[i][1]+1<exonlist[i+1][0]-1:
                if [exonlist[i][1],exonlist[i+1][0]] not in splicelist:
                    intronlist.append([exonlist[i][1]+1,exonlist[i+1][0]-1])
                else:
                    dummyleftnode=(exonlist[i][1]+1+exonlist[i+1][0]-1)/2
                    dummyrightnode=dummyleftnode+1
                    intronlist.append([exonlist[i][1]+1,dummyleftnode])
                    intronlist.append([dummyrightnode,exonlist[i+1][0]-1])
                    novelnodelist.append(dummyleftnode)
#        print exonlist
#        print intronlist
#        print splicelist
#        print startnodelist
#        print endnodelist
        return (exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist)
   


    def getgene2annographstruct(self,gene):
        geneannogrsdict=self._getgeneannogrsdict()
        graphstruct=geneannogrsdict[gene]
        return graphstruct
        
    def getgene2annodivdict(self, gene):
        geneannodivdict=self.getgeneannodivdict()
        divdict=geneannodivdict[gene]
        return divdict
    
    def _getgeneannogrsdict(self):
        if self.geneannogrsdict=={}:
            self._num_gagd_accessed+=1
            geneannogrsdict_file='%s/geneannogrsdict.pck'%self.dir
            if os.path.isfile(geneannogrsdict_file):
                #print geneannogrsdict_file
                geneannogrsdict=cPickle.load(open(geneannogrsdict_file))
            else:
                geneannogrsdict={}
                islanddict=self.getislanddict()
                lnum=0
                for gene in islanddict.keys():
                    lnum+=1
                    if lnum%1000==1:
                        message='Processing Graph Struct construction at gene %d'%lnum
                        common.printstatus(message,'S',common.func_name())
                    geneannogrsdict[gene]=self._genetranscriptdict2graphstruct(gene)
                cPickle.dump(geneannogrsdict,open(geneannogrsdict_file,'w'))  
                message='Completed Processing Graph Struct construction'
                common.printstatus(message,'S',common.func_name())
            if self._num_gagd_accessed>2:
                self.geneannogrsdict=geneannogrsdict
            return geneannogrsdict 
        else:
            return self.geneannogrsdict
        
    def getgeneannogrsdict(self):
        return self._getgeneannogrsdict()
   
    def getgeneannodivdict(self):
        if self.geneannodivdict=={}:
            self._num_gadd_accessed+=1
            geneannodivdict_file='%s/geneannodivdict.pck'%self.dir
            if os.path.isfile(geneannodivdict_file):
                geneannodivdict=cPickle.load(open(geneannodivdict_file))
            else:
                geneannodivdict={}
                islanddict=self.getislanddict()
                lnum=0
                for gene in islanddict.keys():
                    lnum+=1
                    if lnum%1000==1:
                        message='Processing Divergence Struct construction at gene %d:%s'%(lnum,gene)
                        common.printstatus(message,'S',common.func_name())
                    graphstruct=self.getgene2annographstruct(gene)
                    divdict=common.graphstruct2divdict(graphstruct)
                    geneannodivdict[gene]=divdict
                cPickle.dump(geneannodivdict,open(geneannodivdict_file,'w'))  
                message='Completed Processing Divergence Struct construction'
                common.printstatus(message,'S',common.func_name())
            if self._num_gadd_accessed>2:
                self.geneannodivdict=geneannodivdict
            return geneannodivdict
        else:
            return self.geneannodivdict

    def getchrjuncdict(self):
        genetranscriptdict_file='%s/genetranscriptdict.pck'%self.dir
        genetranscriptdict=cPickle.load(open(genetranscriptdict_file))
        chrjuncdict={}
        for gene in genetranscriptdict.keys():
            tx1=genetranscriptdict[gene].keys()[0]
            chr='%s'%genetranscriptdict[gene][tx1][0]
            #chr=chr.replace('chrchr','chr')
            strand=genetranscriptdict[gene][tx1][1]
            graphstruct=self._genetranscriptdict2graphstruct(gene)
            splicelist=graphstruct[2]
            if chr in chrjuncdict.keys():
                if strand in chrjuncdict[chr].keys():
                    chrjuncdict[chr][strand]+=splicelist[0:]
            else:
                chrjuncdict[chr]={}
                chrjuncdict[chr][strand]=splicelist[0:]
                chrjuncdict[chr][1-strand]=[]
        return chrjuncdict
