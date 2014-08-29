import os
import cPickle
import common
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib.patches import Arc
import cvxpy as cvx
        
class actFile:
    '''Class to index and retrieve data from large files
    input is file and chr:start-end and gene for act file
    1) bam : uses samtools to build index and retrieve data
    2) bdg : bedgraph file chr start end wt
    3) jun : junction file chr start end wt
    4) act : actg file : gff format: only indexed by genes
    samtools and bedtools paths are known
    '''
    def __init__(self, filename):
        self.name=filename
        self.ext=filename[-3:]
        self.idxname=filename+'.'+ self.ext[:-1]+'i'
        if not self.checkindex():
            self.buildindex()
        wgsfilename=self.name[0:-3]+'wgs'
        if os.path.isfile(wgsfilename):
            self.wgsdict=cPickle.load(open(wgsfilename))
        else:
            self.wgsdict={}
    
    def checkindex(self):
        if os.path.isfile(self.idxname):
            if os.path.getmtime(self.idxname) > os.path.getmtime(self.name):
                return 1
        return 0
    
    def buildindex(self):
        indexdict={}
        total=0
        for line in open(self.name):
            gene=line.rstrip('\n').split('\t')[8].split(':')[0]
            if gene[0:2] not in indexdict.keys():
                indexdict[gene[0:2]]={}
                indexdict[gene[0:2]][gene]=[total,1]
            else:
                if gene not in indexdict[gene[0:2]].keys():
                    indexdict[gene[0:2]][gene]=[total,1]
                else:
                    indexdict[gene[0:2]][gene][1]+=1
            total += len(line)
        cPickle.dump(indexdict,open(self.idxname,'w'))
        
    def retrieve(self,key):
        outlines=[]
        indexdict=cPickle.load(open(self.idxname))
        gene=key
        if gene[0:2] not in indexdict.keys():
            return []
        elif gene not in indexdict[gene[0:2]].keys():
            return []
        else:
            offset,numlines=indexdict[gene[0:2]][gene]
            filein=open(self.name)
            filein.seek(offset)
            for i in range(numlines):
                line = filein.readline().rstrip('\n').split('\t')
                outlines.append(line)
            return outlines

    def getwtgraphstruct(self,gene):
        '''Text file to interpretation
        '''
        if len(self.wgsdict)!=0:
            return self.wgsdict.get(gene,[])
        actlines=self.retrieve(gene)
        exonlist=[]
        intronlist=[]
        splicelist=[]
        startnodelist=[]
        endnodelist=[]
        novelnodelist=[]
        exonwtlist=[]
        intronwtlist=[]
        splicewtlist=[]
        for actline in actlines:
            source=actline[1]
            feature=actline[2]
            startn=int(actline[3])
            endn=int(actline[4])
            score=float(actline[5])
            if feature=='node':
                if source=='annot':
                    startnodelist.append(startn)
                elif source=='annoe':
                    endnodelist.append(endn)
                elif source=='novel':
                    novelnodelist.append(startn)
            if feature=='exon':
                #if len(exonlist)!=0:
                #    if exonlist[-1][1]<startn-1:
                #        intron=[exonlist[-1][1]+1,startn-1]
                #        if intron not in intronlist:
                #            intronlist.append(intron)
                #            intronwtlist.append(0)
                exon=[startn,endn]
                exonlist.append(exon)
                exonwtlist.append(score)
            if feature=='splice':
                splice=[startn,endn]
                splicelist.append(splice)
                splicewtlist.append(score)        
            if feature=='retint':
                intron=[startn,endn]
                intronlist.append(intron)
                intronwtlist.append(score)      
    #    print zip(exonlist,exonwtlist)
    #    print zip(intronlist,intronwtlist)
    #    print zip(splicelist,splicewtlist)
    #    print startnodelist
    #    print endnodelist
    #    print novelnodelist    
        return (exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist)


    def Toimage(self,gene, imagedir,genomerange=[0,10000000000],highlightnodelist=[],ext='pdf',inwgs=[]):
        
        ts=self.name.split('/')[-1][:-4]
        #image rendering constants
        hconst=30;gwidth=100.0;exonplotdistbasis=10.0;msize=10;hoffset=10;woffset=50
        
        if inwgs==[]:
            wtgraphstruct=self.getwtgraphstruct(gene)    
        else:
            wtgraphstruct=inwgs
        exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist=wtgraphstruct
        
        #print zip(exonlist,exonwtlist)
        #print zip(intronlist,intronwtlist)
        #print zip(splicelist,splicewtlist)
        #print startnodelist,endnodelist,novelnodelist
 
        
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
                awidth=gwidth/gheight*(y2-y1)
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

        
    def getcorrectionwt(self,edgesz, txstend, edgewt, readsz,edgetype):
        """
            Function to compute penalty for weight change: more means expensive to change value
            edgesz= 0 for splice,  total genomic size for exons/retained intron
            txstend = 1 if exon is first or last of a transcript (might change to set start and end appropriately
            readsz = size of read
            weight range = 1 to 4
        """
        k = 1.5
        edgesz=max(0,edgesz-txstend*readsz)
        if edgetype==0 and edgewt<1:
            return 25
        if edgetype==2 and edgewt<4.0:
            return 20.0/max(edgewt,1.0)
        
        return min(4, 2+k*(math.sqrt(edgesz*1.0/readsz))/math.log(2+edgewt,2)) 

    def wgs2problem(self,wgs,readsz=100):
        """
        X (nx1) = Variables are exon wts, splice wts and tx start and end wts
        B (mx1)=  0 for each node comparing inflow and outflow
               = given wt for each exon,splice*correction weight
               = weq/10 * inflow outflow difference for tx start or end
        A (mxn) = set of equations
        """
        weq=100.0
        X=[]; B=[]; A=[]
        
        exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist=wgs
        sparsegraphdict,nodedict=wgs2sparsegraphdict(wgs)
        
        #print sparsegraphdict
        #print nodedict
        
        txstartlist=[nodedict[x][0] for x in nodedict.keys() if nodedict[x][1]==2]
        txendlist=[nodedict[x][0] for x in nodedict.keys() if nodedict[x][1]==3]
        txstartlist.sort()
        txendlist.sort()
        #print txstartlist
        #print txendlist
        
        exonwttuplelist=zip(exonlist+intronlist,exonwtlist+intronwtlist)
        splicewttuplist=zip(splicelist,splicewtlist)
        exonwttuplelist.sort()
        splicewttuplist.sort()
        
        #print exonwttuplelist
        #print splicewttuplist
    
        exondict={}; splicedict={}
        exonidx=0;spliceidx=0
        for node1 in sorted(sparsegraphdict.keys()):
            for node2 in sorted(sparsegraphdict[node1].keys()):
                if sparsegraphdict[node1][node2][1]==1:
                    exondict[(node1,node2)]=exonidx
                    exonidx+=1
                if sparsegraphdict[node1][node2][1]==2:
                    splicedict[(node1,node2)]=spliceidx
                    spliceidx+=1
        
        incomingedgedict={}
        for node1 in sorted(sparsegraphdict.keys()):
            for node2 in sorted(sparsegraphdict[node1].keys()):
                if node2 not in incomingedgedict:
                    incomingedgedict[node2]={}
                incomingedgedict[node2][node1]=sparsegraphdict[node1][node2][0:]
        
        #print exondict
        #print splicedict
        #print incomingedgedict
        
        numvar=len(exonwttuplelist)+len(splicewttuplist)+len(txstartlist)+len(txendlist)
        for exon,wt in exonwttuplelist:
            X.append([2,exon[0],exon[1],wt])
            txstend=0
            if exon[0] in txstartlist:
                txstend+=1
            if exon[1] in txendlist:
                txstend+=1     
            if exon in exonlist:
                edgetype=1
            else:
                edgetype=2  
            w=self.getcorrectionwt(exon[1]-exon[0], txstend, wt, readsz,edgetype)
            row=[0]*numvar
            row[len(X)-1]=w
            A.append(row)
            B.append(w*wt)
            #B.append(w*math.log10(1.0+exon[2]))
        for splice,wt in splicewttuplist:     
            X.append([3,splice[0],splice[1],wt])  
            edgetype=0
            w=self.getcorrectionwt(0, 0, wt, readsz,edgetype)  
            row=[0]*numvar
            row[len(X)-1]=w
            A.append(row)
            B.append(w*wt)   
        for node in txstartlist:
            X.append([1,node,node,0.0])
        for node in txendlist:
            X.append([-1,node,node,0.0])    
    
        nodelist=sorted(nodedict.keys())
        for node in nodelist:
            if (nodelist.index(node)==0) and (nodedict[node][1] not in [2,3]):
                continue
            if (nodelist.index(node)==len(nodelist)-1) and (nodedict[node][1] not in [2,3]):
                continue            
            row=[0]*numvar
            txstflg=0; txendflg=0; inwt=0.0; outwt=0.0
            if nodedict[node][1]==2:
                idx=len(exondict)+len(splicedict)+txstartlist.index(nodedict[node][0])
                row[idx]=-1*weq
                txstflg=1
                txidx=idx
            if nodedict[node][1]==3:
                idx=len(exondict)+len(splicedict)+len(txstartlist)+txendlist.index(nodedict[node][0])
                row[idx]=weq    
                txendflg=1
                txidx=idx
            
            if node in sparsegraphdict:
                for node2 in sparsegraphdict[node]:
                    edge=(node,node2)
                    if edge in exondict:       
                        idx=exondict[edge]
                        row[idx]=weq
                        outwt+=sparsegraphdict[node][node2][0]
                    elif edge in splicedict:
                        idx=len(exondict)+splicedict[edge]
                        row[idx]=weq
                        outwt+=sparsegraphdict[node][node2][0]  
                    else:
                        message='edge not found %d-%d'%(node,node2)
                        common.printstatus(message,'W',common.func_name(),1)         
                
            if node in  incomingedgedict:
                for node2 in incomingedgedict[node]:
                    edge=(node2,node)
                    if edge in exondict:       
                        idx=exondict[edge]
                        row[idx]=-1*weq
                        inwt+=incomingedgedict[node][node2][0]
                    elif edge in splicedict:
                        idx=len(exondict)+splicedict[edge]
                        row[idx]=-1*weq
                        inwt+=incomingedgedict[node][node2][0]
                    else:
                        message='edge not found %d-%d'%(node,node2)
                        common.printstatus(message,'W',common.func_name(),1)                    
            A.append(row)
            B.append(0)
    
    #       Only internal tx start and end need to closer to difference   
            if txstflg==1:
                txwt=max(0.0,outwt-inwt)
            if txendflg==1:
                txwt=max(0.0,inwt-outwt)
            if txstflg==1 or txendflg==1:
            # and inwt>0.01
                row=[0]*numvar
                row[txidx]=5
                A.append(row)
                B.append(5*txwt)
#                if outwt>0.01:
#                    row=[0]*numvar
#                    row[txidx]=weq/20
#                    A.append(row)
#                    B.append(weq*txwt/20)
#                if outwt<=0.01:
#                    row=[0]*numvar
#                    row[txidx]=weq/5
#                    A.append(row)
#                    B.append(weq*txwt/5)   
            #print node, nodedict[node]
            #print row 
        return((A,B,X))             
          
    def ACT2Corrected(self,gene,num_iterations=5):
        """
            Next steps: Some way to preserve flows at divergence nodes
            One way could be reallocate flows at all divergence nodes in the original ratio and fix it 
            Iterate 10 times
        """
        inwgs=self.wgsdict[gene]
        outwgs=inwgs
        component1=1.0
        for iteri in range(num_iterations):
            component1=1.0-iteri*1.0/num_iterations
            wgs=addwgs(inwgs,outwgs,component1)
            A,B,X=self.wgs2problem(wgs)
            Xvar = cvx.variable(len(X),1)    
            A=cvx.matrix(A)
            B=cvx.matrix(B)
            B=B.T
            p = cvx.program(cvx.minimize(cvx.norm2(A*Xvar-B)),[cvx.geq(Xvar,0.0)])
            try:
                p.solve(quiet=1)
            except:
                message='Could not solve for %s'%(gene)
                common.printstatus(message,'W',common.func_name(),1)   
                return (outwgs,100.0)   
            if iteri==0:                          # Get optimal value
                err=cvx.norm2(A*Xvar-B)
            #print err.value/len(X)
            Xval=Xvar.T.value.tolist()[0]
            X_corr= [a[:] for a in X]
            for i in range(len(Xval)):
                X_corr[i][3]=int(Xval[i]*100)/100.0
            
            #print X_corr
            exonlist=[[a[1],a[2]] for a in X_corr if a[0]==2]
            exonwtlist=[a[3] for a in X_corr if a[0]==2]
            #print 'E',exonlist
            intronlist=[]
            intronwtlist=[]
            splicelist=[[a[1],a[2]] for a in X_corr if a[0]==3]
            splicewtlist=[a[3] for a in X_corr if a[0]==3]
            removelist=[]
            for i in range(len(exonlist)):
                exon=exonlist[i]
                if exon in splicelist:
                    exonwt=exonwtlist[i]
                    intronlist.append([exon[0]+1,exon[1]-1])
                    intronwtlist.append(exonwt)
                    removelist.append(i)
            removelist.reverse()
            for i in removelist:
                exonlist.pop(i)
                exonwtlist.pop(i)
                    
            #print 'E',exonlist
            startnodelist=[a[1]for a in X_corr if a[0]==1]
            endnodelist=[a[1]for a in X_corr if a[0]==-1]
            novelnodelist=wgs[5]
            #print exonlist
            #print wgs[0]
            #print intronlist
            #print wgs[1]

            exonwtlist1=[exonwtlist[i] for i in range(len(exonwtlist)) if exonlist[i] in wgs[0]]
            intronwtlist1=[exonwtlist[i] for i in range(len(exonwtlist)) if exonlist[i] in wgs[1]]
            #wgrstuple=(exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist)
            outwgs=(wgs[0],wgs[1],splicelist,wgs[3],wgs[4],novelnodelist,exonwtlist1,intronwtlist1,splicewtlist)
        
        return (outwgs,err.value/len(X))
        
    def genedivdictTogeneflowdict(self,genedivdict):
        '''
        gene:pos:[flgs],[wts,flowvec]
        '''
        geneflowdict={}
        for gene in genedivdict.keys():
            wgrs=self.getwtgraphstruct(gene)
            if wgrs==[]:
                continue
            wgs=wtgraphstruct(wgrs)
            flowdict=wgs.divdictToflow(genedivdict[gene])
            geneflowdict[gene]=flowdict
        return geneflowdict
    
class wtgraphstruct:
    def __init__(self, tuple):
        self.tuple=tuple
        if len(tuple[6])!=0:
            self.wtflg=1
            
            
    def Toactfilelines(self,island):
        chrnm='%s'%island[0]
        gene=island[3]
        if island[4]=='+':
            tdir=1
        else:
            tdir=0
        actlines=[]
        nodelist=[]
        exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist=self.tuple
        for edge in splicelist:
            if edge[0] not in nodelist:
                nodelist.append(edge[0])
            if edge[1] not in nodelist:
                nodelist.append(edge[1])
        for edge in exonlist+intronlist:
            if edge[0] not in nodelist:
                if edge[0]-1 not in nodelist:
                    nodelist.append(edge[0])
            if edge[1] not in nodelist:
                if edge[1]+1 not in nodelist:
                    nodelist.append(edge[1])
        for node in startnodelist:
            if node not in nodelist:
                nodelist.append(node)
        for node in endnodelist:
            if node not in nodelist:
                nodelist.append(node)
        nodelist.sort()
        for node in nodelist:
            if node in startnodelist:
                actlines.append('%s\tannos\tnode\t%d\t%d\t0.00\t%d\t.\t%s'%(chrnm,node,node,tdir,gene))
            elif node in endnodelist:
                actlines.append('%s\tannoe\tnode\t%d\t%d\t0.00\t%d\t.\t%s'%(chrnm,node,node,tdir,gene))
            elif node in novelnodelist:
                actlines.append('%s\tnovel\tnode\t%d\t%d\t0.00\t%d\t.\t%s'%(chrnm,node,node,tdir,gene))
            else:
                actlines.append('%s\tannot\tnode\t%d\t%d\t0.00\t%d\t.\t%s'%(chrnm,node,node,tdir,gene))
        edgelist=exonlist[0:]+intronlist[0:]+splicelist[0:]
        edgelist.sort()
        for edge in edgelist:
            if edge in exonlist:
                wt=exonwtlist[exonlist.index(edge)]
                actlines.append('%s\tannot\texon\t%d\t%d\t%8.2f\t%d\t.\t%s'%(chrnm,edge[0],edge[1],wt,tdir,gene))
            elif edge in splicelist:
                wt=splicewtlist[splicelist.index(edge)]
                if edge[0] in novelnodelist or edge[1] in novelnodelist:
                    actlines.append('%s\tnovel\tsplice\t%d\t%d\t%8.2f\t%d\t.\t%s'%(chrnm,edge[0],edge[1],wt,tdir,gene))
                else:
                    actlines.append('%s\tannot\tsplice\t%d\t%d\t%8.2f\t%d\t.\t%s'%(chrnm,edge[0],edge[1],wt,tdir,gene))
            elif edge in intronlist:
                wt=intronwtlist[intronlist.index(edge)]
                actlines.append('%s\tnovel\tretint\t%d\t%d\t%8.2f\t%d\t.\t%s'%(chrnm,edge[0],edge[1],wt,tdir,gene))
            else:
                message='Orphan edge: %s'%str(edge)
                common.printstatus(message,'W',common.func_name())   
        return actlines

       
    
    def getedgevalue(self,edge,edgetype):
        '''
        edgetype = 10/11,20/21,3 = exon,retained intron,splice
        foundflg=0 if not found
                =1 if exact match
                =2 if one sided match and +-1 match
                =3 if +-1 match on both
        ALERT: exon wt finding should be improved: All overlapping exon weight should be included
        '''
        goodoverlapsize=40
        exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist=self.tuple
        edge=list(edge)
#        message='Edge: %s; Edgetype: %d'%(str(edge), edgetype)
#        common.printstatus(message,'S',common.func_name())   
#        print zip(exonlist,exonwtlist)
#        print zip(intronlist,intronwtlist)
#        print zip(splicelist,splicewtlist)
#        print startnodelist
#        print endnodelist
#        print novelnodelist   
        allexonlist=exonlist+intronlist
        allexonwtlist=exonwtlist+intronwtlist 
        foundflg=0; wt=0.0
        if edgetype in [10,11,20,21]:
            for i in range(len(allexonlist)):
                exon=allexonlist[i]
                if edge[0]-exon[0] in [-1,0,1] and edge[1]-exon[1] in [-1,0,1] :
                    wt=allexonwtlist[i]
                    if edge[0]-exon[0] in [0] and edge[1]-exon[1] in [0]:
                        foundflg=1
                    elif (edge[0]-exon[0] in [0] and edge[1]-exon[1] in [-1,1]) or (edge[0]-exon[0] in [-1,1] and edge[1]-exon[1] in [0]):
                        foundflg=2
                    else:
                        foundflg=3
                    break
            if foundflg==0:
                if edgetype in [10,20]:
                    for i in range(len(allexonlist)):
                        exon=allexonlist[i]
                        if edge[0]-exon[0] in [-1,0,1]:
                            # and (exon[1]>edge[1] or exon[1]-exon[0]>goodoverlapsize):
                            wt=allexonwtlist[i]
                            foundflg=4
                            break
                if edgetype in [11,21]:
                    for i in range(len(allexonlist)):
                        exon=allexonlist[i]
                        if edge[1]-exon[1] in [-1,0,1]:
                        # and (exon[0]<edge[0] or exon[1]-exon[0]>goodoverlapsize):
                            wt=allexonwtlist[i]
                            foundflg=4
                            break                 
        if edgetype==3:  
            for i in range(len(splicelist)):
                splice=splicelist[i]
                if edge[0]-splice[0] in [0] and edge[1]-splice[1] in [0] :
                    wt=splicewtlist[i]
                    foundflg=1
                    break     
        if foundflg==0:
            message='Edge: %s; Type = %d; Edge weight: %10.4f; Found flag: %d'%(str(edge),edgetype,wt,foundflg) 
            common.printstatus(message,'S',common.func_name()) 
            #message='splicelist is %s'%str(splicelist)
            #common.printstatus(message,'S',common.func_name()) 
        return (wt,foundflg)
    
    def divdictToflow(self,divdict):
        '''
        divdict=position:[[incoming/outgoing=0,1,exonstart=0=no,1=yes,2=start transcript,3=end transcript][exon, flowlist]]
        exonstart=0=no,1=yes,2=start transcript and exonstart,3=end transcript and exonstart,4=start transcript and no exonstart(insplice),
        5=end transcript and no exonstart outsplice 
        flowlist=[exon/splicelist]
        wtgraphstruct=(exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist)
        flowlist len=1 for start or end transcript
        getedgevalue edgetype = 10/11,20/21,3 = exon,intron,splice
        flowdict[position]=[[outgoingflg,exonstartflg],[[wt1,wt2],nflowvec]]
        '''
        flowdict={}
        for position in divdict.keys():
            flowvec=[]
            outgoingflg,exonstartflg=divdict[position][0]
            exon,flowlist=divdict[position][1]
            #message='Flowlist : %s'%str(flowlist)
            #common.printstatus(message,'S',common.func_name())      
            edgetype=10+outgoingflg
            wt1=self.getedgevalue(exon,edgetype)[0]
            if exonstartflg==2:
                if len(flowlist)!=1:
                    message='Transcript Start Exon has incoming splice; Flowlist : %s, %s'%(str(exon),common.fl2str(flowlist))
                    common.printstatus(message,'W',common.func_name())     
                else:
                    #incoming flow
                    edgetype=10+outgoingflg
                    #prev exon
                    wtoth=self.getedgevalue(flowlist[0],edgetype)[0]
                    flowvec=[wtoth,max(wt1-wtoth,0)]
                    if wt1-wtoth<0:
                        message='Flow decreases at transcript start exon: %s; Prev Exon: %s; Weight Start=%10.4f, Before=%10.4f'%(str(exon),common.fl2str(flowlist[0]),wt1,wtoth)
                        common.printstatus(message,'W',common.func_name())     
                    #wt2=wt1
                    wt2=sum(flowvec)
            elif exonstartflg==3:
                if len(flowlist)!=1:
                    message='Transcript End Exon has outgoing splice; Flowlist : %s, %s'%(str(exon),common.fl2str(flowlist))
                    common.printstatus(message,'W',common.func_name())     
                else:
                    #outgoing flow
                    edgetype=10+outgoingflg
                    #next exon
                    wtoth=self.getedgevalue(flowlist[0],edgetype)[0]
                    flowvec=[wtoth,max(wt1-wtoth,0)]
                    if wt1-wtoth<0:
                        message='Flow increases at transcript end exon %s: Prev Exon: %s; Weight End=%10.4f, After=%10.4f'%(str(exon),common.fl2str(flowlist[0]),wt1,wtoth)
                        common.printstatus(message,'W',common.func_name())   
                    #wt2=wt1
                    wt2=sum(flowvec)
            elif exonstartflg==1:
                flowvec=[]
                edgetype=10+outgoingflg
                flowvec.append(self.getedgevalue(flowlist[0],edgetype)[0])
                for flowedge in flowlist[1:]:
                    flowvec.append(self.getedgevalue(flowedge,3)[0])
                wt2=sum(flowvec)
            elif exonstartflg==0:
                flowvec=[]
                for flowedge in flowlist:
                    flowvec.append(self.getedgevalue(flowedge,3)[0])
                wt2=sum(flowvec)  
            if len(flowvec)>0:     
                nflowvec=common.normalize_vector(flowvec)
                flowdict[position]=[[outgoingflg,exonstartflg],[[wt1,wt2],nflowvec]]
        return flowdict

def wgrstuple2file(wgrstuple,island,fileptr):
    wtgrs=wtgraphstruct(wgrstuple)
    actfilelines=wtgrs.Toactfilelines(island)
    for actfileline in actfilelines:
        fileptr.write('%s\n'%actfileline) 
    actlines=[x.strip('\n').split('\t') for x in actfilelines]     

def addwgs(wgs1,wgs2,component1):
    if component1==1.0:
        return wgs1
    exonlist1,intronlist1,splicelist1,startnodelist1,endnodelist1,novelnodelist1,exonwtlist1,intronwtlist1,splicewtlist1=wgs1
    exonlist2,intronlist2,splicelist2,startnodelist2,endnodelist2,novelnodelist2,exonwtlist2,intronwtlist2,splicewtlist2=wgs2
    if (exonlist1==exonlist2) and (intronlist1==intronlist2) and (splicelist1==splicelist2) and (startnodelist1==startnodelist2) and (endnodelist1==endnodelist2) and (novelnodelist1==novelnodelist2):
        exonwtlist=[component1*wt1+(1-component1)*wt2 for wt1,wt2 in zip(exonwtlist1,exonwtlist2)]
        intronwtlist=[component1*wt1+(1-component1)*wt2 for wt1,wt2 in zip(intronwtlist1,intronwtlist2)]
        splicewtlist=[component1*wt1+(1-component1)*wt2 for wt1,wt2 in zip(splicewtlist1,splicewtlist2)]
        return [exonlist1,intronlist1,splicelist1,startnodelist1,endnodelist1,novelnodelist1,exonwtlist,intronwtlist,splicewtlist]
    else:
        message='Two wgs are different \n%s \n%s\n %d,%d,%d,%d,%d,%d'%(str(wgs1),str(wgs2),
                                                                    (exonlist1==exonlist2),(intronlist1==intronlist2),(splicelist1==splicelist2),
                                                                    (startnodelist1==startnodelist2),(endnodelist1==endnodelist2),(novelnodelist1==novelnodelist2))
        common.printstatus(message,'W',common.func_name(),1)

def wgs2sparsegraphdict(wgs):
    '''
    Assumes retained intron edge not = splice edge ever
    '''
    exonlist,intronlist,splicelist,startnodelist,endnodelist,novelnodelist,exonwtlist,intronwtlist,splicewtlist=wgs
    nodelist=[]
    for edge in splicelist:
        if edge[0] not in nodelist:
            nodelist.append(edge[0])
        if edge[1] not in nodelist:
            nodelist.append(edge[1])
    exoniclist=exonlist+intronlist
    exonicwtlist=exonwtlist+intronwtlist
    exonictuplist=zip(exoniclist,exonicwtlist)
    exonictuplist.sort()
    for exonictup in exonictuplist:
        if exonictup[0][0] not in nodelist:
            if exonictup[0][0]-1 not in nodelist:
                nodelist.append(exonictup[0][0])
            else:
                exonictup[0][0]-=1
        if exonictup[0][1] not in nodelist:
            if exonictup[0][1]+1 not in nodelist:
                nodelist.append(exonictup[0][1])
            else:
                exonictup[0][1]+=1     
    nodelist.sort()
    nodedict=dict([(id,[val,0]) for (id,val) in enumerate(nodelist)])
    for nodeid in nodedict:
        if nodedict[nodeid][0] not in novelnodelist:
            nodedict[nodeid][1]=1
        if (nodedict[nodeid][0] in startnodelist) or (nodedict[nodeid][0]-1 in startnodelist) or (nodedict[nodeid][0]+1 in startnodelist):
            nodedict[nodeid][1]=2
        if (nodedict[nodeid][0] in endnodelist) or (nodedict[nodeid][0]-1 in endnodelist) or (nodedict[nodeid][0]+1 in endnodelist):
            nodedict[nodeid][1]=3          
    sparsegraphdict={}
    for exonic,wt in exonictuplist:
        node1=exonic[0]; node2=exonic[1]
        node1idx=nodelist.index(node1)
        node2idx=nodelist.index(node2)
        if  node1idx not in sparsegraphdict:
            sparsegraphdict[node1idx]={}
        sparsegraphdict[node1idx][node2idx]=[wt,1]
    for splice,wt in zip(splicelist,splicewtlist):
        node1=splice[0]; node2=splice[1]
        node1idx=nodelist.index(node1)
        node2idx=nodelist.index(node2)
        if  node1idx not in sparsegraphdict:
            sparsegraphdict[node1idx]={}
            message='new splice node added %d'%(node1)
            common.printstatus(message,'W',common.func_name(),1)
        if node2idx not in sparsegraphdict[node1idx]:
            sparsegraphdict[node1idx][node2idx]=[wt,2]
        else:
            message='Two edges between two nodes %d(%d)-%d(%d)'%(node1,node1idx,node2,node2idx)
            common.printstatus(message,'W',common.func_name(),1)
            sparsegraphdict[node1idx][node2idx][0]+=wt        
    return (sparsegraphdict,nodedict)

