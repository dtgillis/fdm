
do_index='zf_graphs'
do_index='pa_u937A'
do_index='pa_u937A2'
do_index='pa_u937A'
do_index='multitx'
do_index='zf_graphs'
do_index='zf_html'
do_index='zf_graphs'
do_index='test1'
do_index='testfdm2'
do_index='mergebam'
do_index='combp'
do_index='mergebam4'
do_index='confighelper'
do_index='testwgsdict'
do_index='testactcorrect'
do_index='testdom'
do_index='combp2'
do_index='testpval'


if do_index=='zf_graphs':
    import actgraph
    import ast
    import os 
    
    infile='/playpen/rootFDM/project/ZF_GZ0/ffull/TOP_SIG__data2.txt'
    worklist=[]
    for lntxt in open(infile):
        ln=lntxt.strip('\n').split('\t')
        if ln[0]=='gene':
            continue
        l2=ast.literal_eval(ln[15])
        l2.append(ast.literal_eval(ln[14]))
        minx=min([min(x) for x in l2])
        maxx=max([max(x) for x in l2])
        worklist.append([max(float(ln[2]),float(ln[5]),float(ln[8])),ln[0],int(ln[1]),[minx,maxx]])
    
    worklist.sort()
    worklist.reverse()
    worklist=worklist[1900:2000]
    tslist=['ZF_CTRL1_A','ZF_RAP_A','ZF_TG3_A']
    project_name='ZF_GZ0'
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/report/image'%project_name
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        for maxfdm,gene,pos,genomerange in worklist:
            if gene in ['CR623802','OGG1 type 1e']:
                continue
            print maxfdm,gene,pos
            act0.Toimage(gene, imagedir,genomerange,[pos],'pdf')
    
    for maxfdm,gene,pos,genomerange in worklist:
        os.system('montage -geometry 760x1600 -tile 3x1 "%s/*%s*%d-%d.pdf" "%s/all_%s_%s-%s_joined.pdf"'%
                  (imagedir,gene.replace('/','_'),genomerange[0],genomerange[1],imagedir,gene.replace('/','_').replace(' ','_'),genomerange[0],genomerange[1]))

if do_index=='zf_html':
    import actgraph
    import ast
    import os 
    import common
    
    infile='/playpen/rootFDM/project/ZF_GZ0/ffull/TOP_SIG__data2.txt'
    csvfile='/playpen/rootFDM/project/ZF_GZ0/report/TOP_SIG_REPORT.txt'
    fout=open(csvfile,'w')
    htmlfile=open('/playpen/rootFDM/project/ZF_GZ0/report/TOP_SIG_REPORT.html','w')
    outtxtlist=[]
    for lntxt in open(infile):
        ln=lntxt.strip('\n').split('\t')
        if ln[0]=='gene':
            fout.write('%s\t%s\n'%(lntxt.strip(),'act graph'))
            continue
        l2=ast.literal_eval(ln[15])
        l2.append(ast.literal_eval(ln[14]))
        minx=min([min(x) for x in l2])
        maxx=max([max(x) for x in l2])
        imgname='image/all_%s_%d-%d_joined.pdf'%(ln[0],minx,maxx)
        urlname='<a href="%s">act graph</a>'%imgname
        outtxtlist.append([2.0-max(float(ln[2]),float(ln[5]),float(ln[8])),'%s\t%s\n'%(lntxt.strip(),urlname)])
    outtxtlist.sort()
    for maxfdm,txt in outtxtlist:
        fout.write(txt)
    fout.close()
    common.csv2html(open(csvfile),htmlfile)
    htmlfile.close()

    
if do_index=='pa_u937A':
    genelist=['NPM1','RUNX2','HOXA9','GAS6','SF1','SF3B2','SIGLEC6','SMYD3','RAVER1','HNRNPM','BAFF','CEBPA',
              'ELANE','GFI1', 'HOXA9', 'NOTCH2NL','ASXL1','U2AF1']
    genelist.sort()
    tslist=['STEM-U937A']
    project_name='STEM2'
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/image'%project_name
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        for gene in genelist:
            print gene
            act0.Toimage(gene, imagedir)
    
    
if do_index=='pa_u937A2':
    import actgraph    
    tslist=['STEM-U937A']
    project_name='STEM2'
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/image'%project_name
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        act0.Toimage('CHD4',imagedir,[6679249,6716551])    

if do_index=='multitx':
    import cPickle        
    x=cPickle.load(open('/playpen/rootFDM/annotation/fdm_human_hg19_knowngene.gtf/genetranscriptdict.pck'))
    print len(x)
    y=[a for a in x if len(x[a])>1]
    print len(y)
    for z in y:
        print z

if do_index=='test1':
    import common
    ln=['a','b','c','a:1','a:2','b:1','b:2','c:1','c:2','ab:3','ac:3','bc:3']
    hdict=common.findheaderdict(ln,['3'])
    print hdict
    
    
if do_index=='testfdm2':
    import common
    group1=[1,1,1,1,1,1]
    group2=[0,0,0,1,0,0,0,1,1]
    group1=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    group2=[0,1,0,0,0,1,0,0,1,1,0,1]    
    print common.getgroupshufflepvalue(group1,group2,numiter=10000)



    
    
if do_index=='mergebam':
    import os
    
    def mergebamfiles(pathsamtools,mergedict):
        for out in mergedict:
            cmd1='mkdir -p %s'%('/'.join(out.split('/')[:-1]))
            print cmd1
            cmd2='%s/samtools merge %s %s'%(pathsamtools,out,' '.join(mergedict[out]))
            print cmd2
            #os.system(cmd1) 
            #os.system(cmd2)
        
    def twoexptojsdfile(in1,in2,out):
        return
    
    def mergesregenmetadata(mergedict,jsd_pair_dict):
        #new expression.txt
        for out in jsd_pair_dict:
            twoexptojsdfile(jsd_pair_dict[out][0],jsd_pair_dict[out][1],out)
    
    pathsamtools='/home/darshan/nextgen/tools/samtools-0.1.8'
    rootdir='/playpen/sregen'
    outdir='%s/TYPE3'%rootdir
    indirs=['%s/TYPE3A'%rootdir,'%s/TYPE3B'%rootdir,'%s/TYPE3C'%rootdir]
    #the first one dictates the set
    dirs=dict([(x,os.listdir('%s/data'%x)) for x in indirs])
    for x in dirs:
        dirs[x]=['%s/data/%s/alignments.bam'%(x,y) for y in dirs[x]]
        dirs[x].sort()
    #print dirs
    allfiles=[ dirs[x] for x in indirs]
    mergefilelist=zip(*allfiles)
    #print mergefilelist
    out=['%s/data/%s/alignments.bam'%(outdir,x.split('/')[-2]) for x in dirs[indirs[0]]]
    #print out
    mergedict=dict(zip(out,mergefilelist))
    print mergedict
    mergebamfiles(pathsamtools,mergedict)
    
if do_index=='mergebam4':
    import os
    
    def mergebamfiles(pathsamtools,mergedict):
        for out in mergedict:
            cmd1='mkdir -p %s'%('/'.join(out.split('/')[:-1]))
            #print cmd1
            cmd2='%s/samtools merge %s %s'%(pathsamtools,out,' '.join(mergedict[out]))
            print cmd2
            os.system(cmd1) 
            os.system(cmd2)            
        
    def twoexptojsdfile(in1,in2,out):
        return
    
    def mergesregenmetadata(mergedict,jsd_pair_dict):
        #new expression.txt
        for out in jsd_pair_dict:
            twoexptojsdfile(jsd_pair_dict[out][0],jsd_pair_dict[out][1],out)
    
    pathsamtools='/home/darshan/nextgen/tools/samtools-0.1.8'
    rootdir='/playpen/sregen'
    outdir='%s/TYPE4'%rootdir
    indirs=['%s/TYPE4%s'%(rootdir,x) for x in ['A','B','C','D','E','F']]
    #the first one dictates the set
    dirs=dict([(x,os.listdir('%s/data'%x)) for x in indirs])
    #print dirs
    for x in dirs:
        dirs[x]=['%s/data/%s/alignments.bam'%(x,y) for y in dirs[x]]
        dirs[x].sort()
    #print dirs
    allfiles=[ dirs[x] for x in indirs]
    set1=allfiles[0]
    set2=allfiles[1][0:5]+allfiles[2][0:20]+allfiles[1][5:]+allfiles[2][20:]
    set3=allfiles[3][0:5]+allfiles[4][0:10]+allfiles[5][0:10]+allfiles[3][5:]+allfiles[4][10:]+allfiles[5][10:]
    #print len(set1),len(set2),len(set3)
    
    sets=[set1,set2,set3]
    #print sets
    mergefilelist=zip(*sets)
    #print mergefilelist
    out=['%s/data/%s/alignments.bam'%(outdir,x.split('/')[-2]) for x in dirs[indirs[0]]]
    #print out
    mergedict=dict(zip(out,mergefilelist))
    #print mergedict
    mergebamfiles(pathsamtools,mergedict)    

if do_index=='combp':
    import common
    pvaluelist=[0.0001,0.85,0.8,0.85,0.05,0.0002,0.9,0.9]
    print common.combinepvalues(pvaluelist,method='fisher')
    print common.combinepvalues(pvaluelist,method='liptak')
    
if do_index=='confighelper':
    rootdir='/playpen/rootFDM'
    sources=['T4T%02dS01'%x for x in range(1,51)]
    sourcepath=['%s/dataout/%s/alignments.bam'%(rootdir,x) for x in sources]
    origsourcepath=['/playpen/sregen/TYPE4/data/T%02dS01/alignments.bam'%x for x in range(1,51)]
    
    mkdircmds=['mkdir -p %s'%'/'.join(x.split('/')[:-1]) for x in sourcepath]
    print '\n'.join(mkdircmds)
    
    lncmds=['ln -s %s %s'%(x,y) for x,y in zip(origsourcepath,sourcepath)]
    print '\n'.join(lncmds)
    
    project_groups_config_line=','.join(sources)
    data_config_lines='\n'.join(['%s=%s'%(x,y) for x,y in zip(sources,sourcepath) ])
    act_config_lines='\n'.join(['%s=0'%x for x in sources])
    
    print 'project_groups=',project_groups_config_line
    print data_config_lines
    print act_config_lines
    
if do_index=='testwgsdict':
    import common
    import cPickle
    wgsdict=cPickle.load(open('/playpen/rootFDM/dataout/DMx/DM_DM1.wgs'))
    genelist=wgsdict.keys()
    genelist=['AF085995']
    for gene in genelist[0:5]:
        wgs=wgsdict[gene]
        print gene
        print wgs
        sparsedict,nodedict=common.wgs2sparsegraphdict(wgs)
        print sparsedict
        print nodedict
        print

if do_index=='testactcorrect':
    import cPickle
    import actgraph
    import os
    import time
    
    ts='DM2'
    
    rpkfile='/playpen/rootFDM/dataout/%s/DMD_%s.rpk'%(ts,ts)
    actfile='/playpen/rootFDM/dataout/%s/DMD_%s.act'%(ts,ts)
    outactfile='/playpen/rootFDM/dataout/%s/DMDCORR_%s.act'%(ts,ts)
    fileptr=open(outactfile,'w')
    
    rpkdict=dict([lntxt.split('\t')[0],float(lntxt.split('\t')[4])] for lntxt in open(rpkfile))
    
    act0=actgraph.actFile(actfile)
    imagedir1='/playpen/rootFDM/dataout/%s/image'%ts
    imagedir2='/playpen/rootFDM/dataout/%s/imagecorr'%ts
    os.system('mkdir -p %s'%imagedir1)
    os.system('mkdir -p %s'%imagedir2)
    
    printall=0
    
    genelist=['FLJ00075','HES4']
    genelist=['CYP1A1','EGF','CIRBP','DAPK3','DHX15','DNMT1','EIF4G2','EREG','FAT2','CYR61','CXorf2','COX7A1','CPA1',
              'CPB1','CTPS','DCT','ECGF1','EYA3','EXOSC10','CLCN5','COL1A1','CYP1B1','CYP2D6','ERCC3','ERCC5','FAH',
              'F12','DAF','CYP1A2','CYP2B6','CYP2F1','CYP4A11','EIF2AK2','ENTPD1','CTNNB1','DSG3','FGB','CX3CL1',
              'EPO','ELAVL1','ERBB4','POLR2A','EZH2','F13A1','CRH','CHGB','EDN1','DEFB4','DGKG','DLG4','ERN1','CSNK2A2',
              'CSK','EPHA7','EIF2AK3','EPHA3','FGR','DCC','CHRNA4','FGF9','CHRNG','FANCG','CHST1','DUSP3','FKBP5','COX15',
              'CSPG3','DFFA','CRLF1','EGFR','E2F2','DAP','CSF1R','CP','CYP3A5','DYSF','CRYGB','FER','FGF2','DNAJC7','CYP2J2',
              'FADD','CRP','DRD5','CYP3A7','CYP2A13','CYP2C9','FBP2','ELK1','DPYD','FABP4','CUL4A','UBC','CKS2','18S','HPRT1',
              'ABCC8','POLR2A','CNIH','ABCB6','ARMET','CELSR1','ABCB1','ABCC6','ABCG2','CENPF','ABCC5','ABCC4','ALKBH',
              'CACNA2D2','CARD4','ARFGEF1','CACNG3','AP4B1','CD86','ADAT1','BRDG1','CLCA4','C11orf9','APH1A','C14orf122','CGI-119',
              'ASB4','ACTL6B','CLST11240','BET1L','BET1L','PHF21A','BAIAP2','C10orf92','BHMT2','C14orf103','C10orf3','CKAP2','AGPAT5',
              'ChGn','ABHD10','CACNA2D3','BTNL2','ABCC1','AGPAT4','C21orf7','CA10','CABC1','ANKRD2','C11orf15','ANGPTL7','BACH2',
              'ABCG5','CENTD3','APG3L','CDCP1','C20orf116','ASPSCR1','APG9L1','CARD14','ADIPOR2','ALG9','C14orf161','CLUAP1',
              'ACTR5','C9orf76','ASRGL1','C13orf23','ATF6','CDH19','C6orf47','CAP1','ABCC11','CGI-115','C2GNT3','BEX1','CCNL2',
              'ADAMTSL1','C15orf16','ABC1','C16orf33','BM88','AD023','C10orf93','C6orf210','SCCPDH','LMBRD1','C9orf121','C6orf78',
              'C3orf18','C12orf8','AICDA','UBC','CHCHD2','18S','HPRT1','DDX41','POLR2A','GNA13','DELGEF','CRISP3','C1orf61','GNAI3',
              'CYP46A1','GALNT6','EML2','FBXO22','DNAJB9','GLS2','EPN1','GGA1','DKFZP586A0522','GABRP','GAPDHS','DHDH','DLG7','DRD1IP',
              'DCXR','GNG13','FKBP11','DKFZp434H2215','FLJ20130','FAM83E','CASZ1','FLJ20449','EPN3','IARS2','FLJ10374','FLJ10560','DPPA4',
              'FLJ10922','GIMAP5','FOXJ2','FLJ11127','GK001','FSCN3','FGF22','DAB1','FLJ13149','FLJ13220','ELAVL4','GAN','EPB41L4A','ELSPBP1',
              'MOSC1','FLJ13855','GDAP1L1','CUEDC2','DCC1','GLB1L','FYCO1','FLJ12650','EFTUD1','UGT2A3','FLJ13848','CSPP','FLJ12443','FLJ13236',
              'FLJ21019','CCDC15','WDR71','ELL3','FOXB1','DKFZp762E1312','DUSP21','FLJ21901','FLJ14627','FLJ14816','DNAJC4','DCD','CRLF3','DIRAS2',
              'C1orf113','CHMP2B','DNAJB1','CYP3A4','FLJ20014','RNF186','FLJ23584','FBXO39','DEFB125','ESX1L','FLJ37549','FLJ35894','FBXL12',
              'DLEU1','FLJ25168','FXYD7','FLJ20244','UBC','18S','HPRT1','TP53','SLC25A4','TNFRSF9','S100A5','SCN9A','SCP2','SLC5A2','SPIB',
              'CLEC3B','UBE2B','RHAG','SLC12A3','TGFBI','TGM1','TYR','SLC10A2','SOD1','SLC2A4','SLC12A2','TFF1','THBS1','RECQL4','POLR2A',
              'RAD52','TRIM21','RORC','TNFRSF8','TNFSF8','TFRC','SCAP1','SRPK1','TK1','TJP2','TIE1','SHH','TIAM1','TPD52','SRF','TCOF1',
              'TAP1','SERPINB5','USP5','SERPINB7','SNX4','TNFSF14','SNAP23','SGCE','SDHA','SNRPD3','TRIP13','TRFP','TRAF4','SEPT2','SFRS10',
              'SPARCL1','SYNGR1','SFRS11','STOML1','SERPINI1','RPL3L','SFPQ','SLC1A6','SLC17A1','SCO2','THRB','STAT4','TFAP2C','RUNX1T1','TP73',
              'SERPINB2','STAT5A','UGCG','TPST2','SLC9A5','SCGB2A2','SLC22A4','THRSP','TJP1','TNNC2','TNNI2','ROM1','TNP1','SLC16A3','S100A8','SUCLG1',
              'TRA1','TGFBR2','STAT6','SRD5A1','TRPC1','TDG','SFRS3','TMSB4Y','UBC','18S','HPRT1','MYC','RAD51','MDM4','MMP7','NEUROD1','NVL',
              'PEX13','PGM1','PLS1','PMP2','POLB','POLE2','PPP1CC','PRPSAP1','PSD','PSMA5','PSMB10','PTPRB','PZP','MYO5A','PEX1','NQO2','P4HB',
              'PNLIP','PRLR','HTRA1','PPFIA2','POLR2A','MCM3','MGMT','MCM3AP','MLNR','PPP3CA','NFKB2','POMC','PGLYRP1','MVK','NTRK3','PCTK2','PIK3C3','PRKCG',
              'MAPKAPK3','MPP2','PLA2G2A','PTCH','PPIC','MX1','NGFR','MTRF1','PEX14','MGAM','PIGB','MRPL33','NDUFAB1','NDUFS1','PDE1C','MYB','NEF3','POU1F1',
              'NR2C2','MDS1','PPARA','NR1I3','MMP1','MMP13','MMP2','MMP9','MIF','PAX3','MDM2','PER1','PPP6C','NME1','MAGEB3','MYH8','OMG','P2RY11','PITX1','PKP4',
              'NUMB','NOL4','PABPC4','NAP1L3','POLR2G','PSAP','RAB11A','OASL','PCNA','PSMD7','RAC2','PFDN1','NAT2','P2RY1','UBC','18S','HPRT1','CCNA2','CFLAR','ACTN2',
              'ADPRH','ALOX12B','ANXA5','AOX1','CD68','AHSG','AMBP','APOF','ARF5','ATP1B2','BTN1A1','CD58','CDA','CENPA','ACO1','ABCD3','ACADM','ACADS','ADA','APOH',
              'ASPA','BTD','ABCC2','AMT','ANXA1','CACNG1','ANGPT2','CDH1','CCR7','CCL20','ARHGDIB','BLM','POLR2A','CD3D','CD81','ADCYAP1','CD80','ATM','CDK4','CDK8',
              'CDK9','BRS3','APC','APOB','ADCY8','ADM','BARD1','ATP6AP1','ABCB11','BAG1','CDC45L','ALDH4A1','BECN1','B2M','BCL2A1','ANXA13','BPHL','CAD','CALR','BRPF1','CALB1','CCS','CDX4','CASP5','CEACAM5',
              'CDK5R1','ATOH1','ADRB1','ATP1A2','ATP6V1B1','CCR3','CCNB2','CDKN1A','BIRC2','C16orf3','ALPI','C11orf13','CDK5','C4BPB','BAG4','BAT1','ADORA2B',
              'ACAT1','BYSL','BRCA2','ALOX15','BMP3','AP1S1','AMD1','ACADVL','UBC','18S','HPRT1','PLAU','PAX5','POLR2A','RBM5','NRL','RRH','NTSR2','NTS','PTK2','PKIA','NRGN',
              'PDGFRA','OPTN','S100A12','PPIF','PLXNC1','S100P','NELL1','ORC2L','PPP2R5B','RPP30','NPC2','PSMC4','RPIP8','PDIA5','RAB32','NPM3','PWP1','PICALM','POLI','MMD',
              'PADI4','PSD4','MRPL13','ORC6L','RAB30','RB1CC1','RNF10','POLR1D','PADI3','PRRX2','DTL','MASK','NUDT11','NEIL3','RSAD1','C1orf103','PBK','POLR1B','NS3TP1','PELI1',
              'NDUFV3','NGB','MRPL46','RAB38','OSBPL11','MOSPD3','MLPH','PSTPIP2','PRKRIP1','NIP30','RUFY1','PUS3','POU6F2','PPARG','RRAGB','NCK1',
              'REV1L','MS4A4A','OR5V1','RNASE6','P2RY5','PAICS','OR1F1','PURG','PNMA3','RAPGEFL1','RLN3R1','SALL4',
              'NEDD8','PES1','MTCBP-1','NKG7','RAD21','RNUT1','RNASE6','PDIA2','NURIT','PPARD','MPP4','RPP38','NONO','UBC','PREI3','18S','HPRT1','HIF1A','GBAS','INPPL1','GFAP',
              'GFPT1','GLP1R','GYS1','HSD17B2','ITIH2','ITIH3','KCNJ3','KCNS3','KIF3C','LIFR','LOXL2','LUM','FURIN','HGD','ICAM1','G6PD','IL4R','IGFBP2','GABRA1','GRM1',
              'GRM3','GSTP1','IL1R1','ITGB7','JAK3','HOXD13','IGF2','INSL4','POLR2A','GRLF1','LPL','GP1BA','KIT','IL6','IL13','IL18R1','KDR','IGFBP5','GABRA6','GABRD','GRIA2',
              'IGF2R','HBEGF','GSTT1','LGI1','KCNK5','KYNU','HABP2','HADH2','KCNC3','GOLGA5','GIF','GRPR','GZMM','HIP2','FLI1','HCFC1','ICAM3','JAK1',
              'LRP1','JAK2','LTA','HS3ST1','GAS1','KCNA6','GPR65','HOXB7','KCNQ2','GPR8','GPR22','HIST1H1D','HSPA6','JUN','KNG1','JUNB','ITIH4','HRG','IGFBP1','KCNC1',
              'HIST1H3A','FSCN1','HIST3H3','HDAC1','GABRA4','GGH','G6PC','GLI3','KALRN','HBZ','UBC','18S',
              'HPRT1','LMNA','GULP1','HUNK','POLR2A','LOXL1','GPR64','IRS1','MAP3K5','LRRC17','MLLT4','IQGAP2',
              'ITGB1BP2','INPP5A','CD180','KLF10','KIF20A','INADL','KNTC2','GZMA','IRF6','MASP2','HEAB','KIF3A','GUCA2B','MYLPF','HSU79274','LR8','IL1F6','JMJD2A','KIAA0773',
              'KIAA0152','KIAA0101','KIAA0652','MELK','KIAA0196','VASH1','CLCC1','LEPROTL1','HS747E2A','LOC51066','LOC51337','LOC55831','MDS032','HCA112','IL17RB','KRT24','IL1F9',
              'LXN','MEPE','LZTFL1','HB-1','LGR7','ITM2B','TINAGL1','LHX5','MFSD1','March7','HTR5A','MGC4504','MAPKAP1','LENG4','KLHL18','IQCG','HEY2','LHX6','C17orf40','IARS','HOXB8',
              'LOC81691','MEGF11','KCNE1L','BTBD15','TAAR2','HYDIN','C16orf24','LRRC7','MGC2963','MIA2','MAB21L1','KCTD5','HYDIN',
              'HOXA5','KIAA1900','HTRA4','LOC132321','HOXA7','HDHD1A','KIF4A','LIN28','MGC17337','KCNE4','UBC','LDHA','LYPLA2','18S','HPRT1','SPDEF','POLR2A','SEMA3A','WIF1','TRAF2',
              'SLC16A2','TRAF1','SLC17A2','TSPAN31','SIX1','TGOLN2','SLC2A1','SLC22A7','SLC38A3','SOX15','STMN2','VIL1','TU3A','SLC6A14','SLCO2B1','TUSC2','SCRG1','TAGLN3','ST8SIA5','SLC16A8','TMEM5','TTLL3','TTC15','VPS24',
              'TRIM17','TRPC4','VRK3','SPTBN5','UBN1','SPOCK3','SURF2','TSC','SLC30A6','TMEM19','SEC5L1','UEV3','TRIM36','XAB2','SLC15A2','SEZ6L','TNFAIP1','TRIB2','UCP1','TEX27','SRR','SLC13A1','TPK1','SLC30A5','SLC25A23',
              'SHCBP1','TREML2','TMC7','WDR19','SPTBN4','WNT4','STMN4','TCF7L1','TFE3','TBR1','TRPS1','TCF8','TNFAIP3','STAT1','SEPT9','SLC4A1AP','SPATA20','WNT2B','TOE1','SYNC1','ST3GAL4','TLE4',
              'SUGT1','USP28','TMSB10','TEX13B','TBN','SHC3','CD300LG','STAG3','THUMPD1','TAIP-2','TFG','ST18','SMAP','SEMA4G','TXNDC10','TMSL8','STARD8','UBC','18S','HPRT1','RB1','BIRC5','CHAD','CHGA','ABCD1','CD3E','VWF','ZNF262',
              'ZNF258','WRN','POLR2A','ZNF146','AFP','IL1B','VCAM1','FLT3','CGA','CHEK1','MAPK14','TGFA','MAP3K14','MAP3K8','WNT1','RARRES3','XK','ZP2','FPGT','ABCD2','TROAP','ABCA1','DMXL1','MLANA','HMGCS2','COG2','COMMD3',
              'GTF3C5','ZNF225','SEC24D','CRELD1','ZBTB20','GTF2A1','ZNF331','ZNF167','MOSPD1','ZNF350','APOBEC3G','NR0B2','ZNF557','SE57-1','COASY','RARA','ZNF140','ZNF192','ABCA2','P2RXL1',
              'FOXA1','ZFP36L2','ZNF358','FMO5','ZNF101','ZIM3','MMP21','ZSCAN4','ZNF434','ZNF548','WNT2','VAMP3','ZMYM1','UBC','18S','HPRT1','SELE','LDLR','IL8','PTGS2','TYMS','DAD1','FABP1','DPP4','MX2','TNFRSF11B','DNCL1','MAP2K6','AES']
    
    genelist=[]
    allislandlist=cPickle.load(open('/playpen/rootFDM/annotation/fdm_human_hg19_knowngene.gtf/islandlist.pck'))
    allislandlist.sort()
    if len(genelist)>0:
        islandlist=[island for island in allislandlist if island[3] in genelist]
        print len(genelist), len(islandlist)
    else:
        islandlist=allislandlist
    islandlist=islandlist
    lncnt=0
    for island in islandlist:
        lncnt+=1
        gene=island[3]
        time1=time.time()
        wgrstuple,error=act0.ACT2Corrected(gene,5)
        timetaken=time.time()-time1
        meanerr=error/(rpkdict[gene]+0.0001)
        print '###\t%d/%d\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%d'%(lncnt,len(islandlist),gene,error,rpkdict[gene],meanerr,timetaken, 
                                                        len(wgrstuple[0])+len(wgrstuple[1])+len(wgrstuple[2])+len(wgrstuple[3])+len(wgrstuple[4]))
        if lncnt%200==1 or (rpkdict[gene]>10 and meanerr<0.1) or printall==1:  
            #print 'printing...',gene 
            act0.Toimage(gene, imagedir1) 
            act0.Toimage(gene, imagedir2,inwgs=wgrstuple)  
        try:
            actgraph.wgrstuple2file(wgrstuple,island,fileptr)
        except:
            print 'Problem in writing act',[len(x) for x in wgrstuple]
    fileptr.close()

if do_index=='testdom':
    import project
    
    flow_file='/playpen/rootFDM/project/TYPE301/flows/TYPE301_ALL_flows.txt'
    adict=project.extractNdominantsplices(flow_file,N=2)
    k=adict.keys()
    print adict[k[0]]

if do_index=='combp2':
    import common
    pvaluelist=[0.00,0.031200,0.000000,0.0,0.393500,0.000000,0.002200,0.00]
    print common.combinepvalues(pvaluelist,method='fisher')
    print common.combinepvalues(pvaluelist,method='liptak')

if do_index=='testpval':    
    import project
    distrlist=[0,0,0,0,0,0,0,0,0,0]
    x=0
    print project.distr2pvalue(distrlist,x)

    
    
    
    
    
    
    
    
    
    
