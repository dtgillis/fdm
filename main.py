import sys
import os

from optparse import OptionParser
import scipy.stats.stats as st
import cvxpy as cvx   
import time

sys.path.insert(0, '..')
import common
import gtffile
import parameter
import gtffile
import bamfile
import junctionfile
import bedgraph
import project
import actgraph
import plot

def testparallelcompletion(run_name,projectdir):
    statusdir='%s/status'%projectdir
    donefilename='%s/%s.done'%(statusdir,run_name)
    os.system('touch %s'%donefilename)
    if run_name in ['p0','p2','p5','p7','p9']:
        numprocs=1
        numdone=1
    time.sleep(5)
    listdir='%s/list'%projectdir
    if run_name[0:2] in ['p1','p3']:
        numprocs=sum(1 for line in open('%s/alllist.txt'%listdir))
        numdone=len([ f for f in os.listdir(statusdir) if f.startswith(run_name[0:2])])
    if run_name[0:2] in ['p4']:
        numprocs=sum(1 for line in open('%s/chrlist.txt'%listdir))
        numdone=len([ f for f in os.listdir(statusdir) if f.startswith(run_name[0:2])])
    if run_name[0:2] in ['p6']:
        numprocs=sum(1 for line in open('%s/allpair.txt'%listdir))
        numdone=len([ f for f in os.listdir(statusdir) if f.startswith(run_name[0:2])])        
    if run_name[0:2] in ['p8']:
        numprocs=sum(1 for line in open('%s/allpair.txt'%listdir))*len([ f for f in os.listdir('%s/splitgene'%listdir)])
        numdone=len([ f for f in os.listdir(statusdir) if f.startswith(run_name[0:2])])    
                
    if numdone==numprocs:
        donefilename='%s/done.%s'%(statusdir,run_name[0:2])
        os.system('touch %s'%donefilename)
                        
    message='%s: %d of %d completed'%(run_name[0:2],numdone,numprocs)
    common.printstatus(message,'S',common.func_name())            
            

if __name__ == "__main__":
    run_name='p0'
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-c', '--config', dest='cfgfile',default='config/fdm_0.cfg')
    parser.add_option('-p', '--preact', dest='preactlist',default='na')
    parser.add_option('-o', '--collate', dest='collate_flg',default='0')    
    parser.add_option('-a', '--act', dest='actlist',default='na')
    parser.add_option('-f', '--extflows', dest='flow_chr',default='chr0')
    parser.add_option('-m', '--catflows',action='store_true', dest='run_concat_flows',default=False)  
    parser.add_option('-g', '--anzgenefile', dest='analysis_genefile',default='')
    parser.add_option('-r', '--projrange', dest='proj_range',default='') 
    parser.add_option('-t', '--fdmfast', dest='fdm_fast_pair',default='') 
    parser.add_option('-w', '--fdmfilter',action='store_true', dest='run_filter_fdm',default=False) 
    parser.add_option('-l', '--fdmfull', dest='fdm_full_pair',default='')    
    parser.add_option('-e', '--report',action='store_true', dest='run_report',default=False)  
    

    # 1: to recreate union junctions, 2 for graph and divdict only
    (options, args) = parser.parse_args()
    config_file=options.cfgfile         
   
    cfg=parameter.runparam(config_file)
    cfg.parse()
    
    if cfg.flow_prefix=='DEFAULT':
        cfg.flow_prefix=cfg.project_name
        
    if cfg.ffast_prefix=='DEFAULT':
        cfg.ffast_prefix=cfg.project_name            

    if cfg.ffull_prefix=='DEFAULT':
        cfg.ffull_prefix=cfg.project_name               
    
    if options.preactlist!='na':
        cfg.setpreactlist(options.preactlist)
        
    if options.collate_flg!='0':
        collate_flg=int(options.collate_flg)
        cfg.setcollateflg(collate_flg)        
        
    if options.actlist!='na':
        cfg.setactlist(options.actlist)
        
    if options.flow_chr!='chr0' and options.flow_chr!='file':
        cfg.setextflowflg(options.flow_chr)
    
    if options.flow_chr=='file' and options.analysis_genefile!='':
        cfg.setextflowgenefile(options.analysis_genefile)
        
    if options.run_concat_flows:
        cfg.run_concat_flows=1
    else:
        cfg.run_concat_flows=0

    if options.run_filter_fdm:
        cfg.run_filter_fdm=1
    else:
        cfg.run_filter_fdm=0        
        
    if options.run_report:
        cfg.run_report_flag=1
    else:
        cfg.run_report_flag=0              
    
    if options.proj_range!='':
        cfg.setprojrange(options.proj_range)
        
    if options.fdm_fast_pair!='':
        cfg.setfdmfastpair(options.fdm_fast_pair)     
        if options.analysis_genefile!='':   
            cfg.setfdmgenefile(options.analysis_genefile)
 
    if options.fdm_full_pair!='':
        cfg.setfdmfullpair(options.fdm_full_pair)     
        if options.analysis_genefile!='':   
            cfg.setfdmgenefile(options.analysis_genefile)
 
    message='Program Started'
    common.printstatus(message,'S',common.func_name())
    
    rootdir=cfg.root_dir
    gtffilename=cfg.annotation_file  

    gtffile=gtffile.gtfFile(gtffilename,0,cfg.annotation_file_type)
    gtffile.getgenetranscriptdict()
   
    geneannodivdict=gtffile.getgeneannodivdict()   
        
    annochrjuncdict=gtffile.getchrjuncdict()
    
    islanddict=gtffile.getislanddict()
    genelist=islanddict.keys()
    chrlist=gtffile.getchrlist()
    islandlist=gtffile.getislandlist()
    chrgenelistdict={}
    for island in islandlist:
        chrnm='%s'%island[0]
        if chrnm not in chrgenelistdict:
            chrgenelistdict[chrnm]=[island[3]]
        else:
            chrgenelistdict[chrnm].append(island[3])
     
    # PRE ACT
    if cfg.run_pre_act_flag==1:
        if len(cfg.pre_act_run_list)==1:
            run_name='p1.%s'%cfg.pre_act_run_list[0]
        else:
            run_name='p1'
        message='Process BAM files for ACT graphs for %s'%str(cfg.pre_act_run_list)
        common.printstatus(message,'S',common.func_name())
        for ts in cfg.pre_act_run_list:
            message='Procesing BAM files for %s'%ts
            common.printstatus(message,'S',common.func_name())
            bamfilename=cfg.data_dict[ts]
            junfilename='%s/%s.jun'%(cfg.datadir_dict[ts],ts)
            bdgfilename='%s/%s.bdg'%(cfg.datadir_dict[ts],ts)
            
            bam0=bamfile.bamFile(bamfilename,[cfg.pathsamtools,cfg.pathbedtools])
            if bam0.checkindex()==0:
                bam0.buildindex()
            bam0.Tobdg(ts)
            bam0.Tojun(ts)
            
            jun0=junctionfile.junFile(junfilename,[cfg.pathucsctools,cfg.chromszfile])
            jun0.Tobigbed(annochrjuncdict)
            jun0.buildindex()
            
            bdg0=bedgraph.bdgFile(bdgfilename,[cfg.pathucsctools,cfg.chromszfile])
            bdg0.Tobigwig() 
            bdg0.buildindex()             
        

    #if max(cfg.run_splice_collate_flag,cfg.run_proj_act_flag,cfg.run_extract_flows,cfg.run_cluster_flag,cfg.run_fdm_fast_flag,cfg.run_fdm_full_flag)==1:
    prj0=project.project(cfg,gtffile)
        
    if cfg.run_splice_collate_flag>=1:
        run_name='p2'
        prj0.setprojmetadict(cfg.run_splice_collate_flag)
        
    # ACT
    if cfg.run_proj_act_flag==1:
        if len(cfg.act_run_list)==1:
            run_name='p3.%s'%cfg.act_run_list[0]
        else:
            run_name='p3'
        message='Creating ACT graphs for %s'%str(cfg.act_run_list)
        common.printstatus(message,'S',common.func_name())
        #islandlist=islandlist[0:30]
        for ts in cfg.act_run_list:
            actfilename=prj0.makeAct(ts, islandlist)
            act0=actgraph.actFile(actfilename)
            act0.buildindex()   
             
    #ext flows
    if cfg.run_extract_flows==1:
        if len(cfg.ext_flow_genelist)==1:
            run_name='p4.%s'%cfg.ext_flow_genelist[0]
        else:
            run_name='p4'        
        message='started Flow computation'
        common.printstatus(message,'S',common.func_name())
        ext_flow_genelist=cfg.ext_flow_genelist
        #print flow_islandlist
        if len(ext_flow_genelist)==1:
            if ext_flow_genelist[0]=='all':
                ext_flow_genelist=genelist
            else:
                if ext_flow_genelist[0] in chrgenelistdict:
                    ext_flow_genelist=chrgenelistdict[ext_flow_genelist[0]]
                else:
                    ext_flow_genelist=[]      
            #elif ext_flow_genelist[0][0:3]=='chr':
            #    ext_flow_genelist=chrgenelistdict[ext_flow_genelist[0]]      
        flow_islandlist=[island for island in islandlist if island[3] in ext_flow_genelist]   
        message='Flow computation number of genes: %d'%len(flow_islandlist)
        common.printstatus(message,'S',common.func_name())
        if len(flow_islandlist)!=0:
            prj0.extractflows(cfg.flow_prefix,flow_islandlist)

    if cfg.run_concat_flows==1:
        run_name='p5'
        message='started Concatenating Flow Files'
        common.printstatus(message,'S',common.func_name())
        prj0.concatchrflows()
        prj0.deletetempflows()

    # FDM Fast Run
    # Create top list on the basis of fast fdm
    if cfg.run_fdm_fast_flag==1:       
        # All pairs
        if cfg.run_type not in [4,5]:
            message='started Fast FDM computation'
            common.printstatus(message,'S',common.func_name())
            
            if len(cfg.ffast_run_list)>1:
                message='This feature is deprecated. Only run one pair at a time'
                common.printstatus(message,'F',common.func_name())
            ffast_tspair=cfg.ffast_run_list[0]
            run_name='p6.%s__%s'%('_'.join(ffast_tspair[0]),'_'.join(ffast_tspair[1]))  
    #        print cfg.ffast_run_list  
    #        print ffast_tspair    
    #        print run_name
    #        print cfg.ffast_prefix
            ffast_genelist=cfg.ffast_genelist
            #print flow_islandlist
            if len(ffast_genelist)==1:
                if ffast_genelist[0]=='all':
                    ffast_genelist=genelist
                elif ffast_genelist[0][0:3]=='chr':
                    ffast_genelist=chrgenelistdict[ffast_genelist[0]]      
            ffast_islandlist=[island for island in islandlist if island[3] in ffast_genelist]   
            message='Fast FDM computation number of genes: %d'%len(ffast_islandlist)
            common.printstatus(message,'S',common.func_name())      
            #prj0.compute_fast_fdm(cfg.ffast_min_cov,cfg.ffast_min_fdm,cfg.flow_prefix,cfg.ffast_prefix,cfg.ffast_genelist,ffast_tspairlist)
            prj0.compute_fast_fdm(cfg.flow_prefix,cfg.ffast_prefix,ffast_tspair)
        else:
            run_name='p6.blank'
    
    if cfg.run_filter_fdm==1:
        #merge all fdm files
        #filter
        #split gene files
        run_name='p7'
        if cfg.run_type not in [4,5]:
            message='started FDM filter'
            common.printstatus(message,'S',common.func_name()) 
            prj0.merge_fast_fdm()       
            prj0.filter_fast_fdm()
        prj0.splitgenelist()  

    #FDM Full Run  
    if cfg.run_fdm_full_flag==1:
        if cfg.run_type not in [4,5]:
            if len(cfg.ffull_run_list)==1:
                run_name='p8.%s__%s.%s'%('_'.join(cfg.ffull_run_list[0][0]),'_'.join(cfg.ffull_run_list[0][1]),options.analysis_genefile.split('/')[-1])
            else:
                run_name='p8'            
            #print run_name
            message='Compute FULL FDM for %s'%str(cfg.ffull_run_list)
            common.printstatus(message,'S',common.func_name())
            ffull_genelist=cfg.ffull_genelist
            if len(ffull_genelist)==1:
                if ffull_genelist[0]=='all':
                    ffull_genelist=genelist
                elif ffull_genelist[0][0:3]=='chr':
                    ffull_genelist=chrgenelistdict[ffull_genelist[0]]
            for tspair in cfg.ffull_run_list:
                message='started fdm run for %s'%str(tspair)
                common.printstatus(message,'S',common.func_name())
                prj0.compute_full_fdm(tspair,ffull_genelist)
        else:
            run_name='p8.blank'   

    
    if cfg.run_report_flag==1:
        run_name='p9'
        prj0.run_report()
            


    projectdir='%s/project/%s'%(cfg.root_dir,cfg.project_name)
    testparallelcompletion(run_name,projectdir)
    message='Done'
    common.printstatus(message,'S',common.func_name())
     
