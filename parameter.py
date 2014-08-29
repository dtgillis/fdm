import ConfigParser

import os
import common

class runparam:
    def __init__(self, config_file):
        if os.path.isfile(config_file):
            self.config_file=config_file
        else:
            message='Config file %s does not exist'%config_file
            common.printstatus(message,'F',common.func_name())   
            
    def _processlist(self,tsliststr):
        if len(tsliststr.split('->'))==1:
            tslist=tsliststr.split(',')
            for ts in tslist:
                if ts not in self.data_dict.keys():
                    message='Command Line has error(s): %s does not exist'%ts
                    common.printstatus(message,'F',common.func_name())  
        elif len(tsliststr.split('->'))==2:
            tsrange=tsliststr.split('->')
            alltslist=self.data_dict.keys()
            tslist=[]
            for ts in tslist:
                if ts>=tsrange[0] and ts<=tsrange[1]:
                    tslist.append(ts)
        return tslist
            
    def _checkdir(self,dirname,dirpath,function_name):
        if not(os.path.exists(dirpath)):
            message='Config file %s has errors: '%self.config_file
            message+='%s does not exist'%dirname
            common.printstatus(message,'F',function_name) 
        else:
            return dirpath
        
    def _checkfile(self,filename,filepath,function_name):
        if not(os.path.isfile(filepath)):
            message='Config file %s has errors: '%self.config_file
            message+='%s does not exist'%filename
            common.printstatus(message,'F',function_name) 
        else:
            return filepath       
        
    def _checknum(self,x,type,function_name):
        if type=='int':
            try:
                int(x)
                return int(x)
            except:
                message='Config file %s has errors: '%self.config_file
                message+='%s is not integer'%x
                common.printstatus(message,'F',function_name) 
        if type=='float':
            try:
                float(x)
                return float(x)
            except:
                message='Config file %s has errors: '%self.config_file
                message+='%s is not float'%x
                common.printstatus(message,'F',function_name) 

    def setfdmfastpair(self,fdm_fast_pair):
        allts=fdm_fast_pair.replace('::',',').replace('+',',').split(',')
        if allts[0]!='BLANK':
            for ts in allts:
                if ts not in self.data_dict.keys():
                    message='Command Line has error(s): %s does not exist'%ts
                    common.printstatus(message,'F',common.func_name())  
            run_pair=fdm_fast_pair.split('::')
            self.ffast_run_list=[[run_pair[0].split('+'),run_pair[1].split('+')]]
            self.ffast_prefix='%s__%s__%s'%(self.ffast_prefix,'_'.join(fdm_fast_pair.split('::')[0].split('+')),'_'.join(fdm_fast_pair.split('::')[1].split('+')))
        self.run_all_flag=0
        self.run_pre_act_flag=0        
        self.run_splice_collate_flag=0
        self.run_proj_act_flag=0
        self.run_extract_flows=0
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=1
        self.run_fdm_full_flag=0
        self.run_cluster_flag=0
        self.run_report_flag=0
    
    def setfdmfullpair(self,fdm_full_pair):
        allts=fdm_full_pair.replace('::',',').replace('+',',').split(',')
        if allts[0]!='BLANK':
            for ts in allts:
                if ts not in self.data_dict.keys():
                    message='Command Line has error(s): %s does not exist'%ts
                    common.printstatus(message,'F',common.func_name())  
            run_pair=fdm_full_pair.split('::')
            self.ffull_run_list=[[run_pair[0].split('+'),run_pair[1].split('+')]]
        self.run_all_flag=0
        self.run_pre_act_flag=0        
        self.run_splice_collate_flag=0
        self.run_proj_act_flag=0
        self.run_extract_flows=0
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=0
        self.run_fdm_full_flag=1
        self.run_cluster_flag=0
        self.run_report_flag=0
    
        
    def setfdmgenefile(self,gene_file):
        self.ffull_prefix='%s_%s'%(self.ffull_prefix,gene_file.split('/')[-1].split('.')[0])
        self.ffast_prefix='%s_%s'%(self.ffast_prefix,gene_file.split('/')[-1].split('.')[0])
        self.ffast_genelist=[ln.rstrip('\n').split('\t')[0] for ln in open(gene_file).readlines()]
        self.ffull_genelist=[ln.rstrip('\n').split('\t')[0] for ln in open(gene_file).readlines()]

    def setpreactlist(self,preactlist):
        tslist=self._processlist(preactlist)
        self.pre_act_run_list=tslist
        self.run_all_flag=0
        self.run_pre_act_flag=1        
        self.run_splice_collate_flag=0
        self.run_proj_act_flag=0
        self.run_extract_flows=0
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=0
        self.run_fdm_full_flag=0
        self.run_cluster_flag=0
        self.run_report_flag=0
        
    def setactlist(self,actlist):
        tslist=self._processlist(actlist)
        self.act_run_list=tslist
        self.run_all_flag=0
        self.run_pre_act_flag=0        
        self.run_splice_collate_flag=0
        self.run_proj_act_flag=1
        self.run_extract_flows=0
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=0
        self.run_fdm_full_flag=0
        self.run_cluster_flag=0
        self.run_report_flag=0
                    
            
    def setcollateflg(self,collate_flg):
        self.run_all_flag=0
        self.run_pre_act_flag=0        
        self.run_splice_collate_flag=collate_flg
        self.run_proj_act_flag=0
        self.run_extract_flows=0
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=0
        self.run_fdm_full_flag=0
        self.run_cluster_flag=0
        self.run_report_flag=0

    def setextflowflg(self,flow_chr):
        self.run_all_flag=0
        self.run_pre_act_flag=0        
        self.run_splice_collate_flag=0
        self.run_proj_act_flag=0
        self.run_extract_flows=1
        self.ext_flow_genelist=[flow_chr]
        self.flow_prefix='%s_%s'%(self.flow_prefix,flow_chr.upper())
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=0
        self.run_fdm_full_flag=0
        self.run_cluster_flag=0
        self.run_report_flag=0         
        
    def setextflowgenefile(self,gene_file):
        self.run_all_flag=0
        self.run_pre_act_flag=0        
        self.run_splice_collate_flag=0
        self.run_proj_act_flag=0
        self.run_extract_flows=1
        self.ext_flow_genelist=[ln.rstrip('\n').split('\t')[0] for ln in open(gene_file).readlines()]
        self.flow_prefix='%s_%s'%(self.flow_prefix,gene_file.split('/')[-1].split('.')[0])
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=0
        self.run_fdm_full_flag=0
        self.run_cluster_flag=0
        self.run_report_flag=0         

    def setprojrange(self,proj_range):
        prleft,pright=proj_range.split('-')
        group1=[x for x in self.project_groups[0] if x>=prleft and x<=pright]
        group2=[x for x in self.project_groups[1] if x>=prleft and x<=pright]
        self.project_groups=[group1,group2]
        self.run_all_flag=0
        self.run_pre_act_flag=0        
        self.run_splice_collate_flag=0
        self.run_proj_act_flag=0
        self.run_extract_flows=1
        self.ext_flow_genelist=['all']
        self.flow_prefix='%s_ALL_%s_%s'%(self.flow_prefix,prleft,pright)
        self.run_cluster_flag=0        
        self.run_fdm_fast_flag=0
        self.run_fdm_full_flag=0
        self.run_cluster_flag=0
        self.run_report_flag=0         
            
    def parse(self):
        config = ConfigParser.SafeConfigParser()
        config.read(self.config_file)
        
        self.pathsamtools=self._checkdir('samtools path',config.get('tools','pathsamtools'),common.func_name())
        self.pathucsctools=self._checkdir('ucsc path',config.get('tools','pathucsctools'),common.func_name())
        self.pathbedtools=self._checkdir('bedtools path',config.get('tools','pathbedtools'),common.func_name())
        
        self.chromszfile =self._checkfile('chromosome size file',config.get('reference','chromszfile'),common.func_name())
        self.chrfaifile =self._checkfile('chromosome index file',config.get('reference','chrfaifile'),common.func_name())
        
        self.annotation_file =self._checkfile('gene annotation file',config.get('reference','annotation_file'),common.func_name())
        self.annotation_file_type = config.get('reference','annotation_file_type')
        
        
        self.root_dir=self._checkdir('root_dir path',config.get('project','root_dir'),common.func_name())  
            
        #Data
        data_files=config.items('Data')
        data_dict1=dict(data_files)
        data_dict={}
        for key in data_dict1.keys():
            data_dict[key.upper()]=data_dict1[key]
        for sample in data_dict.keys():
            if not(os.path.isfile(data_dict[sample])):
                message='Config file has error(s): %s does not exist'%data_dict[sample]
                common.printstatus(message,'F',common.func_name())  
        self.data_dict=data_dict
        
        datakeys=data_dict.keys()
        self.datadir_dict=dict(zip(datakeys,['/'.join(data_dict[key].split('/')[:-1]) for key in datakeys]))
        
        #Project
        self.project_name=config.get('project','project_name')
        project_groups=config.get('project','project_groups')
        self.project_groups=[group.split(',') for group in project_groups.split('::')]
        proj_tslist=[]
        for groups in self.project_groups:
            groups.sort()
            proj_tslist+=groups
        for ts in proj_tslist:
            if ts not in datakeys:
                message='Config file has error(s) in [Project] Groups: %s does not exist in Data'%ts
                common.printstatus(message,'F',common.func_name()) 
        project_type=config.get('project','project_type')
        if project_type not in ['1','2']:
            message='Config file has error(s) in Project Type should be 1 or 2'
            common.printstatus(message,'F',common.func_name()) 
        self.project_type=int(project_type)
        self.run_type=self._checknum(config.get('project','run_type'),'int',common.func_name())        

        self.ffast_min_cov=self._checknum(config.get('compute_params','ffast_min_cov'),'float',common.func_name())
        self.ffast_min_fdm=self._checknum(config.get('compute_params','ffast_min_fdm'),'float',common.func_name())
        
        self.ffull_partition=self._checknum(config.get('compute_params','ffull_partition'),'int',common.func_name())
        self.ffull_permutation=self._checknum(config.get('compute_params','ffull_permutation'),'int',common.func_name())
        self.ffull_pvalue=self._checknum(config.get('compute_params','ffull_pvalue'),'float',common.func_name())
        
        self.cluster_max_dbi=self._checknum(config.get('compute_params','cluster_max_dbi'),'float',common.func_name())
        #self.cluster_min_med_cov=self._checknum(config.get('compute_params','cluster_min_med_cov'),'float',common.func_name())
        #self.cluster_min_med_fdm=self._checknum(config.get('compute_params','cluster_min_med_fdm'),'float',common.func_name())
        
        self.ffull_genesplit_size=self._checknum(config.get('compute_params','ffull_genesplit_size'),'int',common.func_name())
        self.report_top_x=self._checknum(config.get('compute_params','report_top_x'),'int',common.func_name())
        self.graph_top_x=self._checknum(config.get('compute_params','graph_top_x'),'int',common.func_name())
        
        
        #Run Flag
        self.run_all_flag=int(config.get('Runflags','run_all_flag'))
        self.run_pre_act_flag=max(self.run_all_flag,int(config.get('Runflags','run_pre_act_flag')))
        self.run_splice_collate_flag=max(self.run_all_flag,int(config.get('Runflags','run_splice_collate_flag')))
        self.run_proj_act_flag=max(self.run_all_flag,int(config.get('Runflags','run_proj_act_flag')))
        self.run_extract_flows=max(self.run_all_flag,int(config.get('Runflags','run_extract_flows')))
        self.run_fdm_fast_flag=max(self.run_all_flag,int(config.get('Runflags','run_fdm_fast_flag')))
        #todo filter
        self.run_fdm_full_flag=max(self.run_all_flag,int(config.get('Runflags','run_fdm_full_flag')))
        # not necessary
        self.run_cluster_flag=max(self.run_all_flag,int(config.get('Runflags','run_cluster_flag')))
        self.run_report_flag=max(self.run_all_flag,int(config.get('Runflags','run_report_flag')))
        
        #Runpreact
        if self.run_pre_act_flag==1:
            pre_act_run_flags_tdict=dict(config.items('Runpreact'))
            pre_act_run_flags_dict={}
            for key in pre_act_run_flags_tdict.keys():
                pre_act_run_flags_dict[key.upper()]=int(pre_act_run_flags_tdict[key])
            for ts in pre_act_run_flags_dict.keys():
                if ts not in datakeys:
                    message='Config file has error(s) in [Runpreact]: %s does not exist in Data'%ts
                    common.printstatus(message,'F',common.func_name())   
            self.pre_act_run_list=[key for key in pre_act_run_flags_dict.keys() if pre_act_run_flags_dict[key]==1]   
            self.pre_act_run_list.sort()
        else:
            self.pre_act_run_list=[]


        #Runact
        if self.run_proj_act_flag==1:
            act_run_flags_tdict=dict(config.items('Runact'))
            act_run_flags_dict={}
            for key in act_run_flags_tdict.keys():
                act_run_flags_dict[key.upper()]=int(act_run_flags_tdict[key])
            for ts in act_run_flags_dict.keys():
                if ts not in datakeys:
                    message='Config file has error(s) in [Runact]: %s does not exist in Data'%ts
                    common.printstatus(message,'F',common.func_name())   
            self.act_run_list=[key for key in act_run_flags_dict.keys() if act_run_flags_dict[key]==1]   
            self.act_run_list.sort()
        else:
            self.act_run_list=[]
        
        #Extractflows
        ext_flow_gene_file=config.get('Extractflows','ext_flow_gene_file')  
        if ext_flow_gene_file.lower() =='all' or ext_flow_gene_file[0:3].lower()=='chr':
            self.ext_flow_genelist=[ext_flow_gene_file.lower()]
            self.flow_prefix=self.project_name
        else:
            self.ext_flow_genelist=[ln.rstrip('\n').split('\t')[0] for ln in open(ext_flow_gene_file).readlines()]
        self.flow_prefix=config.get('Extractflows','flow_prefix') 
        
        #Runfastfdm   
        self.ffast_prefix=config.get('Runfastfdm','ffast_prefix')   
        ffast_gene_file=config.get('Runfastfdm','ffast_gene_file')
        if ffast_gene_file.lower() =='all' or ffast_gene_file[0:3].lower()=='chr':
            self.ffast_genelist=[ffast_gene_file.lower()]
            self.ffast_prefix=self.project_name
        else:
            self.ffast_genelist=[ln.rstrip('\n').split('\t')[0] for ln in open(ffast_gene_file).readlines()]
        
        tffast_run_dict=dict(config.items('Runfastfdm'))
        ffast_run_dict={}
        for key in tffast_run_dict:
            if key[0:9]=='ffast_run':
                ffast_run_dict[key]=tffast_run_dict[key]
        allitems=[]
        for key in ffast_run_dict.keys():
            allitems+=ffast_run_dict[key].replace('::',',').replace('+',',').replace(':',',').replace('|',',').split(',')
        for item in allitems:
            if item not in datakeys:
                message='Config file has error(s): [Runfastfdm] has incorrect definition;[Data] does not have %s'%(item)
                common.printstatus(message,'W',common.func_name()) 
                
        ffast_run_list=[]
        for key in ffast_run_dict.keys():
            run_str=ffast_run_dict[key]
            if len(run_str.split('::'))==1:
                if len(run_str.split(','))>1:
                    # within group a,b,c,d
                    run_items=run_str.split(',')
                    for i in range(len(run_items)-1):
                        for j in range(i+1,len(run_items)):
                            ffast_run_list.append([[run_items[i]],[run_items[j]]])
                else:
                    # all pairs listed a:b|c:d
                    splitrunlist=[x for x in run_str.split('|')]
                    ffast_run_list=[[[x.split(':')[0]],[x.split(':')[1]]] for x in splitrunlist]
            else:
                run_pair=run_str.split('::')
                if len(run_pair[0].split(','))==1:
                    ffast_run_list.append([run_pair[0].split('+'),run_pair[1].split('+')])
                else:
                    for item1 in run_pair[0].split(','):
                        for item2 in run_pair[1].split(','):
                            ffast_run_list.append([[item1],[item2]])
        ffast_run_list.sort()
        self.ffast_run_list=ffast_run_list

        
        #RunCluster     
#        self.cluster_prefix=config.get('RunCluster','cluster_prefix')           
#        cluster_gene_file=config.get('RunCluster','cluster_gene_file')
#        if cluster_gene_file.lower() =='all' or cluster_gene_file[0:3].lower()=='chr':
#            cluster_gene_file=cluster_gene_file.lower()
#            if cluster_gene_file=='chrx':
#                self.cluster_genelist=['chrX']
#            elif cluster_gene_file=='chry':
#                self.cluster_genelist=['chrY']
#            else:
#                self.cluster_genelist=[cluster_gene_file]            
#        else:
#            self.cluster_genelist=[ln.rstrip('\n') for ln in open(cluster_gene_file).readlines()]
        #self.cluster_unknown_flag=int(config.get('RunCluster','cluster_unknown_flag'))

           
        #Runfullfdm
        self.ffull_prefix=config.get('Runfullfdm','ffull_prefix')     
        ffull_gene_file=config.get('Runfullfdm','ffull_gene_file')
        if ffull_gene_file.lower() =='all' or ffull_gene_file[0:3].lower()=='chr':
            self.ffull_genelist=[ffull_gene_file.lower()]
            self.ffull_prefix=self.project_name
        elif ffull_gene_file.lower() =='none':
            self.ffull_genelist=[]
        else:
            self.ffull_genelist=[ln.rstrip('\n').split('\t')[0] for ln in open(ffull_gene_file).readlines()]
        
     
        tffull_run_dict=dict(config.items('Runfullfdm'))
        ffull_run_dict={}
        for key in tffull_run_dict:
            if key[0:9]=='ffull_run':
                ffull_run_dict[key]=tffull_run_dict[key]
        allitems=[]
        for key in ffull_run_dict.keys():
            allitems+=ffull_run_dict[key].replace('::',',').replace('+',',').replace(':',',').replace('|',',').split(',')
        for item in allitems:
            if item not in datakeys:
                message='Config file has error(s): [Runfullfdm] has incorrect definition;[Data] does not have %s'%(item)
                common.printstatus(message,'W',common.func_name()) 
                
        ffull_run_list=[]
        for key in ffull_run_dict.keys():
            run_str=ffull_run_dict[key]
            if len(run_str.split('::'))==1:
                if len(run_str.split(','))>1:
                    # within group
                    run_items=run_str.split(',')
                    for i in range(len(run_items)-1):
                        for j in range(i+1,len(run_items)):
                            ffull_run_list.append([[run_items[i]],[run_items[j]]])
                else:
                    # all pairs listed a:b|c:d
                    splitrunlist=[x for x in run_str.split('|')]
                    ffull_run_list=[[[x.split(':')[0]],[x.split(':')[1]]] for x in splitrunlist]
            else:
                run_pair=run_str.split('::')
                if len(run_pair[0].split(','))==1:
                    ffull_run_list.append([run_pair[0].split('+'),run_pair[1].split('+')])
                else:
                    for item1 in run_pair[0].split(','):
                        for item2 in run_pair[1].split(','):
                            ffull_run_list.append([[item1],[item2]])
        ffull_run_list.sort()
        if self.run_fdm_fast_flag==1:
            self.ffull_run_list=ffull_run_list
        else:
            self.ffull_run_list=[]

        

