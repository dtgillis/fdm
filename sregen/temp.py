import sys
sys.path.insert(0, '.')

import main
import plot

runid,totreads,readsize='run_22_PE_Large2',5000000,200
runid,totreads,readsize='run_22_PE_Large1',10000000,200

folderdict={'data':[2,2,'/playpen/sregen/%s/data'%runid],'metadata':'/playpen/sregen/%s/metadata'%runid}
main.makemetadataplot(folderdict,totreads,readsize)

print 'done'