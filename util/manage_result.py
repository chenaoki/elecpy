import glob
import commands

def combine_restart( path_org, path_dst, restart=0, verbose=False, safe=False):
    keys = ['phie', 'vmem', 'cell']
    for key in keys:
        
        data_dst =sorted(glob.glob(path_dst+key+'_*'))
        for i, blob in enumerate(data_dst):
            cmd = 'mv {0} {1}'.format( blob, blob+'.org')
            #if key != 'cell': cmd+='.npy'
            if verbose:
                print cmd
            if not safe:
                commands.getoutput(cmd)

        for i, blob in enumerate(data_dst):
            cmd = 'cp -r {0} {1}'.format( blob+'.org', path_dst+key+'_{0:0>4}'.format(restart+i))
            if key != 'cell': cmd+='.npy'
            if verbose:
                print cmd
            if not safe:
                commands.getoutput(cmd)
                
        for i, blob in enumerate(data_dst):
            cmd = 'rm -r {0}'.format( blob+'.org')
            if verbose:
                print cmd
            if not safe:
                commands.getoutput(cmd)
        
        data_org = sorted(glob.glob(path_org+key+'_*'))
        for i, blob in enumerate(data_org):
            if i <= restart:
                cmd = 'cp -r {0} {1}'.format( blob, path_dst )
                if verbose:
                    print cmd
                if not safe:
                    commands.getoutput(cmd)
