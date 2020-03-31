import numpy as np
import os
import glob
import h5py
from optparse import OptionParser

def make_2Darray(src_dir, width, height, out_dir,time):


    with h5py.File(src_dir,"r") as r:
     with h5py.File(os.path.join(out_dir+'out.h5'),'w') as w:
         for param in r['{:0=5}'.format(time)]:
             src_array= r['{0:0=5}/{1}'.format(time,param)]
             out_array=np.ones((height,width))*src_array
             out_array=out_array.flatten()
             w.create_dataset('0000/{}'.format(param),data=out_array)

            
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        '-s','--source_dir',
        dest='src_dir', action='store', type='string', default=None,
        help="source directory name (ex:/Volumes/Recordings/SimulationResults/20170518-2/00001)")
    parser.add_option(
        '-v','--height',
        dest='height', action='store', type='int', default=100,
        help="2D array height")
    parser.add_option(
        '-w','--width',
        dest='width', action='store', type='int', default=100,
        help="2D array width")
    parser.add_option(
        '-o','--out_dir',
        dest='out_dir', action='store', type='string', default=None,
        help="output directory name (ex:/Volumes/Recordings/SimulationResults/20170612-1)")
    (options, args) = parser.parse_args()

    make_2Darray(src_dir=options.src_dir, width=options.width,
                 height=options.height, out_dir=options.out_dir)
