import numpy as np
import os
import glob
from optparse import OptionParser

def make_2Darray(src_dir, width, height, out_dir):
    if not os.path.exists(os.path.join(out_dir, 'cell_0000')):
        os.makedirs(os.path.join(out_dir, 'cell_0000'))
    file_dirs = glob.glob(src_dir + '/*.npy')
    for file_dir in file_dirs:
        file_dir = file_dir.replace('/', os.sep)
        file_name = file_dir.split(('/').replace('/', os.sep))[-1]
        src_array = np.load(file_dir)
        out_array = np.ones((height, width)) * src_array
        np.save(os.path.join(out_dir, 'cell_0000', file_name).replace('/', os.sep), out_array)
        if file_name == 'v.npy':
            np.save(os.path.join(out_dir, 'vmem_0000.npy').replace('/', os.sep), out_array)
    phie_array = np.zeros((height, width))
    np.save(os.path.join(out_dir, 'phie_0000.npy').replace('/', os.sep), phie_array)

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

