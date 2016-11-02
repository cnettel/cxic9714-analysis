#!/usr/bin/env python
"""
A script for phasing a single diffraction pattern and calculate PRTF.

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  August 21, 2016
"""

# Import other modules
import numpy as np
import scipy as sp
import scipy.ndimage as ndimage
import h5py, os, sys, time, argparse, logging

# Path to current directory
curdir = os.path.dirname(os.path.abspath(__file__)) + "/"

# Import modules from src directory
sys.path.append(curdir + "../src")
import cspad

# Import module for phasing
import spimage

# Configuration (Experiment)
distance_back  = 2.4 #[m]
distance_front = 0.497 #[m]
pixelsize  = 110e-6 #[m]
wavelength = 0.2262e-9 #[m]

# Path to file with single strong hit
META = curdir + "../meta/"
filename_signal = META + "single-shot_assembled.h5"

# Parse arguments
parser = argparse.ArgumentParser(prog='phasing.py', description='A script for doing phasing on a single hit.')
parser.add_argument('-o', '--output',  metavar='PATH', type=str, default='./', help='Output path')
args = parser.parse_args()

# Loading the signal (assembled)
with h5py.File(filename_signal, 'r') as f:
    signal_data = np.rot90(f['data/data'][:], k=2)
    signal_mask = ~np.rot90(f['data/mask'][:].astype(np.bool), k=2)  
    signal_diameter = f['diameter'][...]
    signal_intensity = f['intensity'][...]
    signal_id = f['id'][...]

# Downsampling
downsampling = 4
intensities, mask = spimage.downsample(signal_data, downsampling, mask=signal_mask, mode='integrate')

# Cropping
cropping = 154
cx, cy = -2,-6
intensities_cropped = intensities[cropping+cy:-cropping+cy-1,cropping+cx:-cropping+cx-1]
mask_cropped = mask[cropping+cy:-cropping+cy-1,cropping+cx:-cropping+cx-1]

# Phasing parameters
niter_raar = 0
niter_hio  = 5000
niter_er   = 1000
niter_total = niter_raar + niter_hio + niter_er
beta = 0.9
support_size = 24.

# Run phasing with 10000 individual reconstructions
R = spimage.Reconstructor()
R.set_intensities(intensities_cropped)
R.set_mask(mask_cropped)
R.set_number_of_iterations(niter_total)
R.set_number_of_outputs_images(5)
R.set_number_of_outputs_scores(200)
R.set_initial_support(radius=support_size/2)
R.set_support_algorithm("static", number_of_iterations=niter_total)
R.append_phasing_algorithm("raar",beta_init=beta, beta_final=beta, number_of_iterations=niter_raar)
R.append_phasing_algorithm("hio", beta_init=beta, beta_final=beta, number_of_iterations=niter_hio)
R.append_phasing_algorithm("er",  number_of_iterations=niter_er)

# Nr. of reconstructions (N*M total)
N = 100
M = 100

os.system('rm %s' %(args.output + '/phasing.h5'))
for n in range(N):
    output = R.reconstruct_loop(M)
    print "Done Reconstructions: %d/%d" %((n+1)*M, N*M)
    #import matplotlib.pyplot as plt
    #plt.imshow(np.abs(output['real_space_final'][-1]))
    #plt.show()
    with h5py.File(args.output + '/phasing.h5', 'a') as f:
        for k,v in output.iteritems():
            if isinstance(v,dict):
                for kd,vd in v.iteritems():
                    if kd not in f.keys():
                        f.create_dataset(kd, (N*M,), dtype=vd.dtype)
                    f[kd][n*M:(n+1)*M] = vd
            else:
                if k not in f.keys():
                    f.create_dataset(k, (N*M, v.shape[1], v.shape[2]), dtype=v.dtype)
                f[k][n*M:(n+1)*M] = v
