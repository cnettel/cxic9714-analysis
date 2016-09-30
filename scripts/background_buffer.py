#!/usr/bin/env python
"""
A script for analyisis of buffer background (back and front)

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  August 8, 2016
"""

# Import other modules
import numpy as np
import scipy as sp
import h5py, os, sys, time, argparse, logging

# Import modules from src directory
sys.path.append("../src")
import cspad

# Data files
background_filenames  = ["../background/cxic9714-r0174.cxi", "../background/cxic9714-r0175.cxi"]
back_sigma_filename   = "../analysis/detector/back/signal/bg_sigmamap_tmp.h5"
front_sigma_filename  = "../analysis/detector/front/signal/bg_sigmamap_tmp.h5"
back_geometry_filename = "../analysis/detector/back/back_geometry.h5"
front_geometry_filename = "../analysis/detector/front/front_geometry.h5"
#back_mask_filename  = "../analysis/masks/back_mask_halo_reduce.h5"
#front_mask_filename = "../analysis/masks/front_mask.h5"

# CXI keys
BACK_DATA  = 'entry_1/image_1/data'
FRONT_DATA = 'entry_1/image_2/data'

# Parse arguments
parser = argparse.ArgumentParser(prog='background_buffer.py', description='A script for analysis of buffer background')
parser.add_argument('-o', metavar='PATH', type=str, default='./', help='Output path')
args = parser.parse_args()

# CSPAD Detector
cspad_back  = cspad.CSPAD(sigmafile=back_sigma_filename,  geometryfile=back_geometry_filename, 
                          snr=6.,  nx=414,  ny=414, dy=13)
cspad_front = cspad.CSPAD(sigmafile=front_sigma_filename, geometryfile=front_geometry_filename,
                          snr=5.5, nx=1746, ny=1746, pixelsize=110e-6, dx=15, dy=20)

# Initialize counters
back_x   = 0
back_x2  = 0
front_x  = 0
front_x2 = 0
nframes  = 0

# Open CXI file with background from r0174
with h5py.File(background_filenames[0], 'r') as f:

    # Nr. of bacgkround frames
    N = f[BACK_DATA].shape[0]
    
    # Iterate through background events (skipping the first 51)
    for i in range(10500, N):
        sys.stdout.write("\rReading background frame %d/%d" %(i-10500+1,N-10500))
        sys.stdout.flush()

        # Increment back detector photon counts
        bcount = cspad_back.photoncounts(f[BACK_DATA][i])
        back_x  += bcount
        back_x2 += bcount**2

        # Increment front detector photon counts
        fcount = cspad_front.photoncounts(f[FRONT_DATA][i])
        front_x  += fcount
        front_x2 += fcount**2

        # Increment counter
        nframes += 1

# Open CXI file with background from r0175
with h5py.File(background_filenames[1], 'r') as f:

    # Nr. of bacgkround frames
    N = f[BACK_DATA].shape[0]
    
    # Iterate through background events (skipping the first 51)
    for i in range(16000, N):
        sys.stdout.write("\rReading background frame %d/%d" %(i-16000+1,N-16000))
        sys.stdout.flush()

        # Increment back detector photon counts
        bcount = cspad_back.photoncounts(f[BACK_DATA][i])
        back_x  += bcount
        back_x2 += bcount**2

        # Increment front detector photon counts
        fcount = cspad_front.photoncounts(f[FRONT_DATA][i])
        front_x  += fcount
        front_x2 += fcount**2

        # Increment counter
        nframes += 1
        
print "\nCreating background stats"
# Mean/Std of back detector photon counts
back_mean  = back_x / float(nframes)
back_std   = np.sqrt((back_x2 / float(nframes))  - (back_mean**2))

# Mean/Std of back detector photon counts
front_mean = front_x / float(nframes)
front_std  = np.sqrt((front_x2 / float(nframes)) - (front_mean**2))

# Save results to file
with h5py.File(args.o + '/background_buffer_stats.h5', 'w') as f:
    f['back/mean']  = back_mean
    f['back/std']   = back_std
    f['front/mean'] = front_mean
    f['front/std']  = front_std
    f['nframes']    = nframes
