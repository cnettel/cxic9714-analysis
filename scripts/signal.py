#!/usr/bin/env python
"""
A script for analyisis of a hit (right size) with strong signal (back and front)

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  August 8, 2016
"""

# Import other modules
import numpy as np
import scipy as sp
import scipy.ndimage as ndimage
import h5py, os, sys, time, argparse, logging

# Import modules from src directory
sys.path.append("../src")
import cspad

# Data files
signal_filename  = "../data/hits/r0182_20141209_1845/cxic9714-r0182.cxi"
back_sigma_filename  = "../analysis/detector/back/signal/bg_sigmamap_tmp.h5"
front_sigma_filename = "../analysis/detector/front/signal/bg_sigmamap_tmp.h5"
back_geometry_filename = "../analysis/detector/back/back_geometry.h5"
front_geometry_filename = "../analysis/detector/front/front_geometry.h5"
front_mask_filename = "../analysis/masks/front_mask.h5"

# CXI keys
BACK_DATA  = 'entry_1/image_1/data'
BACK_MASK  = 'entry_1/image_1/mask'
FRONT_DATA = 'entry_1/image_2/data'
FRONT_MASK = "entry_1/image_2/mask" 
DIAMETER   = "entry_1/image_1/model/diameterNM"
INTENSITY  = "entry_1/image_1/model/intensityMJUM2"
OFFCENTERX = "entry_1/image_1/model/offCenterX"
OFFCENTERY = "entry_1/image_1/model/offCenterY"
EXP_ID     = "entry_1/experiment_identifier"
GMD        = "cheetah/event_data/gmd1"

# MASK BYTES
IS_BAD = 128
IS_MISSING = 512 

# Parse arguments
parser = argparse.ArgumentParser(prog='signal.py', description='A script for analysis of a strong hit')
parser.add_argument('-o', metavar='PATH', type=str, default='./', help='Output path')
args = parser.parse_args()

# Loading mask (front)
with h5py.File(front_mask_filename, 'r') as f:
    fmask = f['data/data'][:].astype(np.bool)

# CSPAD Detector
cspad_back  = cspad.CSPAD(sigmafile=back_sigma_filename,  geometryfile=back_geometry_filename, 
                          snr=6.,  nx=414,  ny=414, dy=13)
cspad_front = cspad.CSPAD(sigmafile=front_sigma_filename, geometryfile=front_geometry_filename,
                          snr=5.5, nx=1746, ny=1746, pixelsize=110e-6, dx=15, dy=20)

# Open CXI file with hit
with h5py.File(signal_filename, 'r') as f:

    # Read data and count photons for event with index 46
    back  = cspad_back.photoncounts(f[BACK_DATA][46])
    front = cspad_front.photoncounts(f[FRONT_DATA][46])

    # Read mask
    backmask = f[BACK_MASK][46]
    backmask = ~(((backmask & IS_BAD) == IS_BAD) | ((backmask & IS_MISSING) == IS_MISSING))  
    frontmask = f[FRONT_MASK][46]
    frontmask = ~(((frontmask & IS_BAD) == IS_BAD) | ((frontmask & IS_MISSING) == IS_MISSING))  

    # Read meta information for event 46
    diameter = f[DIAMETER][46]
    intensity = f[INTENSITY][46]
    offcenterx = f[OFFCENTERX][46]
    offcentery = f[OFFCENTERY][46]
    expid = f[EXP_ID][46]
    gmd = f[GMD][46]

print "Loaded event %s with diameter (%.2f nm) and intensity (%.2f mJ/um)" %(expid, diameter, intensity)
print "and center position: (x,y) = ", offcenterx, offcentery
print "and pulse energy (GMD): %.2f mJ", gmd

# Save back and front independently
with h5py.File(args.o + '/signal_back.h5', 'w') as f:
    f['data/data'] = back
    f['data/mask'] = backmask
with h5py.File(args.o + '/signal_front.h5', 'w') as f:
    f['data/data'] = front
    f['data/mask'] = frontmask

# Assemble back and front detector
distance_back  = 2.4
distance_front = 0.5
scaling = distance_back / distance_front
offset_x = 4
offset_y = 12

# Make a copy of the back/front detector
assembled = np.copy(front)
backtmp   = np.rot90(np.copy(back), k=2)

# Make a copy of the back/front masks
assembled_mask = np.copy(frontmask) & fmask & (cspad_front.sigma != -1)
backtmp_mask   = np.rot90(np.copy(backmask), k=2).astype(np.float)

# Rebin the back detector (Zooming out)
npixel = np.round((backtmp.shape[1]/scaling))
zoomed = ndimage.zoom(backtmp, npixel/backtmp.shape[1] * 10., mode='nearest',order=0)
binned = zoomed.reshape((npixel,10,npixel,10)).sum(axis=(1,3))
backtmp = binned / scaling
    
# Rebin the back detector mask (Zooming out)
npixel = np.round((backtmp_mask.shape[1]/scaling))
zoomed = ndimage.zoom(backtmp_mask, npixel/backtmp_mask.shape[1] * 10., mode='nearest',order=0)
binned = zoomed.reshape((npixel,10,npixel,10)).astype(np.bool).all(axis=(1,3))
backtmp_mask = binned

# Merge rebinned back into front
fh, fw = assembled.shape
bh, bw = backtmp.shape
assembled[fh/2 - bh/2 + offset_y:fh/2 + bh/2 + offset_y, fw/2 - bw/2 + offset_x:fw/2 + bw/2 + offset_x] = backtmp
assembled_mask[fh/2 - bh/2 + offset_y:fh/2 + bh/2 + offset_y, fw/2 - bw/2 + offset_x:fw/2 + bw/2 + offset_x] = backtmp_mask

# Save assembled
with h5py.File(args.o + '/signal_assembled.h5', 'w') as f:
    f['id'] = expid
    f['diameter'] = diameter
    f['intensity'] = intensity
    f['data/data'] = assembled.astype(np.float)
    f['data/data'].attrs['axes'] = np.array(['experiment_identifier:y:x'])
    f['data/mask'] = 1 - assembled_mask
    f['data/mask'].attrs['axes'] = np.array(['experiment_identifier:y:x'])
