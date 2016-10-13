#!/usr/bin/env python
"""
A script for assembling the background from back and front

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  September 23, 2016
"""

# Import other modules
import numpy as np
import scipy as sp
import scipy.ndimage as ndimage
import h5py, os, sys, time, argparse, logging

# Import modules from src directory
sys.path.append("../src")
import cspad

# Parse arguments
parser = argparse.ArgumentParser(prog='background.py', description='A script for assembling background')
parser.add_argument('-o', metavar='PATH', type=str, default='./', help='Output path')
args = parser.parse_args()

# Load mean background from back and front
with h5py.File('../analysis/background/background_buffer_stats.h5', 'r') as f:
    back = f['back/mean'][:]
    front = f['front/mean'][:]

# Load assembled mask
with h5py.File('../analysis/signal/signal_assembled.h5', 'r') as f:
    assembled_mask = f['data/mask'][:]
    
# Assemble back and front detector
distance_back  = 2.4
distance_front = 0.5
scaling = distance_back / distance_front
offset_x = 4
offset_y = 12

# Make a copy of the back/front detector
assembled = np.copy(front).astype(np.float)
backtmp   = np.rot90(np.copy(back), k=2)

print assembled

# Make a copy of the back/front masks
#assembled_mask = np.copy(frontmask) & fmask & (cspad_front.sigma != -1)
#backtmp_mask   = np.rot90(np.copy(backmask), k=2).astype(np.float)

# Rebin the back detector (Zooming out)
npixel = np.round((backtmp.shape[1]/scaling))
zoomed = ndimage.zoom(backtmp, npixel/backtmp.shape[1] * 10., mode='nearest',order=0)
binned = zoomed.reshape((npixel,10,npixel,10)).sum(axis=(1,3))
backtmp = binned / 100. * (scaling**2)
    
# Rebin the back detector mask (Zooming out)
#npixel = np.round((backtmp_mask.shape[1]/scaling))
#zoomed = ndimage.zoom(backtmp_mask, npixel/backtmp_mask.shape[1] * 10., mode='nearest',order=0)
#binned = zoomed.reshape((npixel,10,npixel,10)).astype(np.bool).all(axis=(1,3))
#backtmp_mask = binned

# Merge rebinned back into front
fh, fw = assembled.shape
bh, bw = backtmp.shape
assembled[fh/2 - bh/2 + offset_y:fh/2 + bh/2 + offset_y, fw/2 - bw/2 + offset_x:fw/2 + bw/2 + offset_x] = backtmp
#assembled_mask[fh/2 - bh/2 + offset_y:fh/2 + bh/2 + offset_y, fw/2 - bw/2 + offset_x:fw/2 + bw/2 + offset_x] = backtmp_mask

# Save assembled
with h5py.File(args.o + '/background_assembled.h5', 'w') as f:
    f['data/data'] = assembled.astype(np.float)
    f['data/data'].attrs['axes'] = np.array(['experiment_identifier:y:x'])
    f['data/mask'] = assembled_mask
    f['data/mask'].attrs['axes'] = np.array(['experiment_identifier:y:x'])
