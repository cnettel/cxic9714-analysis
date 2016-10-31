#!/usr/bin/env python
"""
A script for merging pixel histograms from different runs.

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  October 28, 2016
"""
# --------------
# IMPORT MODULES
# ----------------
import numpy as np
import scipy as sp
import h5py, os, sys, time, argparse, logging

def merge_histograms(args):

    # Extract hsitogram details from original file
    with h5py.File(args.origfile, 'r') as file:
        Min = file['data/histogramMin'][0]
        Binsize = file['data/histogramBinsize'][0]
        Count = file['data/histogramCount'][0]
        Nbins = file['data/histogramNbins'][0]
        offset = file['data/offset'][:]
        NY, NX = file['data/histogram'].shape[:2]
        N = NY * NX

    # Go through histograms and merge them
    histfiles = args.mergefiles
    Histograms = np.zeros((NY, NX, Nbins))
    TotalCount = 0
    for histfile in histfiles:

        # Load flat Histogram
        print 'Loading histogram file: %s' %histfile
        with h5py.File(histfile, 'r') as file:
            Hflat = file['data/histogram'][:].reshape(N, Nbins)
            TotalCount += file['data/histogramCount'][0]

        if args.cm:
            # Correct for common mode
            Hflat[:,-1] = 0
            Hflat[:,0]  = 0
            cutoff = args.c - Min
            cm = (Hflat[:,:cutoff].argmax(axis=1) + Min).reshape(N,1).repeat(Nbins, axis=1)

            H = np.array(range(Nbins)).reshape(Nbins, 1).repeat(N, axis=1)
            coord = np.array(range(N)).reshape(N,1).repeat(Nbins, axis=1)

            Htmp = H.T + cm
            outliers = (Htmp >= Nbins) | (Htmp < 0)
            H = Htmp % Nbins

            Hcorr = Hflat[coord.flatten(),H.flatten()].reshape(NY*NX, Nbins)
            Hcorr[outliers] = 0

            # Merge histogram
            Histograms += Hcorr.reshape(NY, NX, Nbins)
        
        else:
            Histograms += Hflat.reshape(NY,NX, Nbins)


    # Save merged histogram
    with h5py.File(args.o, 'w') as file:
        file.create_group('data')
        file['data/histogram'] = Histograms
        file['data/histogramBinsize'] = np.array([Binsize])
        file['data/histogramCount'] = np.array([TotalCount])
        file['data/histogramMin'] = np.array([Min])
        file['data/histogramNbins'] = np.array([Nbins])
        file['data/offset'] = offset
        file['data/data'] = h5py.SoftLink('data/histogram')
        print 'Saved merged histogram: %s' %file.filename

# ---------------
# PARSE ARGUMENTS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='merge_histograms.py', description='A script merging pixel histograms.')
parser.add_argument('origfile', metavar='FILE', type=str, help='A histogram file (from cheetah histogram module)')
parser.add_argument('mergefiles', metavar='FILES', type=str, nargs='+', help='A bunch of histogram files (from cheetah histogram module)')
parser.add_argument('--cm', metavar='BOOL', type=bool, default=False, help='Common mode correction befor merging')
parser.add_argument('-o', metavar='PATH', type=str, default='merged_histogram.h5', help='Output path to store merge histogram')
parser.add_argument('-c', metavar='INT', type=int, default=10, help='Cutoff given in ADUs right of the histograms origin (for estimate of common mode)')
parser.set_defaults(func=merge_histograms)
args = parser.parse_args()
args.func(args)
