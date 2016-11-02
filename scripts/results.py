#!/usr/bin/env python
"""
Extract results from CXI files.

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  November 3, 2016
"""
# --------------
# IMPORT MODULES
# ----------------
import numpy as np
import scipy as sp
import h5py, os, sys, time, argparse, logging
import spimage

# Path to data in CXI file
EXP_ID   = "entry_1/experiment_identifier"
TIME     = "LCLS/machineTime"
TIME_NS  = "LCLS/machineTimeNanoSeconds"
FIDUCIAL = "LCLS/fiducial"
NFRAMES  = "cheetah/event_data/frameNumber"
NPEAKS   = "cheetah/event_data/peakNpix"
BACK_DATA = "entry_1/image_1/data"
BACK_MASK = "entry_1/image_1/mask"
BACK_TAGS = "entry_1/image_1/tags"
BACK_MODEL_DIAMETER = "entry_1/image_1/model/diameterNM"
BACK_MODEL_INTENSITY = "entry_1/image_1/model/intensityMJUM2"
BACK_MODEL_OFFCENTERX = "entry_1/image_1/model/offCenterX"
BACK_MODEL_OFFCENTERY = "entry_1/image_1/model/offCenterY"
ELECTROJET_VOLTAGE = "entry_1/sample_1/injection_1/voltage"
INJECTOR_POSITION = "entry_1/sample_1/geometry_1/translation"
GMD1 = "cheetah/event_data/gmd1"
FRONT_DATA = "entry_1/image_2/data"
FRONT_MASK = "entry_1/image_2/mask"

# Front CSPAD dimension
FRONT_X = (0, 1746)
FRONT_Y = (0, 1746)

# Mask bits
IS_BAD = 128
IS_MISSING = 512

class Data:
    """A data class for extracting data from CXI file"""

    def __init__(self, cspad=False):
        self.expId = []
        self.time  = []
        self.time_ns = []
        self.fiducial = []
        self.backTags = []
        self.nframes = []
        self.npeaks = []
        self.gmd = []
        self.backModelDiameter = []
        self.backModelIntensity = []
        self.backModelOffCenterX = []
        self.backModelOffCenterY = []
        self.electrojetVoltage = []
        self.injectorPosition = []

        if cspad:
            self.backMask = []
            self.backData = []
            self.frontMask = []
            self.frontData = []
        self.cspad = cspad
            
    def add_run(self, file, filter):
        self.expId.append(file[EXP_ID][filter])
        self.time.append(file[TIME][filter])
        self.time_ns.append(file[TIME_NS][filter])
        self.fiducial.append(file[FIDUCIAL][filter])
        self.nframes.append(file[NFRAMES][:].max() - file[NFRAMES][:].min() + 1)
        self.npeaks.append(file[NPEAKS][filter])
        self.gmd.append(file[GMD1][filter])        
        self.backModelDiameter.append(file[BACK_MODEL_DIAMETER][:len(filter)][filter])
        self.backModelIntensity.append(file[BACK_MODEL_INTENSITY][:len(filter)][filter])
        self.backModelOffCenterX.append(file[BACK_MODEL_OFFCENTERX][:len(filter)][filter])
        self.backModelOffCenterY.append(file[BACK_MODEL_OFFCENTERY][:len(filter)][filter])
        self.electrojetVoltage.append(file[ELECTROJET_VOLTAGE][filter])
        self.injectorPosition.append(file[INJECTOR_POSITION][filter,:])

        if self.cspad:
            self.backMask.append(file[BACK_MASK][filter,:,:])
            self.backData.append(file[BACK_DATA][filter,:,:])
            self.frontMask.append(file[FRONT_MASK][filter,FRONT_Y[0]:FRONT_Y[1],FRONT_X[0]:FRONT_X[1]])
            self.frontData.append(file[FRONT_DATA][filter,FRONT_Y[0]:FRONT_Y[1],FRONT_X[0]:FRONT_X[1]])
            
    def expIdToRunNr(self, expid):
        N = expid.shape[0]
        runs = np.hstack([int(expid[i].split('_')[3][2:]) for i in range(N)])
        return runs
            
    def finish(self):
        self.expId = np.hstack(self.expId)
        self.runNr = self.expIdToRunNr(self.expId)
        self.time  = np.hstack(self.time)
        self.time_ns  = np.hstack(self.time_ns)
        self.fiducial = np.hstack(self.fiducial)
        self.nframes = np.sum(self.nframes)
        self.npeaks = np.hstack(self.npeaks)
        self.gmd = np.hstack(self.gmd)
        self.backModelDiameter = np.hstack(self.backModelDiameter)
        self.backModelIntensity = np.hstack(self.backModelIntensity)
        self.backModelOffCenterX = np.hstack(self.backModelOffCenterX)
        self.backModelOffCenterY = np.hstack(self.backModelOffCenterY)
        self.electrojetVoltage = np.hstack(self.electrojetVoltage)
        self.injectorPosition = np.vstack(self.injectorPosition)

        if self.cspad:
            self.backData = np.vstack(self.backData)
            self.backMask = np.vstack(self.backMask)
            self.backMask = ~(((self.backMask & IS_BAD) == IS_BAD) | ((self.backMask & IS_MISSING) == IS_MISSING))
            self.frontData = np.vstack(self.frontData)
            self.frontMask = np.vstack(self.frontMask)
            self.frontMask = ~(((self.frontMask & IS_BAD) == IS_BAD) | ((self.frontMask & IS_MISSING) == IS_MISSING))
            
def advanced_filter(filter_dict, include, exclude):
    filter_include = [np.where(filter_dict.attrs["headings"] == incl)[0][0] if incl in filter_dict.attrs["headings"] else None for incl in include]
    filter_exclude = [np.where(filter_dict.attrs["headings"] == excl)[0][0] if excl in filter_dict.attrs["headings"] else None for excl in exclude]
    filter = np.ones(filter_dict.shape[1]).astype(bool)
    for i in filter_include:
        if i is not None: filter &= filter_dict[i].astype(np.bool)
    for e in filter_exclude:
        if e is not None: filter = filter & ~filter_dict[e].astype(np.bool)
    return filter

def filter_for_selection(expid):
    selected_ids = ['LCLS_2014_Apr14_r0164_183144_f1e0',
                    'LCLS_2014_Apr14_r0170_192837_1b25e',
                    'LCLS_2014_Apr14_r0195_204805_1e32a',
                    'LCLS_2014_Apr15_r0208_095603_1ad3c',
                    'LCLS_2014_Apr15_r0210_100721_166da',
                    'LCLS_2014_Apr15_r0210_101739_cc2d',
                    'LCLS_2014_Apr14_r0165_183436_1e399',
                    'LCLS_2014_Apr14_r0182_201014_168e4']
    return np.array([(expid[:] == id) for id in selected_ids]).any(axis=0)

def extract_results(args):

    # Load data from the files
    d = Data(cspad=args.selection)
    filters = []
    for filename in args.filenames:
        print "Load data from file: ", filename
        f = h5py.File(filename, "r")
        if args.nofilter:
            filter = np.ones(f[EXP_ID].shape[0]).astype(np.bool)
        elif args.selection:
            filter = filter_for_selection(f[EXP_ID])
        else:
            filter = advanced_filter(f[BACK_TAGS], args.include, args.exclude)
        if filter.sum():
            d.add_run(f, filter)
            filters.append(filter)
        f.close()
    d.finish()
    print "Data loading finished: %d frames (of about %d)" %(d.expId.shape[0], d.nframes)
    print "The average hitrate is: %.2f %%" %(float(d.expId.shape[0])/d.nframes*100)
    print "Pulse energy (gmd), mean = %.3f, std = %.3f" %(d.gmd.mean(), d.gmd.std())

    # Convert center position to mrad
    d.backModelOffCenterXmrad = np.arctan(d.backModelOffCenterX * 110e-3 / 2400.) * 1000. + 0.7
    d.backModelOffCenterYmrad = np.arctan(d.backModelOffCenterY * 110e-3 / 2400.) * 1000. + 0.6

    # Post filtering parameters
    postfilter = np.ones(d.expId.shape[0]).astype(np.bool)
    
    # Sort out sizes
    if args.dmin is None: args.dmin = d.backModelDiameter.min()
    if args.dmax is None: args.dmax = d.backModelDiameter.max()
    postfilter &= (d.backModelDiameter >= args.dmin) & (d.backModelDiameter <= args.dmax)

    # Sort out intensities
    if args.imin is None: args.imin = d.backModelIntensity.min()
    if args.imax is None: args.imax = d.backModelIntensity.max()
    postfilter &= (d.backModelIntensity >= args.imin) & (d.backModelIntensity <= args.imax)

    # Sort out injection conditions
    if args.posxmin is None: args.posxmin = d.injectorPosition[:,0].min()
    if args.posxmax is None: args.posxmax = d.injectorPosition[:,0].max()
    if args.poszmin is None: args.poszmin = d.injectorPosition[:,2].min()
    if args.poszmax is None: args.poszmax = d.injectorPosition[:,2].max()
    if args.voltmin is None: args.voltmin = d.electrojetVoltage.min()
    if args.voltmax is None: args.voltmax = d.electrojetVoltage.max()
    postfilter &= (d.injectorPosition[:,0] >= args.posxmin) & (d.injectorPosition[:,0] <= args.posxmax)
    postfilter &= (d.injectorPosition[:,2] >= args.poszmin) & (d.injectorPosition[:,2] <= args.poszmax)
    postfilter &= (d.electrojetVoltage >= args.voltmin) & (d.electrojetVoltage <= args.voltmax)
    
    # Sort out center positions
    if args.centerxmin is None: args.centerxmin = d.backModelOffCenterXmrad.min()
    if args.centerxmax is None: args.centerxmax = d.backModelOffCenterXmrad.max()
    if args.centerymin is None: args.centerymin = d.backModelOffCenterYmrad.min()
    if args.centerymax is None: args.centerymax = d.backModelOffCenterYmrad.max()
    postfilter &= (d.backModelOffCenterXmrad >= args.centerxmin) & (d.backModelOffCenterXmrad <= args.centerxmax)
    postfilter &= (d.backModelOffCenterYmrad >= args.centerymin) & (d.backModelOffCenterYmrad <= args.centerymax)
        
    # Apply post filter
    d.backModelDiameter   = d.backModelDiameter[postfilter]
    d.backModelIntensity  = d.backModelIntensity[postfilter]
    d.backModelOffCenterXmrad = d.backModelOffCenterXmrad[postfilter]
    d.backModelOffCenterYmrad = d.backModelOffCenterYmrad[postfilter]
    d.injectorPosition  = d.injectorPosition[postfilter,:]
    d.electrojetVoltage = d.electrojetVoltage[postfilter]
    d.time  = d.time[postfilter]
    d.time_ns = d.time_ns[postfilter]
    d.runNr = d.runNr[postfilter]
    d.expId = d.expId[postfilter]
    d.npeaks = d.npeaks[postfilter]
    d.gmd = d.gmd[postfilter]
    if args.selection:
        d.backData = d.backData[postfilter]
        d.backMask = d.backMask[postfilter]
    if args.selection:
        d.frontData = d.frontData[postfilter]
        d.frontMask = d.frontMask[postfilter]
    print "Postfilter applied: %d frames" %d.expId.shape[0]
        
    # Intensity in Nr. photons per um^2
    h = 6.62606957e-34 #Js
    c = 299792458 #m/s
    hc = h*c  #Jm
    d.backModelIntensityNrP  = ((d.backModelIntensity / 1000.) * 0.22621e-9) / (hc) 
    d.backModelIntensityNrPh = d.backModelIntensityNrP * np.pi * (1e-3*d.backModelDiameter/2)**2

    # Intensity in W/cm2
    d.backModelIntensityWcm2 = d.backModelIntensity * 1e8 / (1000 * 50e-15)
    
    # Default min/max values for intensity
    if args.imin is None: args.imin = d.backModelIntensity.min()
    if args.imax is None: args.imax = d.backModelIntensity.max()

    # Min/Max/Mean
    print "Diameter ", d.backModelDiameter.min(), d.backModelDiameter.mean(), d.backModelDiameter.max()
    print "Intensity ", d.backModelIntensity.min(), d.backModelIntensity.mean(), d.backModelIntensity.max()
    print "Intensity Wcm2 ", d.backModelIntensityWcm2.min(), d.backModelIntensityWcm2.mean(), d.backModelIntensityWcm2.max()
    print "Intensity NrP ", d.backModelIntensityNrPh.min(), d.backModelIntensityNrPh.mean(), d.backModelIntensityNrPh.max()
    print "Horiz. Deviation ", d.backModelOffCenterXmrad.min(), d.backModelOffCenterXmrad.mean(), d.backModelOffCenterXmrad.max()
    print "Vert. Deviation ", d.backModelOffCenterYmrad.min(), d.backModelOffCenterYmrad.mean(), d.backModelOffCenterYmrad.max()
    print "Injector X ", d.injectorPosition[:,0].min(), d.injectorPosition[:,0].mean(), d.injectorPosition[:,0].max()
    print "Injector Z ", d.injectorPosition[:,2].min(), d.injectorPosition[:,2].mean(), d.injectorPosition[:,2].max()
    print "Injector voltage ", d.electrojetVoltage.min(), d.electrojetVoltage.mean(), d.electrojetVoltage.max()
    
    # Save selected results to results.h5
    if args.results:
        print "Saving results to: ", args.o + 'results.h5'
        with h5py.File(args.o + '/results.h5', 'w') as file:
            file["id"] = d.expId
            file["time_s"] = d.time
            file["time_ns"] = d.time_ns
            file["run_nr"] = d.runNr
            file["diameter"] = d.backModelDiameter
            file["intensity"] = d.backModelIntensity
            #file["intensity_rescaled"] = d.backModelIntensity*d.backNrPhotons/d.backNrPhotonsModel
            file["centerx"] = d.backModelOffCenterXmrad
            file["centery"] = d.backModelOffCenterYmrad
            file["injector_posx"] = d.injectorPosition[:,0]
            file["injector_posz"] = d.injectorPosition[:,2]
            file["injector_voltage"] = d.electrojetVoltage
            file["npeaks"] = d.npeaks
            
    # Save selection of assembled diffraction patterns to selection.h5
    if args.selection:
        print "Save selection to: ", args.o + 'selection.h5'
    
        # Scaling between back and front detector
        distance_back  = 2.4
        distance_front = 0.497
        scaling = distance_back / distance_front
        offset_x = 4
        offset_y = 12
        #print d.backModelOffCenterY / scaling, d.backModelOffCenterX / scaling
        
        # Make a copy of the back/front detector
        assembled = np.round(np.copy(d.frontData) / 17.)
        assembled[assembled<0] = 0.
        backtmp   = np.array([np.rot90(np.copy(d.backData[i]), k=2) for i in range(d.backData.shape[0])])

        # Make a copy of the back/front masks
        assembled_mask = np.copy(d.frontMask)
        backtmp_mask = np.array([np.rot90(np.copy(d.backMask[i]), k=2).astype(np.float) for i in range(d.backMask.shape[0])])
    
        # Rebin the back detector (Zooming out)
        backtmp_list = []
        for i in range(backtmp.shape[0]):
            npixel = np.round((backtmp[i].shape[1]/scaling))
            zoomed = sp.ndimage.zoom(backtmp[i], npixel/backtmp[i].shape[1] * 10., mode='nearest',order=0)
            binned = zoomed.reshape((npixel,10,npixel,10)).sum(axis=(1,3))
            binned = np.round((1./scaling) * binned / 26.)
            binned[binned < 0] = 0.
            backtmp_list.append(binned)
        backtmp = np.array(backtmp_list)
        
        # Rebin the back detector mask (Zooming out)
        backtmp_mask_list = []
        for i in range(backtmp_mask.shape[0]):
            npixel = np.round((backtmp_mask[i].shape[1]/scaling))
            zoomed = sp.ndimage.zoom(backtmp_mask[i], npixel/backtmp_mask[i].shape[1] * 10., mode='nearest',order=0)
            binned = zoomed.reshape((npixel,10,npixel,10)).astype(np.bool).all(axis=(1,3))
            backtmp_mask_list.append(binned)
        backtmp_mask = np.array(backtmp_mask_list)
    
        # Merge rebinned back into front
        N, fh, fw = assembled.shape
        N, bh, bw = backtmp.shape
        for i in range(N):
            assembled[i,fh/2 - bh/2 + offset_y:fh/2 + bh/2 + offset_y, fw/2 - bw/2 + offset_x:fw/2 + bw/2 + offset_x] = backtmp[i]
            assembled_mask[i,fh/2 - bh/2 + offset_y:fh/2 + bh/2 + offset_y, fw/2 - bw/2 + offset_x:fw/2 + bw/2 + offset_x] = backtmp_mask[i]
            
        # Save assembled
        with h5py.File(args.o + '/selection.h5', 'w') as f:
            f['id'] = d.expId
            f['npeaks'] = d.npeaks
            f['time'] = d.time
            f['diameter'] = d.backModelDiameter
            f['intensity'] = d.backModelIntensity
            f['data/data'] = assembled.astype(np.float)
            f['data/data'].attrs['axes'] = np.array(['experiment_identifier:y:x'])
            f['data/mask'] = 1 - assembled_mask
            f['data/mask'].attrs['axes'] = np.array(['experiment_identifier:y:x'])
            
# ---------------
# PARSE ARGUMENTS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='results.py', description='Extract information (correlation plots, histograms, etc...) and good diffraction patterns.')
parser.add_argument('filenames', metavar='FILES', type=str, nargs='+', help='A list of CXI files')
parser.add_argument('--nofilter',action='store_true', help='Do not filter the data.')
parser.add_argument('--filter',  metavar='FILE', type=str, default=None, help='Load filter from file')
parser.add_argument('--include', metavar='TAGS',  type=str, nargs='+', default="all",     help='A list tag to include into filtering of the data (default = all)')
parser.add_argument('--exclude', metavar='TAGS',  type=str, nargs='+', default="nothing", help='A list tag to exclude from filtering of the data (default = nothing)')
parser.add_argument('--dmin', metavar='FLOAT', type=float, default=None, help='Minimal value for diameter')
parser.add_argument('--dmax', metavar='FLOAT', type=float, default=None, help='Maximal value for diameter')
parser.add_argument('--imin', metavar='FLOAT', type=float, default=None, help='Minimal value for intensity')
parser.add_argument('--imax', metavar='FLOAT', type=float, default=None, help='Maximal value for intensity')
parser.add_argument('--posxmin', metavar='FLOAT', type=float, default=None, help='Minimal value for injector position x')
parser.add_argument('--posxmax', metavar='FLOAT', type=float, default=None, help='Maximal value for injector position x')
parser.add_argument('--poszmin', metavar='FLOAT', type=float, default=None, help='Minimal value for injector position z')
parser.add_argument('--poszmax', metavar='FLOAT', type=float, default=None, help='Maximal value for injector position z')
parser.add_argument('--voltmin', metavar='FLOAT', type=float, default=None, help='Minimal value for injector voltage')
parser.add_argument('--voltmax', metavar='FLOAT', type=float, default=None, help='Maximal value for injector voltage')
parser.add_argument('--centerxmin', metavar='FLOAT', type=float, default=None, help='Minimal value for center position in x')
parser.add_argument('--centerxmax', metavar='FLOAT', type=float, default=None, help='Maximal value for center position in x')
parser.add_argument('--centerymin', metavar='FLOAT', type=float, default=None, help='Minimal value for center position in y')
parser.add_argument('--centerymax', metavar='FLOAT', type=float, default=None, help='Maximal value for center position in y')
parser.add_argument('--selection',  action='store_true', help='Save a selection of assembled images')
parser.add_argument('--results',   action='store_true', help='Save results')
parser.add_argument('-o', metavar='PATH', type=str, default='./', help='The output path to store figures')
parser.add_argument('-s', action='store_true', help='Store the plot in a file rather than showing it on screen.')
parser.set_defaults(func=extract_results)
args = parser.parse_args()
args.func(args)
