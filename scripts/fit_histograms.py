#!/usr/bin/env python
"""
A script fitting gaussians to per-pixel histograms.

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  October 28, 2016
"""
# --------------
# IMPORT MODULES
# ----------------
import numpy as np
import scipy as sp
import math
import h5py, os, sys, time, argparse, logging, warnings

# Import modules from src directory
sys.path.append("../src")
import plotting
from fastloop import FastLoop
from fit import fit_photon_histograms

def fitting_mode(args):

    # -------------------------
    # LOAD HISTOGRAMSs and MASK
    # -----------------------------------------
    histfile = h5py.File(args.histfilename, 'r')
    Hmap = histfile['data/histogram']
    Hbinsize = histfile['data/histogramBinsize'][0]
    #Hbinsize = 0.05263157894736842
    Hcount = histfile['data/histogramCount'][0]
    Hmin = histfile['data/histogramMin'][0]
    Hnbins = histfile['data/histogramNbins'][0]
    dark_offset = histfile['data/offset'][:]
    Hbins = np.arange(Hmin, (Hbinsize*(Hnbins-1) + Hmin) + Hbinsize, Hbinsize)
    NY = Hmap.shape[0]
    NX = Hmap.shape[1]
    SH = (NY, NX)
    NXY = NX * NY
    if args.m is None: mask = np.zeros(SH).astype(np.bool)
    else:
        maskfile = h5py.File(args.m, 'r')
        mask = (1 - maskfile['data/data'][:]).astype(np.bool)
        maskfile.close()

    # ----------------------------
    # TEMPORARLY STORE FLAT ARRAYS
    # -----------------------------------------
    infile = h5py.File(args.t + '/tmpin.h5', 'w')
    infile['histogram'] = Hmap[:].reshape(NXY, Hnbins)
    infile['mask'] = mask.flat
    infile.close()
    histfile.close()

    # ---------------------------------------------
    # INITIALIZE TEMPORARY FILE FOR FITTING RESULTS
    # ---------------------------------------------
    outfile = h5py.File(args.t + '/tmpout.h5', 'w')
    for key in ['bg_offset', 'bg_amp', 'bg_sigma', 'photon_offset', 'photon_amp', 'photon_sigma']:
        dset = outfile.create_dataset(key, (NXY,))
        dset.attrs['shape'] = SH
    dset = outfile.create_dataset('status', (NXY,), dtype=h5py.special_dtype(vlen=str))
    dset.attrs['shape'] = SH
    outfile.close()

    # ----------------------
    # PRINT SOME INFORMATION
    # ---------------------------
    start_time = time.localtime()
    timestamp = str(start_time.tm_year) + '%02d' %start_time.tm_mon + '%02d' %start_time.tm_mday + '_' + '%02d' %start_time.tm_hour + '%02d' %start_time.tm_min
    print 'Running a fitting analysis on pixel histograms, started at: ', time.strftime("%a, %d %b %Y %H:%M:%S", start_time)
    print 'Detector shape: ', SH
    print 'Histogram details: %d bins between %d and %d ADUs' %(Hnbins, Hmin, Hmin + (Hnbins-1)*Hbinsize)
    print 'Nr. of bad pixels: %d/%d = %.2f %%' % (mask.sum(), NXY, float(mask.sum()) / NXY * 100.)

    # ---------------
    # START FAST LOOP
    # ---------------------------
    infile = args.t + '/tmpin.h5'
    outfile = args.t + '/tmpout.h5'
    fastloop = FastLoop(infile, outfile, args.c, NXY, fit_photon_histograms, 1000, Hbins)
    fastloop.start()
    fastloop.write()
    os.remove(args.t + '/tmpin.h5')

    # ----------------------------
    # STORE SOME EXTRA INFORMATION
    # ----------------------------
    outfile = h5py.File(args.t + '/tmpout.h5', 'a')
    outfile['dark_offset'] = dark_offset
    outfile.close()

    # --------
    # CLEAN UP
    # ----------------------------
    os.system('cp ' + args.t + '/tmpout.h5 ' +  args.o + '/fitting_results.h5')
    os.remove(args.t + '/tmpout.h5')
    print 'Running a fitting analysis on pixel histograms, finished at: ', time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())


def compare_mode(args):

    # ---------------
    # LOAD HISTOGRAMS
    # -----------------------------------------
    histfile = h5py.File(args.histfilename, 'r')
    Hmap = histfile['data/histogram']
    Hbinsize = histfile['data/histogramBinsize'][0]
    Hcount = histfile['data/histogramCount'][0]
    Hmin = histfile['data/histogramMin'][0]
    Hnbins = histfile['data/histogramNbins'][0]
    dark_offset = histfile['data/offset'][:]
    Hbins = np.arange(Hmin, (Hbinsize*(Hnbins-1) + Hmin) + Hbinsize, Hbinsize)
    NY = Hmap.shape[0]
    NX = Hmap.shape[1]
    SH = (NY, NX)
    NXY = NX * NY

    # --------------------
    # LOAD FITTING RESULTS
    # ----------------------------------------
    fitfile = h5py.File(args.fitfilename, 'r')
    bg_offset = fitfile['bg_offset'][:].reshape(tuple(fitfile['bg_offset'].attrs['shape']))
    bg_amp = fitfile['bg_amp'][:].reshape(tuple(fitfile['bg_amp'].attrs['shape']))
    bg_sigma = fitfile['bg_sigma'][:].reshape(tuple(fitfile['bg_sigma'].attrs['shape']))
    photon_offset = fitfile['photon_offset'][:].reshape(tuple(fitfile['photon_offset'].attrs['shape']))
    photon_amp = fitfile['photon_amp'][:].reshape(tuple(fitfile['photon_amp'].attrs['shape']))
    photon_sigma = fitfile['photon_sigma'][:].reshape(tuple(fitfile['photon_sigma'].attrs['shape']))
    status = fitfile['status'][:].reshape(tuple(fitfile['status'].attrs['shape']))
    dark = fitfile['dark_offset'][:] + bg_offset
    gain = photon_offset - bg_offset

    # ---------
    # LOAD MASK
    # ----------------------------------------------------
    if args.m is None: mask = np.zeros(SH).astype(np.bool)
    else:
        maskfile = h5py.File(args.m, 'r')
        mask = (1 - maskfile['data/data'][:]).astype(np.bool)
        maskfile.close()
        gain[mask] = np.nan

    # -------------
    # PLOT GAIN MAP
    # ------------------------------------------------------
    plot = plotting.Plot(colorbar=True)
    plot.xlabel = 'x (total width = %d pixel)' %gain.shape[1]
    plot.ylabel = 'y (total height = %d pixel)' %gain.shape[0]
    plot.title_label = '%s - gain' %(args.fitfilename)
    plot.colorbar_label = 'ADUs'
    plot.plotting_a_map(0, gain, cmap=args.cmap, vmin=args.vmin, vmax=args.vmax)
    #plot.show()

    # ---------------------
    # PHOTON DETECTION RATE
    # ---------------------
    if args.d:
        N = (photon_amp * np.sqrt(2*math.pi*photon_sigma**2))
        bg_amp /= N
        photon_amp /= N
        Hmap /= N
        pdr = []
        thresholds = np.linspace(0,2,200)
        for i in range(200):
            threshold = bg_offset + thresholds[i] * gain 
            pdr_t = (math.sqrt(math.pi)/2) * (photon_amp*photon_sigma*sp.special.erfc((threshold - photon_offset) / (math.sqrt(2)*photon_sigma))  - bg_amp*bg_sigma * sp.special.erfc( (threshold - bg_offset)  / (math.sqrt(2)*bg_sigma)) )
            pdr.append(pdr_t)

        pdr = np.dstack(pdr)
        bestth =  thresholds[np.argmax(pdr, axis=2)]
        pdrmap = np.amax(pdr,axis=2)
    

    # -------------------------
    # PLOT PHOTON DETECTABILITY 
    # --------------------------
    gaussian_model = lambda x,p: p[0] * np.exp( - (np.power(x-p[1],2)) / (2*np.power(p[2],2)) ) # the gaussian model
    y,x = args.pixel
    H = Hmap[y,x]
    bg_params = [bg_amp[y, x], bg_offset[y, x], bg_sigma[y, x]]
    photon_params = [photon_amp[y,x], photon_offset[y, x], photon_sigma[y, x]]
    photon2_params = [0.08*photon_amp[y,x], 2*photon_offset[y,x], photon_sigma[y, x]]
    bfit = gaussian_model(Hbins, bg_params)
    pfit = gaussian_model(Hbins, photon_params)
    p2fit = gaussian_model(Hbins, photon2_params)
    sfit = 0
    for i in range(5, int(photon_offset[y,x]) - 5 + 1):
        split_params = [0.5*photon_amp[y,x]/(int(photon_offset[y,x]+1)-9), i, bg_sigma[y,x]]
        sfit += gaussian_model(Hbins, split_params)
    sfit2 = 0
    for i in range(5, int(photon_offset[y,x]) - 5 + 1):
        i += photon_offset[y,x]
        split_params = [0.5*0.08*photon_amp[y,x]/(int(photon_offset[y,x]+1)-5), i, bg_sigma[y,x]]
        sfit2 += gaussian_model(Hbins, split_params)
    print photon_params, bg_params
    plot = plotting.Plot(save_png=True)
    plot.axes[0].plot((Hbins-bg_offset[y,x])/gain[y,x], H, 'r-', label='Data')
    plot.axes[0].plot((Hbins-bg_offset[y,x])/gain[y,x], bfit, 'b-', label='Gaussian fit to 0-ph peak')
    plot.axes[0].plot((Hbins-bg_offset[y,x])/gain[y,x], pfit, 'g-', label='Gaussian fit to 1-ph peak')
    plot.axes[0].plot((Hbins-bg_offset[y,x])/gain[y,x], p2fit, 'g-', label='Gaussian fit to 2-ph peak')
    plot.axes[0].plot((Hbins-bg_offset[y,x])/gain[y,x], sfit, 'm-', label='Gaussian fit to 0/1-ph peak')
    plot.axes[0].plot((Hbins-bg_offset[y,x])/gain[y,x], sfit2, 'm-', label='Gaussian fit to 1/2-ph peak')
    plot.axes[0].plot((Hbins-bg_offset[y,x])/gain[y,x], bfit + pfit + p2fit + sfit + sfit2, 'k--', label='Sum of Gaussians')
    if args.d:
        pdryx = pdr[y,x]
        plot.axes[0].plot(thresholds, pdryx, color='0.5', label='Photon detection rate')
        plot.axes[0].axvline(bestth[y,x], color='k', lw=1, ls='dotted')
    #plot.axes[0].plot(Hbins, H-(bfit + pfit), 'm-')
    plot.axes[0].semilogy()
    plot.axes[0].set_ylim([1,None])
    plot.axes[0].set_xlim([-1,4])
    plot.axes[0].set_xlabel("Signal in nr. of photons")
    plot.axes[0].set_title("Pixel: X=%d, Y=%d" %(x,y))
    plot.axes[0].legend()
    plot.show()
    plot.save("fit.png")
    

    # -------------------------
    # PLOT FITTING TO HISTOGRAM
    # ----------------------
    y,x = args.pixel
    H = Hmap[y,x]
    bg_params = [bg_amp[y, x], bg_offset[y, x], bg_sigma[y, x]]
    photon_params = [photon_amp[y, x], photon_offset[y, x], photon_sigma[y, x]]
    gaussian_model = lambda x,p: p[0] * np.exp( - (np.power(x-p[1],2)) / (2*np.power(p[2],2)) ) # the gaussian model
    bfit = gaussian_model(Hbins, bg_params)
    pfit = gaussian_model(Hbins, photon_params)

    print 'Status: ', status[y,x]
    print '0-photon, mean = %.2f, std = %.2f' %(bg_offset[y,x], bg_sigma[y,x])
    print '1-photon, mean = %.2f, std = %.2f' %(photon_offset[y,x], photon_sigma[y,x])
    print 'Estimated gain: ', gain[y,x]
    #print 'Optimal threshold: ', bestth[y,x]
    #print 'Max. photon detectability rate: ', pdrmap[y,x]

    plot = plotting.Plot(save_pdf=True)
    plot.axes[0].plot(Hbins, H, 'r-', label='Data')
    plot.axes[0].plot(Hbins, bfit, 'b-', label='Gaussian fit to 0-ph peak')
    plot.axes[0].plot(Hbins, pfit, 'g-', label='Gaussian fit to 1-ph peak')
    plot.axes[0].plot(Hbins, bfit + pfit, 'k--', label='Sum of Gaussians')
    #plot.axes[0].plot(Hbins, H-(bfit + pfit), 'm-')
    plot.axes[0].semilogy()
    plot.axes[0].set_ylim([1, None])
    plot.axes[0].set_xlabel("Signal in ADUs")
    plot.axes[0].set_title("Pixel: X=%d, Y=%d" %(x,y))
    plot.axes[0].legend()
    plot.show()
    plot.save("test.pdf")

    #print H.shape, Hbins.shape



def generating_mode(args):

    # -----------------------------
    # LOAD FITTING RESULTS AND MASK
    # ----------------------------------------
    fitfile = h5py.File(args.fitfilename, 'r')
    bg_offset = fitfile['bg_offset'][:].reshape(tuple(fitfile['bg_offset'].attrs['shape']))
    bg_amp = fitfile['bg_amp'][:].reshape(tuple(fitfile['bg_amp'].attrs['shape']))
    bg_sigma = fitfile['bg_sigma'][:].reshape(tuple(fitfile['bg_sigma'].attrs['shape']))
    photon_offset = fitfile['photon_offset'][:].reshape(tuple(fitfile['photon_offset'].attrs['shape']))
    photon_amp = fitfile['photon_amp'][:].reshape(tuple(fitfile['photon_amp'].attrs['shape']))
    photon_sigma = fitfile['photon_sigma'][:].reshape(tuple(fitfile['photon_sigma'].attrs['shape']))
    status = fitfile['status'][:].reshape(tuple(fitfile['status'].attrs['shape']))
    dark = fitfile['dark_offset'][:] + bg_offset
    gain = photon_offset - bg_offset
    SH = status.shape
    NXY = SH[0] * SH[1]
    if args.m is None: mask = np.zeros(SH).astype(np.bool)
    else:
        maskfile = h5py.File(args.m, 'r')
        mask = (1 - maskfile['data/data'][:]).astype(np.bool)
        maskfile.close()

    # ----------------------
    # PRINT SOME INFORMATION
    # ------------------------------------------
    mask_bad_pixel = (status == 'is_bad') | mask
    mask_error = (status == 'fit_error') | (status == 'hist_error') | (status == 'other_error')
    mask_ok = (status == 'ok') & (~mask_bad_pixel) & (~mask_error)

    print 'Nr. of bad pixels: %d/%d = %.2f %%' % (mask_bad_pixel.sum(), NXY, float(mask_bad_pixel.sum()) / NXY * 100.)
    print 'Nr. of pixels with errors: %d/%d = %.2f %%' % (mask_error.sum(), NXY, float(mask_error.sum()) / NXY * 100.)
    print 'Nr. of pixels to be used (bad and error excluded): %d/%d = %.2f %%' % (mask_ok.sum(), NXY, float(mask_ok.sum()) / NXY * 100.)
    mask = mask_bad_pixel | mask_error

    # ------------------------------------------------------------
    # CHECK FOR PIXELS WITH TOO WEAK/STRONG PHOTON PEAK AMPLITUDES
    # ------------------------------------------------------------
    if args.pa is not None:
        mask |= (photon_amp < args.pa[0]) | (photon_amp > args.pa[1])
        mask |= np.isnan(photon_amp)
        print 'Nr. of pixels to be used (%.2f < photon amp < %.2f): %d/%d = %.2f %%' % (args.pa[0], args.pa[1], (~mask).sum(), NXY, float((~mask).sum()) / NXY * 100.)

    # ------------------------------------------------------
    # CHECK FOR PIXELS WITH UNREASONABLE PHOTON SIGMA VALUES
    # ------------------------------------------------------
    if args.ps is not None:
        mask |= (photon_sigma < args.ps[0]) | (photon_sigma > args.ps[1])
        print 'Nr. of pixels to be used (%.2f < photon sigma < %.2f): %d/%d = %.2f %%' % (args.ps[0], args.ps[1], (~mask).sum(), NXY, float((~mask).sum()) / NXY * 100.)

    # ----------------------------------------------------------------
    # CHECK FOR PIXELS WITH TOO STRONG/WEAK BACKGROUND PEAK AMPLITUDES
    # ----------------------------------------------------------------
    if args.ba is not None:
        mask |= (bg_amp < args.ba[0]) | (bg_amp > args.ba[1])
        print 'Nr. of pixels to be used (%.2f < bg amp < %.2f): %d/%d = %.2f %%' % (args.ba[0], args.ba[1], (~mask).sum(), NXY, float((~mask).sum()) / NXY * 100.)

    # ----------------------------------------------------------
    # CHECK FOR PIXELS WITH UNREASONABLE BACKGROUND SIGMA VALUES
    # ----------------------------------------------------------
    if args.bs is not None:
        mask |= (bg_sigma < args.bs[0]) | (bg_sigma > args.bs[1])
        print 'Nr. of pixels to be used (%.2f < bg sigma < %.2f): %d/%d = %.2f %%' % (args.bs[0], args.bs[1], (~mask).sum(), NXY, float((~mask).sum()) / NXY * 100.)

    # ----------------------------------------------
    # CHECK FOR PIXELS WITH UNREASONABLE GAIN VALUES
    # ----------------------------------------------
    if args.g is not None:
        mask |= (gain < args.g[0]) | (gain > args.g[1])
        print 'Nr. of pixels to be used (%.2f < gain < %.2f): %d/%d = %.2f %%' % (args.g[0], args.g[1], (~mask).sum(), NXY, float((~mask).sum()) / NXY * 100.)

    # -------------------------------------
    # SHOW HISTOGRAMS OF FITTING PARAMETERS
    # -------------------------------------
    if args.s:# and plotting.plotting_is_enabled:
        params = [photon_amp[~mask], photon_sigma[~mask], bg_amp[~mask],bg_sigma[~mask], gain[~mask]]
        titles = ['Histogram of %s values' %p for p in ['photon amp', 'photon sigma', 'bg amp', 'bg sigma', 'gain']]
        for p in params:
            print p.min(), p.max()
        plotting.plot_fitting_parameter_histograms(params, titles, bins=args.b)

    # ---------------------
    # PHOTON DETECTION RATE
    # ---------------------
    if args.d:
        N = (photon_amp * np.sqrt(2*math.pi*photon_sigma**2))
        bg_amp /= N
        photon_amp /= N
        pdr = []
        thresholds = np.linspace(0,2,200)
        for i in range(200):
            threshold = bg_offset + thresholds[i] * gain 
            pdr_t = (math.sqrt(math.pi)/2) * (photon_amp*photon_sigma*sp.special.erfc((threshold - photon_offset) / (math.sqrt(2)*photon_sigma))  - bg_amp*bg_sigma * sp.special.erfc( (threshold - bg_offset)  / (math.sqrt(2)*bg_sigma)) )
            pdr.append(pdr_t)
        pdr = np.dstack(pdr)
        bestth = thresholds[np.argmax(pdr, axis=2)]
        pdrmap = np.amax(pdr,axis=2)

    # --------------------
    # OVERWRITE BAD PIXELS
    # --------------------
    gain[mask] = 1.
    bg_sigma[mask] = 0.
    gain *= args.r
    if args.i: gain = 1./gain

    # ----------------------------
    # SAVE GAIN MAP AND OTHER MAPS
    # ---------------------------------------------------
    gainfile = h5py.File(args.o + '/gainmap.h5', 'w')
    gainfile['data/data'] = gain
    gainfile.close()
    darkfile = h5py.File(args.o + '/darkcal.h5', 'w')
    darkfile['data/data'] = dark
    darkfile.close()
    bgsigmafile = h5py.File(args.o + '/bg_sigmamap.h5', 'w')
    bgsigmafile['data/data'] = bg_sigma
    bgsigmafile.close()
    phsigmafile = h5py.File(args.o + '/ph_sigmamap.h5', 'w')
    phsigmafile['data/data'] = photon_sigma
    phsigmafile.close()
    if args.d:
        pdrfile = h5py.File(args.o + '/pdrmap.h5', 'w')
        pdrfile['data/data'] = pdrmap
        pdrfile.close()
        bestthfile = h5py.File(args.o + '/bestth.h5', 'w')
        bestthfile['data/data'] = bestth
        bestthfile.close()
    gainmaskfile = h5py.File(args.o + '/included-pixels.h5', 'w')
    gainmaskfile['data/data'] = 1 - mask
    if args.pa is not None: gainmaskfile['data/limits/photon_amp'] = args.pa
    if args.ps is not None: gainmaskfile['data/limits/photon_sigma'] = args.ps
    if args.ba is not None: gainmaskfile['data/limits/bg_amp'] = args.ba
    if args.bs is not None: gainmaskfile['data/limits/bg_sigma'] = args.bs
    if args.g is not None: gainmaskfile['data/limits/gain'] = args.g
    gainmaskfile.close()


# ==========================================================
# ==========================================================
# -------
# LOGGING
# --------------------------
logging.captureWarnings(True)

# ---------------
# PARSE ARGUMENTS
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='fit_histograms.py', description='A script fitting gaussians to per-pixel histograms.')
parser.add_argument('mode', action='store_true', help='fit, compare or generate')
subparsers = parser.add_subparsers()

fit_parser = subparsers.add_parser('fit', help='Fitting Gaussians to the histogram using Scipys leastq')
fit_parser.add_argument('histfilename', metavar='FILE', type=str, help='A histogram file (as created by the cheetah histogram module)')
fit_parser.add_argument('-m', metavar='FILE', type=str, help='A mask file in order to exclude bad pixels from the fitting')
fit_parser.add_argument('-o', metavar='PATH', type=str, default='./', help='Path to the output directory')
fit_parser.add_argument('-t', metavar='PATH', type=str, default='./', help='Path to the directory where things will be stored during runtime')
fit_parser.add_argument('-c', metavar='INT', type=int, default=1, help='Nr. of CPUs to be used')
fit_parser.set_defaults(func=fitting_mode)

compare_parser = subparsers.add_parser('compare', help='Compare the fitted Gaussian to the original histograms')
compare_parser.add_argument('histfilename', metavar='FILE', type=str, help='A histogram file (as created by the cheetah histogram module)')
compare_parser.add_argument('fitfilename', metavar='FILE', type=str, help='A fitting file (as created by the pyGainmap fit module)')
compare_parser.add_argument('pixel', metavar='Y X', type=int, nargs=2, help='y x coordinates to specify which pixel to compare')
compare_parser.add_argument('-m', metavar='FILE', type=str, help='A mask file in order to exclude bad pixels')
compare_parser.add_argument('--cmap', metavar='STR', type=str, default='winter', help='Colormap, default is \'winter\'')
compare_parser.add_argument('--vmin', metavar='FLOAT', type=float, default=None, help='Minimal value')
compare_parser.add_argument('--vmax', metavar='FLOAT', type=float, default=None, help='Maximal value')
compare_parser.add_argument('-d', action='store_true', help='Do photon detectability analysis')
compare_parser.set_defaults(func=compare_mode)

generate_parser = subparsers.add_parser('generate', help = 'Generate a gain map from fitting results')
generate_parser.add_argument('fitfilename', metavar='FILE', type=str, help='A fitting file (as created by the pyGainmap fit module)')
generate_parser.add_argument('-m', metavar='FILE', type=str, help='A mask file in order to exclude bad pixels for the gain map generation')
generate_parser.add_argument('-o', metavar='PATH', type=str, default='./', help='Path to the output directory')
generate_parser.add_argument('-s', action='store_true', help='Show histograms/maps for diagnostic reasons')
generate_parser.add_argument('-b', type=int, metavar='INT', default=100, help='Nr. of bins for showing histograms/maps for diagnostic reasons')
generate_parser.add_argument('-ba', type=float, metavar='FLOAT', nargs=2, help='Minimal and Maximal allowed values for the amplitude of the background peak')
generate_parser.add_argument('-bs', type=float, metavar='FLOAT', nargs=2, help='Minimal and Maximal allowed values for the offset of the background peak')
generate_parser.add_argument('-pa', type=float, metavar='FLOAT', nargs=2, help='Minimal and Maximal allowed values for the amplitude of the photon peak')
generate_parser.add_argument('-ps', type=float, metavar='FLOAT', nargs=2, help='Minimal and Maximal allowed values for the sigma of the photon peak')
generate_parser.add_argument('-g', type=float, metavar='FLOAT', nargs=2, help='Minimal and Maximal allowed values for the gain')
generate_parser.add_argument('-r', type=float, metavar='FLOAT', nargs=1, default=1.,  help='Rescale gain values (important if gainmap results from fluoresence data')
generate_parser.add_argument('-i', action='store_true', help='Save the inverse of the gainmap')
generate_parser.add_argument('-d', action='store_true', help='Do photon detectability analysis')
generate_parser.set_defaults(func=generating_mode)

args = parser.parse_args()
args.func(args)
