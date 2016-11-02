"""
A module fitting gaussian's to a histogram

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last update:  March 5, 2014
"""
import numpy as np
import scipy as sp
from scipy.optimize import leastsq
__all__ = ['fit_photon_histograms']

## TO SPEED UP A LITTLE BIT, CHANGE TO FUNCTIONS IMPORTED FROM CMODULES (THIS REQUIRES CYTHON INSTALLED AND: python setup.py build_ext --inplace)
#from Cmodules import gaussian_model, errfunc
gaussian_model = lambda x,p: p[0] * np.exp( - (np.power(x-p[1],2)) / (2*np.power(p[2],2)) ) # the gaussian model
errfunc = lambda p,x,y : gaussian_model(x,p) - y # the error metric

class FitError(Exception):
    def __init__(self, value): self.parameter = value
    def __str__(self): return repr(self.parameter)

class HistError(Exception):
    def __init__(self, value): self.parameter = value
    def __str__(self): return repr(self.parameter)

def fit_photon_histograms(input):

    # fetch data
    i = input['i']
    hist = input['data']['histogram']
    is_bad_pixel = input['data']['mask']
    bins, = input['args']

    # default values
    bgA, bgmu, bgsigma = 3 * [np.nan]
    A, mu, sigma = 3 * [np.nan]

    status = 'ok' # Data is ok
    if ~is_bad_pixel:
        try:
            ##################################
            # 1. Fitting the background peak #
            ##################################

            # some preliminary estimates
            index_bg_peak = 1 + np.argmax(hist[1:-1])

            # initial parameters
            bgA_start = hist[index_bg_peak]
            bgmu_start = bins[index_bg_peak] + 0.001
            bgsigma_start = 3.5
            p0 = np.hstack([bgA_start, bgmu_start, bgsigma_start])

            # mask out the photon peaks
            bghmask = np.ones(bins.shape[0]).astype(np.bool)
            bghmask[index_bg_peak + int(4*bgsigma_start):] = False
            bghmask[0] = False

            # use leastsq without bounds
            p, ier = leastsq(errfunc, p0, args=(bins[bghmask], hist[bghmask]), maxfev=20000)
            if ier > 4: raise FitError("scipy.optimize.leastsq did not converge")
            bgA, bgmu, bgsigma = p

#            if True: raise HistError("The first bin in the histogram is higher than the fitted bgA parameter")
            if hist[0] >= bgA : raise HistError("The first bin in the histogram is higher than the fitted bgA parameter")
            else:

                # estimate the threshold (between background and one photon)
                index_threshold = np.where(bins == np.round(bgmu + 4. * np.abs(bgsigma)))[0][0]
                delta = 5
                valley = hist[index_threshold - delta:index_threshold + delta]
                threshold = np.where(valley == valley.min())[0][-1] + index_threshold - delta

                ##################################
                # 2. Fitting the one photon peak #
                ##################################

                # Estimate the location of the one photon peak
                peak = threshold + np.argmax(hist[threshold:-1])
                P = bins[peak]
                #mu_bound_min = bgmu + P - delta
                #mu_bound_max = bgmu + P + delta

                # initial parameters
                A_start = hist[peak]
                mu_start = P
                sigma_start = bgsigma
                p0 = np.hstack([A_start, mu_start, sigma_start])

                # mask out the background
                cut_min =  int((peak + threshold + 2) / 2.)
                #print 'vorher'
                #print index_threshold, cut_min
                #cut_min = int(threshold)
                #print 'nachher'
                #print cut_min
                phmask = np.ones(bins.shape[0]).astype(np.bool)
                phmask[:cut_min] = False
                phmask[-1] = False

                # use leastsq without bounds
                p, ier = leastsq(errfunc, p0, args=(bins[phmask], hist[phmask]), maxfev=20000)
                if ier > 4: raise FitError("scipy.optimize.leastsq did not converge")
                A, mu, sigma = p

        except FitError:
            status = 'fit_error'
        except HistError:
            status = 'hist_error'
        except:
            status = 'other_error'

    else: status = 'is_bad'
    return {'i':i, 'bg_offset':bgmu, 'bg_amp':bgA,  'bg_sigma':bgsigma, 'photon_offset':mu, 'photon_amp':A, 'photon_sigma':sigma, 'status':status}

# Testing the code and its performance
if __name__ == '__main__':
    hist, bins = np.load('testing/test_histogram.npy')
    out = fit_photon_histograms({'i':0, 'data':{'histogram':hist, 'mask':False}, 'args':(bins,)})
    print out
