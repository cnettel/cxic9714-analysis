#!/usr/bin/env python
"""
A module for the CSPAD detector (back and front)

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last change:  August 8, 2016
"""
import h5py
import numpy as np

class CSPAD:
    def __init__(self, snr=1., nx=1000, ny=1000, alpha=0.7, pixelsize=1., dx=0, dy=0, sigmafile=None, geometryfile=None, maskfile=None):

        # Properties
        self.snr = snr
        self.nx = nx
        self.ny = ny
        self.sh = (ny,nx)
        self.alpha = alpha

        # Pixelsize
        self.pixelsize=pixelsize

        # Geometry
        self._x, self._y   = (None, None)
        self._dx, self._dy = (dx,dy)
        if geometryfile is not None:
            self._load_geometryfile(geometryfile)

        # Mask
        self.mask = np.ones(self.sh).astype(np.bool)
        if maskfile is not None:
            self._load_maskfile(maskfile)
        
        # Sigma
        self.sigma = None
        if sigmafile is not None:
            self._load_sigmafile(sigmafile)

    def _load_geometryfile(self, filename):
        with h5py.File(filename, 'r') as f:
            x = np.floor(f['x'][:]/self.pixelsize).astype(np.int)
            y = np.floor(f['y'][:]/self.pixelsize).astype(np.int)
            x -= x.min()
            y -= y.min()
            self._x, self._y = x,y

    def _load_maskfile(self, filename):
        with h5py.File(filename, 'r') as f:
            self.mask = self.assemble(f['data/data'][:]).astype(np.bool)
            
    def _load_sigmafile(self, filename):
        with h5py.File(filename, 'r') as f:
            self.sigma = self.assemble(f['data/data'][:])
        self.sigma[self.sigma == 0] = -1.
        self.sigma[~self.mask]      = -1.
        assert (self.sh == self.sigma.shape), "Shape mismatch!, (%d,%d) != (%d,%d)" %(self.sh[0], self.sh[1], self.sigma.shape[0], self.sigma.shape[1])

    def assemble(self, array):
        if (self._x is None) or (self._y is None):
            return array
        assembled = np.zeros(self.sh).astype(array.dtype)
        assembled[self._y+self._dy,self._x+self._dx] = array
        return assembled
        
    def gainmap(self):
        if self.sigma is None:
            return np.ones(self.sh) * self.snr
        return self.sigma * self.snr

    def photoncounts(self, array):
        assert (self.sh == array.shape), "Shape mismatch!, (%d,%d) != (%d,%d)" %(self.sh[0], self.sh[1], array.shape[0], array.shape[1])
        pcounts = array / self.gainmap()
        pcounts[pcounts < self.alpha] = 0.
        return np.round(pcounts).astype(np.int64)
