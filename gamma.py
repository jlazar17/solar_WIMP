import numpy as np
import solar_position_calc as sc
from controls import *
from file_gen import File_Gen
from os import path
import weight_MC as wmc
import h5py

class BackgroundGamma():
      
    def __init__(self, mcfile, _skip=1):
        self._skip  = _skip
        self.slc    = slice(None, None, _skip)
        self.mcfile = mcfile
        self.h5f    = h5py.File(mcfile, "r")
        self.mcfg   = File_Gen(self.mcfile)
        self.flux   = np.load(self.mcfg.get_mc_dn_dz_path('conv-numu'))
        self._seed  = abs(hash(self.mcfg.get_mc_dn_dz_path(self.flux))) % (2**32)
        self.set_mc_quantities()

    def set_mc_quantities(self):
        self.oneweight = self.get_oneweight()
        self.nu_e      = self.h5f["NuEnergy"]["value"][self.slc]
        self.nu_zen    = self.h5f["NuZenith"]["value"][self.slc]
        self.nu_az     = self.h5f["NuAzimuth"]["value"][self.slc]
        self.reco_e    = self.h5f["MuExEnergy"]["value"][self.slc]
        self.reco_zen  = self.h5f["MuExZenith"]["value"][self.slc]
        self.reco_az   = self.h5f["MuExAzimuth"]["value"][self.slc]

        
    def get_oneweight(self):
        if path.exists(self.mcfg.get_ow_path()):
            ow = np.load(self.mcfg.get_ow_path())[self.slc]
        else:
            try:
                ow = self.h5f["oneweight"]["value"][slc]*1e-4
            except:
                ow = wmc.weight_mc(self.mcfg)[self.slc]*1e-4
        return ow

    def gamma_calc(self):
        hist = np.zeros(len(gamma_bins)-1, len(e_bins)-1)
        for i, jd in enumerate(jds):
            x = sc.nParameter(jd)
            obl = sc.solarObliquity(x)
            L = sc.L(x)
            G = sc.g(x)
            lamb = sc.solarLambda(L,G)
            rad = sc.solarR(G)
            solar_solid_angle = 2*np.pi*(1-np.cos(rSun/rad))
            zenith = sc.equatorialZenith(obl, lamb)

            gammaCut = np.arctan(rSun / rad)

            reco_gamma = gaz.opening_angle(self.reco_zen, self.reco_az, zenith, azimuths[i])
            
            hist = np.histogram2d(reco_gamma, reco_e, bins=[gammaBins, eBins], weights=self.flux*self.oneweight)
            numGammaTheta += hist[0]

        return numGammaTheta*_skip
