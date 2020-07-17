import numpy as np
import solar_position_calc as sc
from controls import *
from path_gen import PathGen
from os import path
import weight_MC as wmc
import gen_rescale_az_zen as gaz
import h5py


class BaseGamma:
      
    def __init__(self, mcpath, fluxtype, _skip=1):
        self._skip  = _skip
        self.slc    = slice(None, None, _skip)
        self.mcpath = mcpath
        self.fluxtype = fluxtype
        self.h5f    = h5py.File(mcpath, "r")
        self.mcfg   = PathGen(self.mcpath)
        self.flux   = np.load(self.mcfg.get_mc_dn_dz_path(self.fluxtype))[self.slc]
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
            ow = np.load(self.mcfg.get_ow_path())[self.slc]*1e-4
        else:
            try:
                ow = self.h5f["oneweight"]["value"][slc]*1e-4
            except:
                ow = wmc.weight_mc(self.mcfg)[self.slc]*1e-4
        return ow

    def set_gamma_hist(self):
        pass

    def save_gamma_hist(self):
        pass
