import numpy as np
import h5py
from os import path

from path_gen import PathGen
import weight_MC as wmc
from gen_rescale_az_zen import gen_new_zen_az
class MCReader():

    def __init__(self, mcpath, slc, rescale=False):
        self.mcpath = mcpath
        self.slc    = slc
        self.mcfg   = PathGen(self.mcpath)
        self.h5f    = h5py.File(mcpath, "r")
        self.set_mc_quantities()    

    def set_mc_quantities(self, rescale=False):
        self.nu_e      = self.h5f["NuEnergy"]["value"][self.slc]
        self.nu_zen    = self.h5f["NuZenith"]["value"][self.slc]
        self.nu_az     = self.h5f["NuAzimuth"]["value"][self.slc]
        self.reco_e    = self.h5f["MuExEnergy"]["value"][self.slc]
        if not rescale:
            self.reco_zen  = self.h5f["MuExZenith"]["value"][self.slc]
            self.reco_az   = self.h5f["MuExAzimuth"]["value"][self.slc]
        else:
            self.reco_az, self.reco_zen = gen_rescale_az_zen(self.nu_e)
        self.set_oneweight()
        
    def set_oneweight(self):
        if path.exists(self.mcfg.get_ow_path()):
            ow = np.load(self.mcfg.get_ow_path())[self.slc]*1e-4
        else:
            ow = wmc.weight_mc(self.mcfg)[self.slc]*1e-4
        self.oneweight = ow

    def set_compare():
        self.oneweight = np.where(self.nu_e>3000, 0, self.oneweight)
        

if __name__=='__main__':
    from  sys import argv as args
    mcpath = args[1]
    print(args[2])
    slc = slice(None, None, int(args[2]))
    mcr = MCReader(mcpath, slc)
