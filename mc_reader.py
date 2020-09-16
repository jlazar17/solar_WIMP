import numpy as np
import h5py
from os import path

from path_gen import PathGen
import weight_MC as wmc
from gen_rescale_az_zen import gen_new_zen_az
class MCReader():

    def __init__(self, mcpath, slc, options='00'):
        self.mcpath = mcpath
        self._rescale  = int(options[0])
        self._scramble = int(options[1])
        self.num_files = float(mcpath.split('_')[-1][:-3])/5.
        self.slc    = slc
        self.mcfg   = PathGen(self.mcpath)
        self.h5f    = h5py.File(mcpath, "r")
        self.set_mc_quantities()

    def set_mc_quantities(self):
        self.nu_e      = self.h5f["NuEnergy"]["value"][self.slc]
        if self._scramble:
            print('Scramble mode: activated')
            delta_az = np.random.rand(len(self.nu_e))*2*np.pi
        else:
            delta_az = np.zeros(len(self.nu_e))
        self.nu_zen    = self.h5f["NuZenith"]["value"][self.slc]
        self.nu_az     = self.h5f["NuAzimuth"]["value"][self.slc]
        self.reco_e    = self.h5f["MuExEnergy"]["value"][self.slc]
        if not self._rescale:
            self.reco_zen  = self.h5f["MuExZenith"]["value"][self.slc]
            _ = self.h5f["MuExAzimuth"]["value"][self.slc]
            self.reco_az   = np.mod(_+delta_az, 2*np.pi)
        else:
            print('Rescale mode: activated')
            reco_az, reco_zen = gen_new_zen_az(self.h5f["NuEnergy"]["value"])
            self.reco_az = np.mod(reco_az[self.slc]+delta_az,2*np.pi)
            self.reco_zen = reco_zen[self.slc]
        self.set_oneweight()
        
    def set_oneweight(self):
        if path.exists(self.mcfg.get_ow_path()):
            ow = np.load(self.mcfg.get_ow_path())[self.slc]/self.num_files
        else:
            ow = wmc.weight_mc(self.mcfg)[self.slc]/self.num_files
        self.oneweight = ow

    def set_compare():
        self.oneweight = np.where(self.nu_e>3000, 0, self.oneweight)
        

if __name__=='__main__':
    from  sys import argv as args
    mcpath = args[1]
    print(args[2])
    slc = slice(None, None, int(args[2]))
    mcr = MCReader(mcpath, slc)
