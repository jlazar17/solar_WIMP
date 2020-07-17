import numpy as np
import solar_position_calc as sc
from controls import *
from path_gen import PathGen
from os import path
import weight_MC as wmc
from mc_reader import MCReader
import gen_rescale_az_zen as gaz
import h5py

class BackgroundGamma():
      
    def __init__(self, mcpath, _skip=1):
        print('bg gamma')
        self._skip  = _skip
        self.slc    = slice(None, None, _skip)
        self.mcpath = mcpath
        self.h5f    = h5py.File(mcpath, "r")
        self.mcfg   = PathGen(self.mcpath)
        self.fluxtype = 'conv-numu'
        self.flux   = np.load(self.mcfg.get_mc_dn_dz_path(self.fluxtype))[self.slc]
        self._seed  = abs(hash(self.mcfg.get_mc_dn_dz_path(self.flux))) % (2**32)
        self.mcr    = MCReader(self.mcpath, self.slc)

#    def set_mc_quantities(self):
#        self.oneweight = self.get_oneweight()
#        self.nu_e      = self.h5f["NuEnergy"]["value"][self.slc]
#        self.nu_zen    = self.h5f["NuZenith"]["value"][self.slc]
#        self.nu_az     = self.h5f["NuAzimuth"]["value"][self.slc]
#        self.reco_e    = self.h5f["MuExEnergy"]["value"][self.slc]
#        self.reco_zen  = self.h5f["MuExZenith"]["value"][self.slc]
#        self.reco_az   = self.h5f["MuExAzimuth"]["value"][self.slc]
#
#        print(self.nu_e[:100])
#        print(self.oneweight[:100])
#        
#    def get_oneweight(self):
#        if path.exists(self.mcfg.get_ow_path()):
#            ow = np.load(self.mcfg.get_ow_path())[self.slc]*1e-4
#        else:
#            try:
#                ow = self.h5f["oneweight"]["value"][self.slc]*1e-4
#            except:
#                ow = wmc.weight_mc(self.mcfg)[self.slc]*1e-4
#        return ow

    def set_gamma_hist(self):
        hist = np.zeros((len(gamma_bins)-1, len(e_bins)-1))
        for i, jd in enumerate(jds):
            x = sc.nParameter(jd)
            obl = sc.solarObliquity(x)
            L = sc.L(x)
            G = sc.g(x)
            lamb = sc.solarLambda(L,G)
            rad = sc.solarR(G)
            solar_solid_angle = 2*np.pi*(1-np.cos(r_sun/rad))
            zenith = sc.equatorialZenith(obl, lamb)

            reco_gamma = gaz.opening_angle(self.mcr.reco_zen, self.mcr.reco_az, zenith, azimuths[i])
            
            h = np.histogram2d(reco_gamma, self.mcr.reco_e, bins=[gamma_bins, e_bins], weights=self.flux*self.mcr.oneweight)
            hist += h[0]

        self.gamma_hist = hist*self._skip

    def save_gamma_hist(self):
        np.save(self.mcfg.get_e_d_theta_path(self.fluxtype), self.gamma_hist)

class SignalGamma():
      
    def __init__(self, mcpath, fluxtype, _skip=1):
        print('aignal gamma')
        self._skip  = _skip
        self.slc    = slice(None, None, _skip)
        self.mcpath = mcpath
        self.fluxtype = fluxtype
        self.h5f    = h5py.File(mcpath, "r")
        self.mcfg   = PathGen(self.mcpath)
        self.flux   = np.load(self.mcfg.get_mc_dn_dz_path(self.fluxtype))[self.slc]
        self._seed  = abs(hash(self.mcfg.get_mc_dn_dz_path(self.flux))) % (2**32)
        self.mcr    = MCReader(self.mcpath, self.slc)

#    def set_mc_quantities(self):
#        self.nu_e      = self.h5f["NuEnergy"]["value"][self.slc]
#        self.nu_zen    = self.h5f["NuZenith"]["value"][self.slc]
#        self.nu_az     = self.h5f["NuAzimuth"]["value"][self.slc]
#        self.reco_e    = self.h5f["MuExEnergy"]["value"][self.slc]
#        self.reco_zen  = self.h5f["MuExZenith"]["value"][self.slc]
#        self.reco_az   = self.h5f["MuExAzimuth"]["value"][self.slc]
#        self.oneweight = self.get_oneweight()
#
#        
#    def get_oneweight(self):
#        if path.exists(self.mcfg.get_ow_path()):
#            ow = np.load(self.mcfg.get_ow_path())[self.slc]*1e-4
#        else:
#            try:
#                ow = self.h5f["oneweight"]["value"][self.slc]*1e-4
#            except:
#                ow = wmc.weight_mc(self.mcfg, save=True)[self.slc]*1e-4
#        return ow

    def set_gamma_hist(self):
       hist = np.zeros((len(gamma_bins)-1, len(e_bins)-1))
       for i, jd in enumerate(jds):
           x = sc.nParameter(jd)
           obl = sc.solarObliquity(x)
           L = sc.L(x)
           G = sc.g(x)
           lamb = sc.solarLambda(L,G)
           rad = sc.solarR(G)
           solar_solid_angle = 2*np.pi*(1-np.cos(r_sun/rad))
           zenith = sc.equatorialZenith(obl, lamb)
 
           gamma_cut = np.arctan(r_sun / rad)
           zmax     = zenith+gamma_cut
           zmin     = zenith-gamma_cut
           amax     = azimuths[i]+gamma_cut
           amin     = azimuths[i]-gamma_cut
           m1 = np.logical_and(self.mcr.nu_zen>zmin, self.mcr.nu_zen<zmax)
           m2 = np.logical_and((self.mcr.nu_az>amin%(2*np.pi)), self.mcr.nu_az<amax%(2*np.pi))
           m  = np.logical_and(m1, m2)
 
           nu_gamma   = gaz.opening_angle(self.mcr.nu_zen[m], self.mcr.nu_az[m], zenith, azimuths[i])
           reco_gamma = gaz.opening_angle(self.mcr.reco_zen[m], self.mcr.reco_az[m], zenith, azimuths[i])
           
           n = np.where(nu_gamma <= gamma_cut,
                        self.flux[m] *                          \
                        self.mcr.oneweight[m] *                             \
                        (1. / solar_solid_angle) *       \
                        (1. / (4*np.pi*np.power(rad, 2))),
                        0
                       )
           h = np.histogram2d(reco_gamma, self.mcr.reco_e[m], bins=[gamma_bins, e_bins], weights=n)
           hist += h[0]
 
       self.gamma_hist = hist*self._skip

    def save_gamma_hist(self):
        np.save(self.mcfg.get_e_d_theta_path(self.fluxtype), self.gamma_hist)

def main(mcpath, fluxtype, _skip=1):
    if fluxtype=='conv-numu':
        gamma = BackgroundGamma(mcpath, _skip=_skip)
    else:
        gamma = SignalGamma(mcpath, fluxtype, _skip=_skip)
    print('setting hist')
    gamma.set_gamma_hist()
    print("saving hist")
    gamma.save_gamma_hist()

if __name__=='__main__':
    from sys import argv as args
    main(args[1], args[2], _skip=10)
