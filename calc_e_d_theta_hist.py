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
      
    def __init__(self, mcpath, fluxtype, rescale, _skip=1):
        print('bg gamma')
        self._skip  = _skip
        self.slc    = slice(None, None, _skip)
        self.mcpath = mcpath
        self.h5f    = h5py.File(mcpath, "r")
        self.mcfg   = PathGen(self.mcpath)
        self.fluxtype = fluxtype
        self.flux = np.load(self.mcfg.get_mc_dn_dz_path(self.fluxtype))[self.slc]
        self._seed  = abs(hash(self.mcfg.get_mc_dn_dz_path(self.flux))) % (2**32)
        self.mcr    = MCReader(self.mcpath, self.slc, rescale=rescale)


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

    def save_gamma_hist(self, rescale):
        np.save(self.mcfg.get_e_d_theta_path(self.fluxtype, rescale), self.gamma_hist)

class SignalGamma():
      
    def __init__(self, mcpath, fluxtype, rescale, _skip=1):
        print('signal gamma')
        self._skip  = _skip
        self.slc    = slice(None, None, _skip)
        self.mcpath = mcpath
        self.fluxtype = fluxtype
        self.h5f    = h5py.File(mcpath, "r")
        self.mcfg   = PathGen(self.mcpath)
        if self.fluxtype=='solar-atm':
            self.flux = np.load('/data/user/jlazar/solar_WIMP/data/solar_atm/AtIceCube/interped/SIBYLL2.3_Nominal.npy')[self.slc]
        else:
            self.flux   = np.load(self.mcfg.get_mc_dn_dz_path(self.fluxtype))[self.slc]
        self._seed  = abs(hash(self.mcfg.get_mc_dn_dz_path(self.flux))) % (2**32)
        self.mcr    = MCReader(self.mcpath, self.slc, rescale=rescale)


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

            if self.fluxtype=='solar-atm':
                gamma_cut = np.arctan(1.2*r_sun / rad)
            else: 
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
            if self.fluxtype=='solar-atm':
                n = np.where(nu_gamma <= gamma_cut,
                             self.flux[m] *                          \
                             self.mcr.oneweight[m],
                             0
                            )
            else:
                n = np.where(nu_gamma <= gamma_cut,
                             self.flux[m] *                          \
                             self.mcr.oneweight[m] *                             \
                             (1. / solar_solid_angle) *       \
                             (1. / (4*np.pi*np.power(rad, 2))),
                             0
                            )
            h = np.histogram2d(reco_gamma, self.mcr.reco_e[m], bins=[gamma_bins, e_bins], weights=n)
            hist += h[0]
            
            #nu_gamma   = gaz.opening_angle(self.mcr.nu_zen, self.mcr.nu_az, zenith, azimuths[i])
            #reco_gamma = gaz.opening_angle(self.mcr.reco_zen, self.mcr.reco_az, zenith, azimuths[i])
            #if self.fluxtype=='solar-atm':
            #    n = np.where(nu_gamma <= gamma_cut,
            #                 self.flux[m] *                          \
            #                 self.mcr.oneweight,
            #                 0
            #                )
            #else:
            #    n = np.where(nu_gamma <= gamma_cut,
            #                 self.flux *                          \
            #                 self.mcr.oneweight *                             \
            #                 (1. / solar_solid_angle) *       \
            #                 (1. / (4*np.pi*np.power(rad, 2))),
            #                 0
            #                )
            #h = np.histogram2d(reco_gamma, self.mcr.reco_e, bins=[gamma_bins, e_bins], weights=n)
            #hist += h[0]
 
        self.gamma_hist = hist*self._skip

    def save_gamma_hist(self, rescale):
        np.save(self.mcfg.get_e_d_theta_path(self.fluxtype, rescale), self.gamma_hist)

def main(mcpath, fluxtype, rescale, _skip=1):
    if fluxtype=='conv-numu':
        gamma = BackgroundGamma(mcpath, fluxtype, rescale, _skip=_skip)
    else:
        gamma = SignalGamma(mcpath, fluxtype, rescale, _skip=_skip)
    print('setting hist')
    gamma.set_gamma_hist()
    print("saving hist")
    gamma.save_gamma_hist(rescale)

if __name__=='__main__':
    from sys import argv as args
    main(args[1], args[2], bool(args[3]), _skip=1000)
