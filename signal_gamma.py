import numpy as np
import solar_position_calc as sc
from controls import *
from path_gen import PathGen
from os import path
import weight_MC as wmc
import gen_rescale_az_zen as gaz
import h5py

from base_gamma import BaseGamma

class SignalGamma(BaseGamma):
      
    def set_gamma_hist(self):
        if path.exists(self.mcfg.get_e_d_theta_path(self.fluxtype)):
            self.gamma_hist = np.load(self.mcfg.get_e_d_theta_path(self.fluxtype))
        else:
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
                m1 = np.logical_and(self.nu_zen>zmin, self.nu_zen<zmax)
                m2 = np.logical_and((self.nu_az>amin%(2*np.pi)), self.nu_az<amax%(2*np.pi))
                m  = np.logical_and(m1, m2)
    
                nu_gamma   = gaz.opening_angle(self.nu_zen[m], self.nu_az[m], zenith, azimuths[i])
                reco_gamma = gaz.opening_angle(self.reco_zen[m], self.reco_az[m], zenith, azimuths[i])
                
                n = np.where(nu_gamma <= gamma_cut,
                             self.flux[m] *                          \
                             self.oneweight[m] *                             \
                             (1. / solar_solid_angle) *       \
                             (1. / (4*np.pi*np.power(rad, 2))),
                             0
                            )
                h = np.histogram2d(reco_gamma, self.reco_e[m], bins=[gamma_bins, e_bins], weights=n)
                hist += h[0]
    
            self.gamma_hist = hist*self._skip

    def save_gamma_hist(self):
        np.save(self.mcfg.get_e_d_theta_path(self.fluxtype), self.gamma_hist)

