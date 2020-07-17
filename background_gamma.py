import numpy as np
import solar_position_calc as sc
from controls import *
from path_gen import PathGen
from os import path
import weight_MC as wmc
import gen_rescale_az_zen as gaz
import h5py

from base_gamma import BaseGamma

class BackgroundGamma(BaseGamma):
      
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

                reco_gamma = gaz.opening_angle(self.reco_zen, self.reco_az, zenith, azimuths[i])
                
                h = np.histogram2d(reco_gamma, self.reco_e, bins=[gamma_bins, e_bins], weights=self.flux*self.oneweight)
                hist += h[0]

            self.gamma_hist = hist*self._skip

    def save_gamma_hist(self):
        np.save(self.mcfg.get_e_d_theta_path(self.fluxtype), self.gamma_hist)

