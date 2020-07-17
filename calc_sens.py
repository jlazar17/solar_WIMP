import numpy as np
from path_gen import PathGen
from controls import si_ann_rate_dict, sd_ann_rate_dict, delta_t
from scipy.special import gamma
from scipy.optimize import broyden1

# convert cross sections to pb for computational ease
si_anns = {key:val*1.0e-3 for key,val in si_ann_rate_dict.items()}
sd_anns = {key:val*1.0e-3 for key,val in sd_ann_rate_dict.items()}

def factorial(x):
    x = np.asarray(x)
    return gamma(x+1)

class SensCalc():

    def __init__(self, sigflux):
        self.sigflux = sigflux
        self.mass   = int(self.sigflux.split('m')[-1])
        self.mub    = None
        self.mus    = None

    def add_mub(self, mcpath):
        pg = PathGen(mcpath)
        if self.mub is None:
            self.mub = delta_t*np.load(pg.get_e_d_theta_path("conv-numu"))
        else:
            self.mub += delta_t*np.load(pg.get_e_d_theta_path("conv-numu"))

    def add_mus(self, mcpath):
        pg = PathGen(mcpath)
        factor = 1./self.mass * sd_anns[self.mass] * delta_t
        if self.mus is None:
            print(pg.get_e_d_theta_path(self.sigflux))
            self.mus = factor*np.load(pg.get_e_d_theta_path(self.sigflux))
        else:
            print(pg.get_e_d_theta_path(self.sigflux))
            self.mus += factor*np.load(pg.get_e_d_theta_path(self.sigflux))

    def neg_llh(self, ns, nb):
        """
        ns  float : observed number of signal events
        nb  float : observed number of bacground events
        """
        nb  = np.asarray(nb)
        ns  = np.asarray(ns)
        mu = self.mus*ns + self.mub*nb
        lh = -np.log(np.exp(-mu) * np.power(mu, self.mub) / factorial(self.mub))
        lh = lh[np.where(~np.isnan(lh))]
        lh = lh[np.where(~np.isinf(lh))]
        return np.sum(lh)

    def set_sens(self, chi2=2.71, init_guess=0.5):
        nb = 1
        nllh0 = self.neg_llh(0, nb)
        delta_llh = lambda n:2*(self.neg_llh(n, nb) - nllh0)-chi2
        try:
            sens = broyden1(delta_llh, init_guess) * 1.e-39 # cm^{-2}
        except:
            sens = broyden1(delta_llh, 50) * 1.e-39 # cm^{-2}
        self.sens = sens
