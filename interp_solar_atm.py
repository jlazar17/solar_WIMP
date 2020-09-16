import numpy as np
from scipy.interpolate import splev, splrep
import matplotlib.pyplot as plt

class Interpolator():

    def __init__(self, csv_path):
        self.paper_e     = np.genfromtxt(csv_path, delimiter=',')[:,0]
        self.paper_flux  = np.genfromtxt(csv_path, delimiter=',')[:,1]
        self._mak_e_pad_arrs()
        self._tck         = splrep(np.log(self._e_pad), np.log(self._flux_pad))
        
    def _mak_e_pad_arrs(self):
        el_pad = [self.paper_e[0]/2**n for n in np.linspace(5,1,5)]
        er_pad = [self.paper_e[-1]*2**n for n in np.linspace(1,5,5)]

        self._e_pad    = np.hstack([el_pad, self.paper_e,er_pad])
        self._flux_pad = np.hstack([np.full(5, self.paper_flux[0]), self.paper_flux, np.full(5, self.paper_flux[-1])])

    def interp_flux(self,e):
        le   = np.log(e)
        val  = splev(le, self._tck)
        flux = np.exp(val)*np.power(e, -3)
        return flux

    def make_check_plot(self):
        ee = np.logspace(1,4,41)
        vals = np.exp(splev(np.log(ee), self._tck))
        plt.plot(self.paper_e, self.paper_flux, label='True')
        plt.plot(ee, vals, ls='--', label='spline')
        plt.loglog()
        plt.savefig('solar_atm_flux_interp_test.png')
        plt.close()

if __name__=='__main__':
    i = Interpolator('/data/user/jlazar/solar_WIMP/data/solar_atm_spectrum.csv')
    i.make_check_plot()
