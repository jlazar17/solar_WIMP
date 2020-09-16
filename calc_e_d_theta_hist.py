import time
t0 = time.time()
print(t0)
print('importing numpy')
import numpy as np
print('took %f seconds' % (time.time()-t0))
t = time.time()
print('importing solar_position_calc')
import solar_position_calc as sc
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing controls')
from controls import *
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing PathGen')
from path_gen import PathGen
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing path')
from os import path
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing argparse')
import argparse
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing weighter')
import weight_MC as wmc
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing MCReader')
from mc_reader import MCReader
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing h5py')
import h5py
print('took %f seconds' % (time.time()-t))
t = time.time()
print('importing rescale thing')
import gen_rescale_az_zen as gaz
print('took %f seconds' % (time.time()-t))

def initialize_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mcfile',
                        type=str,
                        help='path to mcfile to be used'
                       )
    parser.add_argument('--fluxtype',
                        type=str,
                        help='fluxtype to be used, e.g. conv-numu or ch5-m1000'
                       )
    parser.add_argument('--options',
                        type=str,
                        default='00'
                       )
    
    args = parser.parse_args()
    return args

################################################################################
class BackgroundGamma():
      
    def __init__(self, mcpath, fluxtype, options, _skip=1):
        print('bg gamma')
        self._skip  = _skip
        self.slc    = slice(None, None, _skip)
        self.mcpath = mcpath
        self.h5f    = h5py.File(mcpath, "r")
        self.mcfg   = PathGen(self.mcpath)
        self.fluxtype = fluxtype
        self.flux = np.load(self.mcfg.get_mc_dn_dz_path(self.fluxtype))[self.slc]
        self._seed  = abs(hash(self.mcfg.get_mc_dn_dz_path(self.flux))) % (2**32)
        self.mcr    = MCReader(self.mcpath, self.slc, options=options)


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

    def save_gamma_hist(self, options):
        np.save('data/test_this', self.gamma_hist)
        print('%s_%s' % (self.mcfg.get_e_d_theta_path(self.fluxtype), options))
        #np.save(self.mcfg.get_e_d_theta_path(self.fluxtype), self.gamma_hist)

################################################################################
class SignalGamma():
      
    def __init__(self, mcpath, fluxtype, options, _skip=1):
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
        self.mcr    = MCReader(self.mcpath, self.slc, options=options)


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
                             self.flux[m] * \
                             self.mcr.oneweight[m],
                             0
                            )
            else:
                n = np.where(nu_gamma <= gamma_cut,
                             self.flux[m] *             \
                             self.mcr.oneweight[m] *    \
                             (1. / solar_solid_angle) * \
                             (1. / (4*np.pi*np.power(rad, 2))),
                             0
                            )
            h = np.histogram2d(reco_gamma, self.mcr.reco_e[m], bins=[gamma_bins, e_bins], weights=n)
            hist += h[0]
            
 
        self.gamma_hist = hist*self._skip

    def save_gamma_hist(self, options):
        print('%s_%s' % (self.mcfg.get_e_d_theta_path(self.fluxtype), options))
        #np.save('%s_%s' % (self.mcfg.get_e_d_theta_path(self.fluxtype), options), self.gamma_hist)
        np.save('data/yes_please_test_me', self.gamma_hist)

def main(mcpath, fluxtype, options, _skip=1):
    if fluxtype=='conv-numu':
        print('bg')
        gamma = BackgroundGamma(mcpath, fluxtype, options, _skip=_skip)
    else:
        print('signal')
        gamma = SignalGamma(mcpath, fluxtype, options, _skip=_skip)
    print('setting hist')
    gamma.set_gamma_hist()
    print("saving hist")
    gamma.save_gamma_hist(options)

if __name__=='__main__':
    args=initialize_args()
    main(args.mcfile, args.fluxtype, args.options, _skip=1000)
