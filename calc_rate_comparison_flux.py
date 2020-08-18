import numpy as np
import os
import argparse

import propa
import charon
from physicsconstants import PhysicsConstants

param = PhysicsConstants()

def initialize_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ch",
                        type=int,
                        help="WIMPSim channel number bb:5, WW:8, tautau:11"
                       )
    parser.add_argument("-m",
                        type=int,
                        help="Dark matter mass"
                       )
    args = parser.parse_args()
    return args

ch_dict = {5:'bb', 8:'WW', 11:'tautau'}

theta_12   = 33.82 # degrees
theta_23   = 48.3 # degrees
theta_13   = 8.61 # degrees
delta_m_12 = 7.39e-5 # eV^2
delta_m_13 = 2.523e-3 # eV^2
delta      = 222 # degrees

Emin     = 10 # GeV
nodes     = 200

if __name__=='__main__':
    args = initialize_args()
    ch = args.ch
    m  = args.m
    nfBRW = propa.NuFlux(ch_dict[ch], m, 200, Emin=Emin, Emax=m, process='ann', bins=200,
                         theta_12=theta_12, theta_23=theta_23, theta_13=theta_13,
                         delta_m_12=delta_m_12, delta_m_13=delta_m_13)
    nfPYT =  propa.NuFlux(ch_dict[ch], m, 200, Emin=Emin, Emax=m, process='ann', bins=200,
                          theta_12=theta_12, theta_23=theta_23, theta_13=theta_13,
                          delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                          pathFlux='/data/user/qliu/DM/DMFlux/Pythia/no_EW/Sun/results/%s_%d_Sun.dat' % (ch_dict[ch],m)) 
    pyt_au = nfPYT.Sun('1AU')
    np.save('test.npy',pyt_au)
    brw_au = nfBRW.Sun('1AU')
    flux = np.zeros((4, 200))
    flux[0,:] = brw_au['nu_mu']
    flux[1,:] = brw_au['nu_mu_bar']
    flux[2,:] = pyt_au['nu_mu']*float(m)
    flux[3,:] = pyt_au['nu_mu_bar']*float(m)
    np.save('/data/user/jlazar/solar_WIMP/data/rate_comparison/ch%d_m%d_rate_comp.npy' % (ch, m), flux)
