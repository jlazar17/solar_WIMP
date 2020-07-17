import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import os
import argparse

import config
from physicsconstants import PhysicsConstants

param = PhysicsConstants()
qr_ch_dict = {5:"bb", 8:"WW", 11:"tautau"}
ws_cn_dict = {"bb":5, "WW":8, "tautau":11}
xsec_path = '/data/user/qliu/DM/GOLEMTools/sources/nuSQuIDS/data/xsections/'

def initialize_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ch",
                        type=int,
                        help="WIMPSim channel number bb:5, WW:8, tautau:11"
                       )
    parser.add_argument("-m",
                        type=int,
                        help="Dark matter mass"
                       )
    
    parser.add_argument("--theta_12",
                        type=float,
                        default=33.82
                       )
    parser.add_argument("--theta_23",
                        type=float,
                        default=48.3
                       )
    parser.add_argument("--theta_13",
                        type=float,
                        default=8.61
                       )
    parser.add_argument("--delta_m_12",
                        type=float,
                        default=7.39e-5
                       )
    parser.add_argument("--delta_m_13",
                        type=float,
                        default=2.523e-3
                       )
    parser.add_argument("--delta",
                        type=float,
                        default=222
                       )
    args = parser.parse_args()
    return args

def calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, params):
    e_max = m
    dn_dz = np.zeros((2, nodes))
    #f = config.NuFlux("Pythia", e_min, e_max, nodes, qr_ch_dict[ch], m, param,
    #                       theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
    #                       delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
    #                       interactions=True, xsec=xsec_path, location='Sunsfc')
    f = config.NuFlux("Pythia", e_min, e_max, nodes, qr_ch_dict[ch], m, param,
                           theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                           delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                           interactions=True, location='Sunsfc')
    nu_mu_dn_dz     = np.asarray([tup[2] for tup in f]) * float(m)
    nu_mu_bar_dn_dz = np.asarray([tup[5] for tup in f]) * float(m)
    dn_dz[0][:] = nu_mu_dn_dz
    dn_dz[1][:] = nu_mu_bar_dn_dz
    return dn_dz

if __name__=='__main__':
    args       = initialize_parser()
    e_min      = 10 # GeV
    ch         = args.ch
    m          = args.m
    theta_12   = args.theta_12
    theta_13   = args.theta_13
    theta_23   = args.theta_23
    delta_m_12 = args.delta_m_12
    delta_m_13 = args.delta_m_13
    delta      = args.delta
    nodes      = 200
    if delta_m_13>0:
        ordering = 'no'
    else:
        ordering = 'io'
    fname = '%d_%d_%f_%f_%f_%f_%f_%f.npy' % (ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta)
    print(fname)
    dn_dz      = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param)

    np.save('/data/user/jlazar/solar_WIMP/data/param_uncertainties//%s' % (ch,m,ordering,fname), dn_dz)
