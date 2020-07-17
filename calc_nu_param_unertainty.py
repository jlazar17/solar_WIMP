import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import os
import argparse

import config
from physicsconstants import PhysicsConstants
from uncertainty_controls import qr_ch_dict, ws_ch_dict
param = PhysicsConstants()

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
    parser.add_argument("--ordering",
                        type=str,
                        help="Mass ordering. Normal (no) or inverted (io)"
                       )
    parser.add_argument("-n",
                        type=int,
                        help="Number of times to run",
                        default=100
                       )
    args = parser.parse_args()
    return args

def make_param_arrays(ordering, n=100):
    if ordering=='no':
        th12s  = np.random.normal(33.8, 0.77, n)
        th23s  = np.random.normal(48.3, 1.5, n)
        th13s  = np.random.normal(8.61, 1.3, n)
        m12s   = np.random.normal(7.39e-5, 0.205e-5, n)
        m13s   = np.random.normal(2.523e-3, 0.031e-3, n)
        deltas = np.random.normal(222, 33, n)
    elif ordering=='io':
        th12s  = np.random.normal(33.8, 0.77, n)
        th23s  = np.random.normal(48.3, 1.3, n)
        th13s  = np.random.normal(8.61, 0.125, n)
        m12s   = np.random.normal(7.39e-5, 0.205e-5, n)
        m13s   = np.random.normal(-2.509e-3, 0.031e-3, n)
        deltas = np.random.normal(285, 25, n)
    else:
        print('invalid mass ordering. must be either "no" or "io"')
        quit()
    return np.vstack([th12s, th23s, th13s, m12s, m13s, deltas]).T

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

def main(ch, m, ordering, n):
    neu_params = make_param_arrays(ordering, n=n)
    savedir = '/data/user/jlazar/solar_WIMP/data/param_uncertainties/'
    for (th12, th23, th13, m12, m13, delta) in neu_params:
        savefname = '%d_%d_%s_%f_%f_%f_%f_%f_%f.npy' % (ch, m, ordering, th12, th23, th13, m12, m13, delta)
        print(savefname)
        flux = calc_flux(ch, m, th12, th23, th13, m12, m13, delta, 10, 200, param)
        np.save("%s/%s" % (savedir, savefname), flux)

if __name__=='__main__':
    args = initialize_parser()
    ch       = args.ch
    m        = args.m
    ordering = args.ordering
    n        = args.n
    seed = hash(str(ch)+str(m)+ordering) % 2**32
    np.random.seed(seed)
    main(ch, m, ordering, n)
