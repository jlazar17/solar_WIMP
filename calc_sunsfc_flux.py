import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import argparse
import os

import config_copy
from physicsconstants import PhysicsConstants

param = PhysicsConstants()

qr_ch_dict = {5:"bb", 8:"WW", 11:"tautau"}
ws_chan_dict = {"bb":5, "WW":8, "tautau":11}

e_min     = 10 # GeV
nodes     = 200
xsec_path = '/data/user/qliu/DM/GOLEMTools/sources/nuSQuIDS/data/xsections/'
data_path = "/Users/jlazar/Documents/IceCube/data/"

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
                       )
    parser.add_argument("--delta",
                        type=float,
                        default=2.523e-3
                       )
    args = parser.parse_args()
    return args


if __name__=='__main__':
    args       = initialize_parser()
    e_min      = 10 # GeV
    ch         = args.ch
    m          = args.m
    e_max     = m # GeV
    theta_12   = args.theta_12
    theta_13   = args.theta_13
    theta_23   = args.theta_23
    delta_m_12 = args.delta_m_12
    delta_m_13 = args.delta_m_13
    delta      = args.delta
    nodes      = 200
    fname = '%d_%d_%f_%f_%f_%f_%f_%f.npy' % (ch, m, theta_12, theta_13, theta_23, delta_m_12, delta_m_13, delta)
    print(fname)
    #dn_dz = config_copy.NuFlux_Solar("Pythia", e_min, e_max, nodes, qr_ch_dict[ch], m, param,
    #                            theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
    #                            delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
    #                            interactions=True, xsec=xsec_path, location="detector")
    dn_dz = config_copy.NuFlux_Solar("Pythia", e_min, e_max, nodes, qr_ch_dict[ch], m, param,
                                theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                                delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                                interactions=True, location="detector")

    np.save('/data/user/jlazar/solar_WIMP/data/param_uncertainty/%s' % fname, dn_dz)

