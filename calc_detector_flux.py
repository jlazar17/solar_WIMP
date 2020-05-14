import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import config_copy
from physicsconstants import PhysicsConstants
import os
import argparse

param = PhysicsConstants()

##### SET UP COMMAND LINE ARGUMENTS #####
parser = argparse.ArgumentParser()
parser.add_argument("-f",
                    type=float,
                    default=1.,
                    help="factor to rescale MEOWS"
                   )
parser.add_argument("--ch",
                    type=int,
                    help="WIMPSim channel number bb:5, WW:8, tautau:11"
                   )
parser.add_argument("-m",
                    type=int,
                    help="Dark matter mass"
                   )
args = parser.parse_args()
n_zen = 101
zens  = np.linspace(80, 180, n_zen)

qr_ch_dict = {5:"bb", 8:"WW", 11:"tautau"}
ws_cn_dict = {"bb":5, "WW":8, "tautau":11}

theta_12   = 33.82 # degrees
theta_23   = 48.3 # degrees
theta_13   = 8.61 # degrees
delta_m_12 = 7.39e-5 # eV^2
delta_m_13 = 2.523e-3 # eV^2
delta      = 222 # degrees

e_min     = 10 # GeV
e_max     = args.m # GeV
nodes     = 200
xsec_path = '/data/user/qliu/DM/GOLEMTools/sources/nuSQuIDS/data/xsections/'

def calc_flux(ch, m, params):
    dn_dz = np.zeros((2, n_zen, nodes))
    for i, zen in enumerate(zens):
        f = config_copy.NuFlux_Solar("Pythia", e_min, e_max, nodes, qr_ch_dict[ch], m, param,
                                theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                                delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                                interactions=True, xsec=xsec_path, angle=zen)
        nu_mu_dn_dz     = np.asarray([tup[2] for tup in f]) * float(m)
        nu_mu_bar_dn_dz = np.asarray([tup[5] for tup in f]) * float(m)
        dn_dz[0][i][:] = nu_mu_dn_dz
        dn_dz[1][i][:] = nu_mu_bar_dn_dz
    return dn_dz

#def calc_flux(ch, m, params):
#    qr_ch = qr_ch_dict[ch]
#    dn_dz = np.zeros((2, n_zen * nodes))
#    for i, zen in enumerate(zens):
#        f = config_copy.NuFlux_Solar("Pythia", e_min, e_max, nodes, qr_ch, m, params,
#                                theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
#                                delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
#                                interactions=True, xsec=xsec_path, angle=zen, logscale=True)
#        nu_mu_dn_dz     = np.asarray([tup[2] for tup in f]) * float(m)
#        nu_mu_bar_dn_dz = np.asarray([tup[5] for tup in f]) * float(m)
#        dn_dz[0][i*nodes:(i+1)*nodes] = nu_mu_dn_dz
#        dn_dz[1][i*nodes:(i+1)*nodes] = nu_mu_bar_dn_dz
#    return dn_dz.T

if __name__=="__main__":
    dn_dz = calc_flux(args.ch, args.m, param)
    np.save("/data/user/jlazar/solar_WIMP/data/qr_dn_dz/ch%d-m%d_dn_dz.npy" % (args.ch, args.m), dn_dz)
