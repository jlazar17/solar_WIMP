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
    parser.add_argument("--whichgen",
                        type=str,
                        default='BRW'
                       )
    args = parser.parse_args()
    return args

n_zen = 11
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
nodes     = 200
#xsec_path = '/home/jlazar/programs/golemsolarwimp_py2/sources/nuSQuIDS/data/'
#xsec_path = '/home/jlazar/programs/GOLEM_SOLAR_WIMP_py2/sources/nuSQuIDS/data/xsections/csms.h4'


def calc_flux(ch, m, gen, params):
    dn_dz = np.zeros((2, len(zens), nodes))
    if gen=='BRW':
        print('cool')
        flux = propa.NuFlux(qr_ch_dict[ch], m, nodes, Emin=e_min, Emax=m, bins=nodes,
                             process='ann', theta_12=theta_12, theta_13=theta_13, 
                             theta_23=theta_23, delta=delta, delta_m_12=delta_m_12, 
                             delta_m_13=delta_m_13, interactions=True)
    elif gen=='pythia':
        print('uncool')
        ch_dict = {5:'bb', 8:'WW', 11:'tautau'}
        flux = propa.NuFlux(qr_ch_dict[ch], m, nodes, Emin=e_min, Emax=m, bins=nodes,
                             process='ann', theta_12=theta_12, theta_13=theta_13, 
                             theta_23=theta_23, delta=delta, delta_m_12=delta_m_12, 
                             delta_m_13=delta_m_13, interactions=True, pathFlux='/data/user/qliu/DM/DMFlux/Pythia/no_EW/Sun/results/%s_%d_Sun.dat' % (ch_dict[ch], m))
    else:
        print('wrong gen')
        quit()
    for i, zen in enumerate(zens):
        det_flux = flux.Sun('detector', zenith=zen)
        dn_dz[0][i][:] = det_flux['nu_mu']
        dn_dz[1][i][:] = det_flux['nu_mu_bar']
    return dn_dz

if __name__=="__main__":
    args = initialize_args()
    print('Doing calc...')
    dn_dz = calc_flux(args.ch, args.m, args.whichgen, param)
    if args.whichgen=='BRW':
        print('saving file to /data/user/jlazar/solar_WIMP/data/qr_dn_dz/ch%d-m%d_dn_dz.np' % (args.ch, args.m))
        np.save("/data/user/jlazar/solar_WIMP/data/1AU/ch%d-m%d_dn_dz.npy" % (args.ch, args.m), dn_dz)
    else:
        np.save("/data/user/jlazar/solar_WIMP/data/1AU/ch%d-m%d_dn_dz_pythia.npy" % (args.ch, args.m), dn_dz * float(args.m))
