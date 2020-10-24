import numpy as np
import os
import argparse

import propa
import charon
from physicsconstants import PhysicsConstants
from controls_2 import theta_23, theta_12, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes


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
    parser.add_argument('--where',
                       default='1AU'
                      )
    args = parser.parse_args()
    return args

qr_ch_dict = {5:"bb", 8:"WW", 11:"tautau"}
ws_cn_dict = {"bb":5, "WW":8, "tautau":11}

class FluxCalculator():
    
    def __init__(self, ch, m, where, gen):
        self.ch    = ch
        self.m     = m
        self.where = where
        self.gen   = gen
        self.evolved_flux = None

    def calc_flux(self):
        if self.gen=='BRW':
            flux = propa.NuFlux(qr_ch_dict[self.ch], self.m, nodes, Emin=e_min, Emax=self.m, bins=nodes,
                                 process='ann', theta_12=theta_12, theta_13=theta_13, 
                                 theta_23=theta_23, delta=delta, delta_m_12=delta_m_12, 
                                 delta_m_13=delta_m_13, interactions=True)
        elif self.gen=='pythia':
            flux = propa.NuFlux(qr_ch_dict[self.ch], self.m, nodes, Emin=e_min, Emax=self.m, bins=nodes,
                                 process='ann', theta_12=theta_12, theta_13=theta_13, 
                                 theta_23=theta_23, delta=delta, delta_m_12=delta_m_12, 
                                 delta_m_13=delta_m_13, interactions=True, pathFlux='/data/user/qliu/DM/DMFlux/Pythia/no_EW/Sun/results/%s_%d_Sun.dat' % (ch_dict[ch], m))
        else:
            print('wrong gen')
            quit()
        evolved_flux = flux.Sun(self.where)
        self.evolved_flux = evolved_flux
        return evolved_flux
    
    def save_flux(self):
        if self.evolved_flux is None:
            print('Flux not set yet. Setting flux..')
            self.calc_flux()
        np.save("/data/user/jlazar/solar_WIMP/data/charon_fluxes/ch%d-m%d_%s_%s_dn_dz.npy" % (self.ch, self.m, self.where, self.gen), self.evolved_flux)
            

if __name__=="__main__":
    args = initialize_args()
    fc = FluxCalculator(args.ch, args.m, args.where, args.whichgen)
    fc.calc_flux()
    fc.save_flux()
