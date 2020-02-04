import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import config_copy
import argparse
from physicsconstants import PhysicsConstants

ws_ch_dict = {"bb":5,
              "WW":8,
              "tautau":11
             }

##### SET UP COMMAND LINE ARGUMENTS #####
parser = parser = argparse.ArgumentParser()
parser.add_argument("--ch",
                    type=str,
                    help="channel name bb, WW, tautau"
                   )
parser.add_argument("-m",
                    type=int,
                    help="Dark matter mass"
                   )
args = parser.parse_args()
param = PhysicsConstants()

xsec_path = '/data/user/qliu/DM/GOLEMTools/sources/nuSQuIDS/data/xsections/'

theta_12   = 33.82 # degrees
theta_23   = 48.3 # degrees
theta_13   = 8.61 # degrees
delta_m_12 = 7.39e-5 # eV^2
delta_m_13 = 2.523e-3 # eV^2
delta      = 222 # degrees

m     = args.m
ch    = args.ch
nodes = 200
e_min = 10 # GeV
e_max = m

f_pyth = config_copy.NuFlux_Solar("Pythia",e_min, e_max, nodes, ch, m, param, 
                                  theta_12=theta_12, theta_23=theta_23, theta_13=theta_13,
                                  delta_m_12=delta_m_12, delta_m_13=delta_m_13, delta=delta, 
                                  interactions=True,location = 'Sunsfc',xsec=xsec_path
                                 )
#f_ws   = config_copy.NuFlux_Solar("WimpSim",e_min, e_max, nodes, ch, m, param, 
#                                  theta_12=theta_12, theta_23=theta_23, theta_13=theta_13,
#                                  delta_m_12=delta_m_12, delta_m_13=delta_m_13, delta=delta, 
#                                  interactions=True,location = 'Sunsfc',xsec=xsec_path
#                                 )


#flux = np.zeros((2, nodes))
#nu_dn_dz_ws    = np.asarray([tup[2] for tup in f_ws]) * float(m)
#nubar_dn_dz_ws = np.asarray([tup[5] for tup in f_ws]) * float(m)
#
#flux[0,:] = nu_dn_dz_ws
#flux[1,:] = nubar_dn_dz_ws
#
#np.save("data/sunsfc_dn_dz/ch%d_m%d_wimpsim_sunsfc.npy" % (ws_ch_dict[ch], m), flux)

flux = np.zeros((2, nodes))
nu_dn_dz_pyth    = np.asarray([tup[2] for tup in f_pyth]) * float(m)
nubar_dn_dz_pyth = np.asarray([tup[5] for tup in f_pyth]) * float(m)

flux[0,:] = nu_dn_dz_pyth
flux[1,:] = nubar_dn_dz_pyth

np.save("data/sunsfc_dn_dz/ch%d_m%d_pythia_sunsfc.npy" % (ws_ch_dict[ch], m), flux)
