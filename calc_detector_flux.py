import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import config_copy
from sys import argv as args
from physicsconstants import PhysicsConstants
import os

param = PhysicsConstants()

ch = args[1]
m  = int(args[2])

n_zen = 101
zens  = np.linspace(80, 180, n_zen)

ws_chan_dict = {"bb":5, "WW":8, "tautau":11}

theta_12   = 33.82 # degrees
theta_23   = 48.3 # degrees
theta_13   = 8.61 # degrees
delta_m_12 = 7.39e-5 # eV^2
delta_m_13 = 2.523e-3 # eV^2
delta      = 222 # degrees

e_min     = 10 # GeV
e_max     = m # GeV
nodes     = 200
xsec_path = '/data/user/qliu/DM/GOLEMTools/sources/nuSQuIDS/data/xsections/'

def set_data_path():
    import re
    import os
    global data_path
    r = re.compile('cobalt.*.icecube.wisc.edu')
    if os.popen('hostname').readline().rstrip("\n")=='dyn-8-50.icecube.wisc.edu':
        data_path = "/Users/jlazar/Documents/IceCube/data/"
    else:
        data_path = "/data/user/jlazar/solar_WIMP/data/"

set_data_path()


ch_num = ws_chan_dict[ch]
dn_dz = np.zeros((2, n_zen, nodes))
for i, zen in enumerate(zens):
    f = config_copy.NuFlux_Solar("Pythia", e_min, e_max, nodes, ch, m, param,
                            theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                            delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                            interactions=True, xsec=xsec_path, angle=zen)
    nu_mu_dn_dz     = np.asarray([tup[2] for tup in f]) * float(m)
    nu_mu_bar_dn_dz = np.asarray([tup[5] for tup in f]) * float(m)
    dn_dz[0][i][:] = nu_mu_dn_dz
    dn_dz[1][i][:] = nu_mu_bar_dn_dz
#if ~os.path.isdir("%s/qr_dn_dz" % data_path):
#    os.mkdir("%s/qr_dn_dz" % data_path)
np.save("%s/qr_dn_dz/ch%d_m%d_dn_dz_new.npy" % (data_path, ch_num, m), dn_dz)
