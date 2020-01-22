import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import config
from sys import argv as args
from physicsconstants import PhysicsConstants

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

data_path = "/data/user/jlazar/solar_WIMP/data/"

flux = config.NuFlux_Solar("Pythia", e_min, e_max, nodes, ch, m, param,
                           theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                           delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                           interactions=True, xsec=xsec_path, angle=90
                          )

print(np.asarray([tup[2] for tup in flux]) * float(m))
