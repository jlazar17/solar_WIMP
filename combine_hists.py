import numpy as np
from sys import argv

ch = int(argv[1])
m  = int(argv[2])

ns  = np.linspace(0,99,100)

shape = (112, 60)

def any_nan(arr, desc_str):
    if np.any(np.isnan(arr)):
        print("nan found in %s" % desc_str)
        return True
    else:
        return False

def set_data_path():
    import re
    import os
    global data_path
    r = re.compile('cobalt.*.icecube.wisc.edu')
    if os.popen('hostname').readline().rstrip("\n")=='dyn-8-50.icecube.wisc.edu':
        data_path = "/Users/jlazar/Documents/IceCube/data/e_d_theta_hist/"
    else:
        data_path = "/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/"


set_data_path()
bad_ones = []

nu_flux     = np.zeros(shape)
nu_bar_flux = np.zeros(shape)
print("ch==%d" % ch)
print("m==%d" % m)
for n in ns:
    nu_flux     += np.load("%s/partial_hists/ch%d_m%d_nu_%d_energy_delta_theta_hist.npy" %(data_path, ch, m, n))
    #if any_nan(nu_flux, "nu_flux"):
    if np.any(np.isnan(nu_flux)):
        bad_ones.append((ch, m, n, "nu"))
    nu_bar_flux += np.load("%s/partial_hists/ch%d_m%d_nuBar_%d_energy_delta_theta_hist.npy" %(data_path, ch, m, int(n)))
    if any_nan(nu_bar_flux, "nu_bar_flux"):
        bad_ones.append((ch, m, n, "nu_bar"))
    
nu_flux     = nu_flux.T[::-1]
nu_bar_flux = nu_bar_flux.T[::-1]

nu_flux[np.where(nu_flux<1e-50)]         = 0 # round small fluxes to zero
nu_bar_flux[np.where(nu_bar_flux<1e-50)] = 0

tot_flux    = nu_flux + nu_bar_flux

np.save("%s/ch%d_m%d_nu_e_d_theta_hist.npy" % (data_path, ch, m), nu_flux)
np.save("%s/ch%d_m%d_nu_bar_e_d_theta_hist.npy" % (data_path, ch, m), nu_bar_flux)
np.save("%s/ch%d_m%d_tot_e_d_theta_hist.npy" % (data_path, ch, m), tot_flux)

