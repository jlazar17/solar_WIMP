import numpy as np
# from scipy.stats import gaussian_kde as kde
from scipy.interpolate import griddata
import os
from sys import argv as args

ch    = int(args[1])
m     = int(args[2])
zens  = np.linspace(80, 180, 101)
e_min = 10

#import matplotlib.pyplot as plt
#plt.style.use("paper.mplstyle")

def set_data_path():
    import re
    global data_path
    r = re.compile('cobalt.*.icecube.wisc.edu')
    if r.match(os.popen('hostname').readline().rstrip("\n")) is not None:
        data_path = "/data/user/jlazar/solar_WIMP/data/"
    elif os.popen('hostname').readline().rstrip("\n")=='MBP-FVFXC6EKHV2D.local':
        data_path = "/Users/jlazar/Documents/IceCube/data/"
    else:
        print("Machine not recognized")
        quit()


def interp_dn_dz(ch, m):

    # Load data files
    dN_dz  = np.load("%s/qr_dn_dz/ch%d_m%d_dn_dz.npy" % (data_path, ch, m))
    mc     = np.load("%s/mcRecarray.npy" % data_path)

    # Get all quantities from dN_dE file arranged
    nu_flux     = dN_dz[0]
    nu_bar_flux = dN_dz[1]
    n_zen       = nu_flux.shape[0]
    n_en        = nu_flux.shape[1]
    es          = np.tile(np.linspace(e_min, m, n_en), n_zen)
    zs          = np.concatenate([np.full(n_en, zen) for zen in zens])

    # Get all mc quantities ready
    nu_i           = np.where(mc["i"]==14)[0]
    nu_bar_i       = np.where(mc["i"]==-14)[0]
    nu_e           = mc["nuE"][nu_i]
    nu_bar_e       = mc["nuE"][nu_bar_i]
    nu_zen         = np.degrees(mc["nuZen"][nu_i])
    nu_bar_zen     = np.degrees(mc["nuZen"][nu_bar_i])


    # Prepare points and values arrays
    points        = np.vstack([es,zs]).T
    nu_values     = np.log10(np.concatenate(nu_flux))
    nu_bar_values = np.log10(np.concatenate(nu_bar_flux))

    # Replace -infinities from log(0) with large neg values
    nu_values[np.where(np.isinf(nu_values))]         = -500
    nu_bar_values[np.where(np.isinf(nu_bar_values))] = -500

    # Interpolate
    nu_gd          = griddata(points, nu_values,     (nu_e, nu_zen),         method="linear")
    nu_bar_gd      = griddata(points, nu_bar_values, (nu_bar_e, nu_bar_zen), method="linear")
    nu_interp      = np.power(10, nu_gd)
    nu_bar_interp  = np.power(10, nu_bar_gd)

    nu_interp[np.where(nu_e>m)[0]]         = 0
    nu_bar_interp[np.where(nu_bar_e>m)[0]] = 0

    # Save interpolated fluxes
    mc_flux           = np.zeros(len(mc["i"]))
    mc_flux[nu_i]     = nu_interp
    mc_flux[nu_bar_i] = nu_bar_interp
    mc_flux[np.where(np.isnan(mc_flux))] = 0
    np.save("%s/mc_dn_dz/ch%d_m%d_mc_dn_dz.npy" % (data_path, ch, m), mc_flux)
    return mc_flux

def main():
    set_data_path()
    interp_dn_dz(ch, m)

main()
