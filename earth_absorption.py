import numpy as np
from scipy.interpolate import interp2d

e = np.genfromtxt("/data/user/jlazar/solar_WIMP/data/fluxes_flavor.txt")[:,2]
n         = len(set(e))
e         = e[:n]
cz        = np.genfromtxt("/data/user/jlazar/solar_WIMP/data/fluxes_flavor.txt")[:,1][::n]
nu_abs    = np.genfromtxt("/data/user/jlazar/solar_WIMP/data/fluxes_flavor.txt")[:,4].reshape(len(cz), len(e))
nubar_abs = np.genfromtxt("/data/user/jlazar/solar_WIMP/data/fluxes_flavor.txt")[:,7].reshape(len(cz), len(e))

nu_interp    = interp2d(np.log10(e), np.asarray(cz), nu_abs)
nubar_interp = interp2d(np.log10(e), np.asarray(cz), nubar_abs)

def calc_absorption(event):
    """
    Args:
        event (tuple): tuple containing neutrino energy (GeV), nuetrino cos(zenith), and neutrino type (14 or -14)
    Returns:
        attenuation factor for neutrino type of input energy and cos(zenith)
    """
    e, cz, i = event
    if i==14:
        interp = nu_interp
    elif i==-14:
        interp = nubar_interp
    else:
        print("particle type %d is not valid. Must be 14 or -14" % i)
    return interp(np.log10(e), cz)
