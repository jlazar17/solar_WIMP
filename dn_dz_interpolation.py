import numpy as np
import h5py
# from scipy.stats import gaussian_kde as kde
from scipy.interpolate import griddata
import os
import argparse

##### SET UP COMMAND LINE ARGUMENTS #####
parser = parser = argparse.ArgumentParser()
parser.add_argument("--ch",
                    type=int,
                    help="WIMPSim channel number bb:5, WW:8, tautau:11"
                   )
parser.add_argument("-m",
                    type=int,
                    help="Dark matter mass"
                   )
parser.add_argument("--mcfile",
                    type=str,
                    help="Path to monte carlo h5 file to be used"
                   )
parser.add_argument("--fluxfile",
                    type=str,
                    help="Path to flux file to be interpolated"
                   )

def get_desc_name(path):
    return path.split("/")[-1].split("_")[0]

def interp_flux(mcfile, fluxfile):

    nzen    = 100
    ne      = 350
    points  = [(i,j) for i in np.linspace(-1,0.2,nzen)
              for j in np.logspace(2,6,ne)]

    vals = np.log10(np.load(fluxfile))
    vals = np.where(np.isinf(vals), -500, vals)

    h5f = h5py.File(mcfile, "r")

    mcflux = np.power(10, np.where(h5f["PrimaryType"]["value"]==14,
                                   # If condition is true use numu flux values
                                   griddata(points, vals[0], list(zip(np.cos(h5f["NuZenith"]["value"]),h5f["NuEnergy"]["value"])), method="linear"),
                                   # If condition is false use numubar flux values
                                   griddata(points, vals[1], list(zip(np.cos(h5f["NuZenith"]["value"]),h5f["NuEnergy"]["value"])), method="linear"),
                                  )
                     )
    return mcflux

def main(args):
    mcflux    = interp_flux(args.mcfile, args.fluxfile)
    desc_name = get_desc_name(args.fluxfile)
    print(desc_name)
    np.save("/data/user/jlazar/solar_WIMP/data/mc_flux/%s_mc_flux.npy" % desc_name, mcflux)

if __name__=="__main__":
    main(parser.parse_args())
