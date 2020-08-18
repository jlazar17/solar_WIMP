import tables
import numpy as np
# from scipy.stats import gaussian_kde as kde
from scipy.interpolate import griddata
import os
import argparse

##### SET UP COMMAND LINE ARGUMENTS #####
parser = parser = argparse.ArgumentParser()
parser.add_argument("--fluxfile",
                    type=str,
                    help="path to flux file to be interpolated"
                   )
parser.add_argument("--mcfile",
                    type=str,
                    help="path to Monte Carlo file to be used"
                   )
parser.add_argument("-m",
                    type=int
                   )

zens  = np.linspace(80, 180, 101)
e_min = 10

def chonk_h5file(h5file, chonk_size=1e5):
    h5f = tables.File(h5file)
    n_items = h5f.root.FinalStateX.shape[0]
    chonks = int(np.ceil(n_items/chonk_size))
    chonk_slices = []
    for chonk_number in range(chonks):
        a = int(np.floor(float(n_items) / float(chonks)))
        b = int(np.ceil(float(n_items) / float(chonks)))
        x = n_items - a*chonks
        if chonk_number < x:
            n = b
            n0 = n*chonk_number
        else:
            n = a
            n0 = b*x + a*(chonk_number - x)
        n1 = min(n0 + n, n_items)
        if n0 >= n_items:
            continue
        chonk_slices.append(slice(n0,n1))
    h5f.close()
    return chonk_slices

def get_mcname(mcfile):
    if mcfile.split("/")[7]=="Systematics":
        mcname = mcfile.split("/")[8]
    else:
        mcname = "Nominal"
    return mcname

def interp_dn_dz(args, savepath):
    mc = tables.File(args.mcfile)
    chunks  = chonk_h5file(args.mcfile)
    mc_flux = np.zeros(mc.root.PrimaryType.shape[0])
    dN_dz  = np.load(args.fluxfile)
    nu_flux     = dN_dz[0]
    nu_bar_flux = dN_dz[1]
    n_zen       = nu_flux.shape[0]
    if args.m is not None:
        n_en        = nu_flux.shape[1]
        es          = np.tile(np.linspace(e_min, args.m, n_en), n_zen)
    else:
        e = np.genfromtxt('data/solar_atm/PostPropagation/SIBYLL2.3_pp_CombinedGHAndHG_H4a_nu.txt')[:,0]
        n_en = len(e)
        print(e)
        es = np.tile(e, n_zen)
    zs          = np.concatenate([np.full(n_en, zen) for zen in zens])
    # Prepare points and values arrays
    points        = np.vstack([es,zs]).T
    nu_values     = np.log10(np.concatenate(nu_flux))
    nu_bar_values = np.log10(np.concatenate(nu_bar_flux))
    for chunk in chunks:
        # Load data files

        # Get all mc quantities ready
        nu_i           = np.where(mc.root.PrimaryType[chunk]["value"]==14)[0]
        nu_bar_i       = np.where(mc.root.PrimaryType[chunk]["value"]==-14)[0]
        nu_e           = mc.root.NuEnergy[chunk]["value"][nu_i]
        nu_bar_e       = mc.root.NuEnergy[chunk]["value"][nu_bar_i]
        nu_zen         = np.degrees(mc.root.NuZenith[chunk]["value"][nu_i])
        nu_bar_zen     = np.degrees(mc.root.NuZenith[chunk]["value"][nu_bar_i])



        # Replace -infinities from log(0) with large neg values
        nu_values[np.where(np.isinf(nu_values))]         = -500
        nu_bar_values[np.where(np.isinf(nu_bar_values))] = -500

        # Interpolate
        nu_gd          = griddata(points, nu_values,     (nu_e, nu_zen),         method="linear")
        nu_bar_gd      = griddata(points, nu_bar_values, (nu_bar_e, nu_bar_zen), method="linear")
        nu_interp      = np.power(10, nu_gd)
        nu_bar_interp  = np.power(10, nu_bar_gd)
        
        if args.m is not None:
            nu_interp[np.where(nu_e>args.m)[0]]         = 0
            nu_bar_interp[np.where(nu_bar_e>args.m)[0]] = 0

        # Save interpolated fluxes
        mc_flux[chunk][nu_i]     = nu_interp
        mc_flux[chunk][nu_bar_i] = nu_bar_interp
        mc_flux[np.where(np.isnan(mc_flux))] = 0
    return mc_flux

def main(args):
    fluxname = args.fluxfile.split("/")[-1].split("_")[0]
    mcname = get_mcname(args.mcfile)
    if fluxname=='SIBYLL2.3':
        savepath = "/data/user/jlazar/solar_WIMP/data/solar_atm/AtIceCube/interped/%s_%s.npy" % (fluxname, mcname)
    else:
        savepath = "/data/user/jlazar/solar_WIMP/data/mc_dn_dz/%s_%s_dn_dz.npy" % (fluxname, mcname)
    mc_flux = interp_dn_dz(args, savepath)
    np.save(savepath,  mc_flux)

main(parser.parse_args())
