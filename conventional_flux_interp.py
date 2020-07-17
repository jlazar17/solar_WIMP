import numpy as np
from scipy.interpolate import griddata
import earth_absorption
from file_gen import File_Gen
import h5py
airs_file = np.genfromtxt("/data/user/jlazar/solar_WIMP/data/AIRS_flux_sib_HG.dat").T
airs_czs = airs_file[0]
airs_e   = airs_file[1]
airs_numu = airs_file[4]
airs_numubar = airs_file[5]

def get_det_flux(ptype):
    
    tups = [(airs_e[i],airs_czs[i], ptype) for i in range(len(airs_e))]
    absorption = np.asarray([earth_absorption.calc_absorption(tup) for tup in tups]).T[0]
    
    if ptype==14:
        det_flux = airs_numu*absorption
    elif ptype==-14:
        det_flux = airs_numubar*absorption
    else:
        print(f"ptype=={ptype} not valid. Must  be 14 or -14")
    return det_flux

def interp_flux(czen, e, ptype):
    vals = np.log10(get_det_flux(ptype))
    vals = np.where(np.isinf(vals),-500, vals)
    points = np.asarray([airs_e, airs_czs]).T
    interp = np.power(10, griddata(points, vals, (e, czen), method="linear"))
    interp[np.where(np.isnan(interp))] = 0
    return interp

def main(mcfile):
    print(1)
    fg = File_Gen(mcfile)
    mc = h5py.File(mcfile, "r")
    ptype      = mc["PrimaryType"][:]['value']
    mc_flux    = np.zeros(len(ptype))
    numu_i     = np.where(ptype==14)[0]
    numubar_i  = np.where(ptype==-14)[0]
    mc_nu_czen = np.cos(mc['NuZenith'][:]["value"])[numu_i]
    mc_nu_e    = mc['NuEnergy'][:]["value"][numu_i]
    mc_nubar_czen = np.cos(mc['NuZenith'][:]["value"])[numubar_i]
    mc_nubar_e    = mc['NuEnergy'][:]["value"][numubar_i]
    nu_flux    = interp_flux(mc_nu_czen, mc_nu_e, 14)
    nubar_flux    = interp_flux(mc_nubar_czen, mc_nubar_e, -14)
    mc_flux[numu_i] = nu_flux
    mc_flux[numubar_i] = nubar_flux
    save_path =f"/data/user/jlazar/solar_WIMP/data/mc_dn_dz/conv-numu_{fg.get_mcname()}_dn_dz.npy"
    print(save_path)
    np.save(save_path, mc_flux)
if __name__=="__main__":
    import sys
    main(sys.argv[1])
    
