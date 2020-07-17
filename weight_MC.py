import numpy as np
import argparse
import LeptonWeighter as LW
import tables
from path_gen import PathGen

default_xs_location = "/data/user/jlazar/solar_WIMP/data/xs_splines/"
default_flux_location = "/data/user/jlazar/solar_WIMP/data/fluxes/"
dsdxdy_nu_CC    = default_xs_location+"dsdxdy_nu_CC_iso.fits"
dsdxdy_nubar_CC = default_xs_location+"dsdxdy_nubar_CC_iso.fits"
dsdxdy_nu_NC    = default_xs_location+"dsdxdy_nu_NC_iso.fits"
dsdxdy_nubar_NC = default_xs_location+"dsdxdy_nubar_NC_iso.fits"
nuSQFluxKaon    = default_flux_location+"kaon_atmospheric.hdf5"
nuSQFluxPion    = default_flux_location+"kaon_atmospheric.hdf5"
nuSQFluxConv    = default_flux_location+"atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5"
def initialize_argparse():
    parser = argparse.ArgumentParser(description='Run LeptonWeighter python example.')
    parser.add_argument('--mcfile')
    parser.add_argument('--dsdxdy_nu_CC', default=default_xs_location+"dsdxdy_nu_CC_iso.fits")
    parser.add_argument('--dsdxdy_nubar_CC', default=default_xs_location+"dsdxdy_nubar_CC_iso.fits")
    parser.add_argument('--dsdxdy_nu_NC', default=default_xs_location+"dsdxdy_nu_NC_iso.fits")
    parser.add_argument('--dsdxdy_nubar_NC', default=default_xs_location+"dsdxdy_nubar_NC_iso.fits")
    parser.add_argument('--nuSQFluxKaon', default=default_flux_location+"kaon_atmospheric.hdf5")
    parser.add_argument('--nuSQFluxPion', default=default_flux_location+"kaon_atmospheric.hdf5")
    parser.add_argument("--nuSQFluxConv", default=default_flux_location+"atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5")
    return parser.parse_args()

# Load h5 file in chunks for memory limited situations
def chonk_h5file(mcgf, chonk_size=1e5):
    h5file  = tables.File(mcgf.mcpath)
    n_items = h5file.root.FinalStateX.shape[0]
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
#    h5file.close()
    return chonk_slices


# Set up LeptonWeighter
def initialize_weighter(mcgf):
    simulation_generators = LW.MakeGeneratorsFromLICFile(mcgf.get_licfile_path())
    
    #kaon_nusquids_flux = LW.nuSQUIDSAtmFlux(args.nuSQFluxKaon)
    #pion_nusquids_flux = LW.nuSQUIDSAtmFlux(args.nuSQFluxPion)
    conv_nusquids_flux = LW.nuSQUIDSAtmFlux(nuSQFluxConv)
    xs = LW.CrossSectionFromSpline(dsdxdy_nu_CC,dsdxdy_nubar_CC,
                                   dsdxdy_nu_NC,dsdxdy_nubar_NC)
    
    weighter = LW.Weighter([conv_nusquids_flux],xs,simulation_generators)
    return weighter

def weight_mc(mcgf,save=False):
    weighter = initialize_weighter(mcgf)
    mc  = tables.File(mcgf.mcpath)
    weights = np.zeros(mc.root.FinalStateX.shape[0])
    chunks  = chonk_h5file(mcgf)
    for chunk in chunks:
        chunk_weights = np.zeros(int(chunk.stop-chunk.start))
        events = zip(
                     mc.root.FinalStateX[chunk]['value'],
                     mc.root.FinalStateY[chunk]['value'],
                     mc.root.FinalType0[chunk]['value'],
                     mc.root.FinalType1[chunk]['value'],
                     mc.root.NuAzimuth[chunk]['value'],
                     mc.root.NuZenith[chunk]['value'],
                     mc.root.NuEnergy[chunk]['value'],
                     mc.root.PrimaryType[chunk]['value'],
                     mc.root.TotalColumnDepth[chunk]['value'],
                    )
         
        for i, event in enumerate(events):
            LWevent                        = LW.Event()
            LWevent.interaction_x          = event[0]
            LWevent.interaction_y          = event[1]
            LWevent.final_state_particle_0 = LW.ParticleType(event[2])
            LWevent.final_state_particle_1 = LW.ParticleType(event[3])
            LWevent.azimuth                = event[4]
            LWevent.zenith                 = event[5]
            LWevent.energy                 = event[6]
            LWevent.primary_type           = LW.ParticleType(event[7])
            LWevent.total_column_depth     = event[8]
            LWevent.x                      = 0.
            LWevent.y                      = 0.
            LWevent.z                      = 0.
    
            chunk_weights[i] = weighter.get_oneweight(LWevent)/2.
    
        weights[chunk.start:chunk.stop] = chunk_weights
    
    
    if save:
        np.save("/data/user/jlazar/solar_WIMP/data/mc_oneweight/%s" % mcgf.get_mcname(), weights)
#    mc.close()
    return weights

if __name__=="__main__":
    args = initialize_argparse()
    savefilename = args.mcfile.split("/")[-1][:-3]+ ".npy"
    mcgf = PathGen(args.mcfile) 
    weights = weight_mc(mcgf)
    np.save("/data/user/jlazar/solar_WIMP/data/mc_oneweight/%s" % mcgf.get_mcname(), weights)
