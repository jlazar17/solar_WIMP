import numpy as np
import argparse
import h5py
import LeptonWeighter as LW
import tables

# parser = argparse.ArgumentParser(description='Run LeptonWeighter python example.')

# default_xs_location = "/net/cvmfs_users/carguelles/sterile/resources/CrossSections/"
#default_flux_location = "/net/cvmfs_users/carguelles/sterile/resources/Fluxes/conventional/"
##default_flux_location = "/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/LeptonWeighter/resources/data/"
#parser.add_argument("--fluxfile")
#parser.add_argument('--LICFile')
#parser.add_argument('--mcfile')
#parser.add_argument('--dsdxdy_nu_CC', default=default_xs_location+"dsdxdy_nu_CC_iso.fits")
#parser.add_argument('--dsdxdy_nubar_CC', default=default_xs_location+"dsdxdy_nubar_CC_iso.fits")
#parser.add_argument('--dsdxdy_nu_NC', default=default_xs_location+"dsdxdy_nu_NC_iso.fits")
#parser.add_argument('--dsdxdy_nubar_NC', default=default_xs_location+"dsdxdy_nubar_NC_iso.fits")
##parser.add_argument('--nuSQFluxKaon', default=default_flux_location+"HondaGaisserKaon.hdf5")
##parser.add_argument('--nuSQFluxPion', default=default_flux_location+"HondaGaisserPion.hdf5")
#parser.add_argument('--nuSQFluxKaon', default=default_flux_location+"kaon_atmospheric.hdf5")
#parser.add_argument('--nuSQFluxPion', default=default_flux_location+"kaon_atmospheric.hdf5")
#parser.add_argument("--nuSQFluxConv", default="/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes/Flux_AIRS_sib_HG_th24_dm2/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5")
#parser.add_argument('--OutFile', default="out")

# parser.set_defaults(use_cache=True)
# args=parser.parse_args()

# Load h5 file in chunks for memory limited situations
def chonk_h5file(h5file, chonk_size=1e5):

    n_items = h5file.root.FinalStateX.shape[0]
    n_chonks = int(np.ceil(n_items/chonk_size))

    for chonk_number in range(n_chonks):
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

		# in if statement?
		yield slice(n0, n1)

# Set up LeptonWeighter
# put these in a dictionary? 
def initialize_weighter(lic_file, nuSQFFluxConv, dsdxdy_nu_CC, dsdxdy_nubar_CC, dsdxdy_nu_NC, dsdxdy_nubar_NC):
    simulation_generators = LW.MakeGeneratorsFromLICFile(lic_file)

    conv_nusquids_flux = LW.nuSQUIDSAtmFlux(nuSQFluxConv)
    
    xs = LW.CrossSectionFromSpline(dsdxdy_nu_CC,dsdxdy_nubar_CC,
                                   dsdxdy_nu_NC,dsdxdy_nubar_NC)
    
    weighter = LW.Weighter([conv_nusquids_flux],xs,simulation_generators)

    return weighter

def weight_mc(mc_filename, lic_file, nuSQFFluxConv, dsdxdy_nu_CC, dsdxdy_nubar_CC, dsdxdy_nu_NC, dsdxdy_nubar_NC): 
    weighter = initialize_weighter(lic_file, nuSQFFluxConv, dsdxdy_nu_CC, dsdxdy_nubar_CC, dsdxdy_nu_NC, dsdxdy_nubar_NC)

    mc  = tables.File(mc_filename)

    weights = np.zeros(mc.root.FinalStateX.shape[0])
    chunks  = chonk_h5file(mc_filename)

    for chunk in chonk_h5file(args.mcfile):
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

	mc.close()
    return weights

if __name__=="__main__":
	mc_filename = ""
	lic_filename = "" 
	nuSQFFluxConv = ""
	dsdxdy_nu_CC = ""
	dsdxdy_nubar_CC = ""
	dsdxdy_nu_NC = ""
	dsdxdy_nubar_NC = ""

    savefilename = args.mcfile.split("/")[-1][:-3]+ ".npy"
    print(savefilename)

    weights = weight_mc(mc_filename, lic_file, nuSQFFluxConv, dsdxdy_nu_CC, dsdxdy_nubar_CC, dsdxdy_nu_NC, dsdxdy_nubar_NC)

    np.save("/data/user/jlazar/solar_WIMP/data/mc_oneweight/%s" % savefilename, weights)
