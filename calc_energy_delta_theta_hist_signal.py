import tables
import numpy as np
import solar_position_calc as sc
from sys import argv as args
import mc
import os
import argparse
import h5py
import gen_rescale_az_zen as gaz
import weight_MC as wmc
from file_gen import File_Gen

SKIP=1

def initialize_parser():
    ##### SET UP COMMAND LINE ARGUMENTS #####
    default_xs_location = "/data/user/jlazar/solar_WIMP/data/xs_splines/"
    default_flux_location = "/data/user/jlazar/solar_WIMP/data/fluxes/"
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",
                        type=float,
                        default=1.,
                        help="factor to rescale MEOWS"
                       )
    parser.add_argument("--ch",
                        type=int,
                        help="WIMPSim channel number bb:5, WW:8, tautau:11"
                       )
    parser.add_argument("-m",
                        type=int,
                        help="Dark matter mass"
                       )
    parser.add_argument("--binning",
                        type=str,
                        help="bin width (f=0.5 degree, f=1.6 degrees)"
                       )
    parser.add_argument("--mcfile",
                        type=str,
                        help="path to MC file"
                       )
    parser.add_argument("--flux",
                        type=str,
                        help="path to flux file"
                       )
    parser.add_argument("--fluxfile",
                        type=str,
                        help="path to flux file"
                       )
    parser.add_argument("--LICFile",
                        type=str,
                        help="path to LIC file used in generation of MC"
                       )
    parser.add_argument('--dsdxdy_nu_CC', default=default_xs_location+"dsdxdy_nu_CC_iso.fits")
    parser.add_argument('--dsdxdy_nubar_CC', default=default_xs_location+"dsdxdy_nubar_CC_iso.fits")
    parser.add_argument('--dsdxdy_nu_NC', default=default_xs_location+"dsdxdy_nu_NC_iso.fits")
    parser.add_argument('--dsdxdy_nubar_NC', default=default_xs_location+"dsdxdy_nubar_NC_iso.fits")
    parser.add_argument('--nuSQFluxKaon', default=default_flux_location+"kaon_atmospheric.hdf5")
    parser.add_argument('--nuSQFluxPion', default=default_flux_location+"kaon_atmospheric.hdf5")
    parser.add_argument("--nuSQFluxConv", default="/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes/Flux_AIRS_sib_HG_th24_dm2/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5")
    parser.add_argument('--OutFile', default="out")
    return parser.parse_args()

rSun        = 6.9e10  # radius of sun in cm
#solarZenPdf = np.loadtxt("/data/user/jlazar/solar_WIMP/data/solar_zenith_pdf.txt")
#z           = np.radians(np.linspace(0, 179, 1800))

deltaT = 60. * 30. # sec
jdStart = 2455349.5
jdEnd = 2457906.5
nDays = jdEnd - jdStart
nStep = deltaT / (3600. * 24.)
jds = np.linspace(jdStart, jdEnd, int(nDays / nStep) + 1)
azimuths = np.random.rand(len(jds)) * 2 * np.pi # solar azimuth positions
#zeniths = np.random.choice(z, len(jds), p=solarZenPdf) # solar zenith positions

gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree

binsPerDecade = 10
highExp       = 6.5
lowExp        = 0.5
eBins         = np.logspace(lowExp, highExp, int(binsPerDecade * (highExp-lowExp) + 1))

def gammaCalc(oneweight, dn_dz, mcfile,  SKIP):
    mc     = tables.File(mcfile, "r")
    nu_e   = mc.root.NuEnergy[:]["value"][::SKIP]
    nu_zen = mc.root.NuZenith[:]["value"][::SKIP]
    nu_az  = mc.root.NuAzimuth[:]["value"][::SKIP]
    reco_e = mc.root.MuExEnergy[:]["value"][::SKIP]
    reco_zen = mc.root.MuExZenith[:]["value"][::SKIP]
    reco_az  = mc.root.MuExAzimuth[:]["value"][::SKIP]
   # else:
   #     reco_az, reco_zen = gaz.gen_new_zen_az(rescale_factor)
   #     reco_az = reco_az[::SKIP]
   #     reco_zen = reco_zen[::SKIP]
    gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree
    gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree
    gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree
    gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree
    gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree
    numGammaTheta    = np.zeros((len(gammaBins)-1, len(eBins)-1))
    for i, jd in enumerate(jds):
        x = sc.nParameter(jd)
        obl = sc.solarObliquity(x)
        L = sc.L(x)
        G = sc.g(x)
        lamb = sc.solarLambda(L,G)
        rad = sc.solarR(G)
        solar_solid_angle = 2*np.pi*(1-np.cos(rSun/rad))
        zenith = sc.equatorialZenith(obl, lamb)

        gammaCut = np.arctan(rSun / rad)
        zmax     = zenith+gammaCut
        zmin     = zenith-gammaCut
        amax     = azimuths[i]+gammaCut
        amin     = azimuths[i]-gammaCut

        m1 = np.logical_and(nu_zen>zmin, nu_zen<zmax)
        m2 = np.logical_and((nu_az>amin%(2*np.pi)), nu_az<amax%(2*np.pi))
        m  = np.logical_and(m1, m2)

        nu_gamma   = gaz.opening_angle(nu_zen[m], nu_az[m], zenith, azimuths[i])
        reco_gamma = gaz.opening_angle(reco_zen[m], reco_az[m], zenith, azimuths[i])
        n = np.where(nu_gamma <= gammaCut,
                     dn_dz[m] *                          \
                     oneweight[m] *                             \
                     (1. / solar_solid_angle) *       \
                     (1. / (4*np.pi*np.power(rad, 2))),
                     0
                    )
        
        
        hist = np.histogram2d(reco_gamma, reco_e[m], bins=[gammaBins, eBins], weights=n)
        numGammaTheta += hist[0]

    return numGammaTheta*SKIP

#def get_mcname(mcfile):
#    if mcfile.split("/")[7]=="Systematics":
#        mcname = mcfile.split("/")[8]
#    else:
#        mcname = "Nominal"
#    return mcname

def main(args, oneweight):
    mcgf =  File_Gen(args.mcfile) 
    dn_dz = np.load(mcgf.get_mc_dn_dz_path(args.flux))[::SKIP]
    numGammaTheta = gammaCalc(oneweight, dn_dz, args.mcfile, SKIP)
    fluxname = mcgf.get_mc_dn_dz_path(args.flux).split("/")[-1].split("_")[0]
    mcname   = mcgf.get_mcname()
    np.save("/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/%s_%s_e_d_theta.npy" % (fluxname, mcname), numGammaTheta.T[::-1])


if __name__=="__main__":

    
    args           = initialize_parser()
    #rescale_factor = args.f
    #ch             = args.ch
    #m              = args.m
    #binning        = args.binning
    #mc_file        = args.mcfile
    #flux_file      = args.fluxfile
    mcgf = File_Gen(args.mcfile) 
    SEED = abs(hash(mcgf.get_mc_dn_dz_path(args.flux))) % (2**32)
    
    np.random.seed(SEED)
    
    mc     = h5py.File(args.mcfile, "r")
    try:
        ow = mc["oneweight"]["value"][::SKIP]*1e-4
    except:
        ow = wmc.weight_mc(mcgf)[::SKIP]*1.e-4
    
    mc     = tables.File(args.mcfile, "r")
    
    nu_e   = mc.root.NuEnergy[:]["value"][::SKIP]
    nu_zen = mc.root.NuZenith[:]["value"][::SKIP]
    nu_az  = mc.root.NuAzimuth[:]["value"][::SKIP]
    reco_e = mc.root.MuExEnergy[:]["value"][::SKIP]
    reco_zen = mc.root.MuExZenith[:]["value"][::SKIP]
    reco_az  = mc.root.MuExAzimuth[:]["value"][::SKIP]
    #if rescale_factor==1:
    #    print("Using MC reco angles")
    #    reco_zen = mc.root.MuExZenith[:]["value"][::SKIP]
    #    reco_az  = mc.root.MuExAzimuth[:]["value"][::SKIP]
    #else:
    #    print("Generating rescaled reco angles")
    #    reco_az, reco_zen = gaz.gen_new_zen_az(rescale_factor)
    #    reco_zen = reco_zen[::SKIP]
    #    reco_az  = reco_az[::SKIP]

    main(args,ow)
