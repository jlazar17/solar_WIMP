import numpy as np
import solar_position_calc as sc
import sys
import mc

#args = sys.argv
#
#nRun       = int(args[1])
#bgParticle = args[2]
#nuType     = args[3]
#binning    = args[4]
#h5File     = "/data/user/jlazar/solar_WIMP/data/mcRecarray.npy"

import argparse
import gen_rescale_az_zen as gaz

##### SET UP COMMAND LINE ARGUMENTS #####
parser = parser = argparse.ArgumentParser()
parser.add_argument("--bgparticle",
                    type=str,
                    help="type of background particle ('kaon', 'kaon_bar', 'pion', 'pion_bar')"
                   )
parser.add_argument("-f",
                    type=float,
                    default=1.,
                    help="factor to rescale MEOWS"
                   )
parser.add_argument("-n",
                    type=int,
                    help="run number"
                   )
parser.add_argument("--nt",
                    type=str,
                    help="neutrino type (nu or nuBar)"
                   )
parser.add_argument("--binning",
                    type=str,
                    help="bin width (f=0.5 degree, f=1.6 degrees)"
                   )

args           = parser.parse_args()
rescale_factor = args.f
nuType         = args.nt
nRun           = args.n
binning        = args.binning
bgParticle     = args.bgparticle

SKIP = 10
SEED = 111255

np.random.seed(SEED)

dataPath = "/data/user/jlazar/solar_WIMP/data"
#dataPath = "/Users/jlazar/Documents/IceCube/data"
auToCm = 1.496e13  # conversion from au to cm
rSun = 6.9e10  # radius of sun in cm
rProductionRegion = 1e9  # radius of signal production region in cm
solarZenPdf = np.loadtxt("%s/solar_zenith_pdf.txt" % dataPath)
z = np.radians(np.linspace(0, 179, 1800))

deltaT = 60. * 60. * 4.  # seconds in one day
jdStart = 2455349.5
jdEnd = 2457906.5
nDays = jdEnd - jdStart
step = deltaT / (3600. * 24.)
jds = np.linspace(jdStart, jdEnd, int(nDays / step) + 1)
azimuths = np.random.rand(len(jds)) * 2 * np.pi
zeniths = np.random.choice(z, len(jds), p=solarZenPdf) + np.radians(np.random.rand(len(jds)))

binsPerDecade = 10
highExp       = 6.5
lowExp        = 0.5
eBins         = np.logspace(lowExp, highExp, int(binsPerDecade * (highExp-lowExp) + 1))
if binning=="c":
    gammaBins     = np.linspace(0, np.pi, 113)# 113 picked to best match binning on MEOWS wiki page
elif binning=="f":
    gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree


def loadFlux(bgParticle):
    return np.load("%s/%s_flux_marjon.npy" % (dataPath, bgParticle))


def loadMC(h5File, nuType):
    return mc.MonteCarlo(h5File, nuType)


def truncateMC(monteCarlo, nRun):
    monteCarlo.setNuZen(monteCarlo.nuZen[0 + nRun::SKIP])
    monteCarlo.setNuAz(monteCarlo.nuAz[0 + nRun::SKIP])
    monteCarlo.setNuE(monteCarlo.nuE[0 + nRun::SKIP])
    monteCarlo.setRecoZen(monteCarlo.recoZen[0 + nRun::SKIP])
    monteCarlo.setRecoAz(monteCarlo.recoAz[0 + nRun::SKIP])
    monteCarlo.setRecoE(monteCarlo.recoE[0 + nRun::SKIP])
    monteCarlo.setOneWeight(monteCarlo.oneWeight[0 + nRun::SKIP])


def gammaCalc(flux, monteCarlo):
    numGammaTheta = np.zeros((len(gammaBins)-1, len(eBins)-1))
    for i in range(len(jds)):
        monteCarlo.setRecoGamma(zeniths[i], azimuths[i])

        n = flux * monteCarlo.oneWeight * deltaT

        hist = np.histogram2d(monteCarlo.recoGamma, monteCarlo.recoE, bins=[gammaBins, eBins], weights=n)
        numGammaTheta += hist[0]

    return numGammaTheta


def main():
    flux = loadFlux(bgParticle)[0 + nRun::SKIP]
    monteCarlo = loadMC(h5File, nuType)
    truncateMC(monteCarlo, nRun)
    numGammaTheta = gammaCalc(flux, monteCarlo)
    np.save("%s/e_d_theta_hist/partial_hists/%s_%s_%d_energy_delta_theta_hist_bg.npy" % (dataPath, bgParticle, binning, nRun), numGammaTheta)


main()
