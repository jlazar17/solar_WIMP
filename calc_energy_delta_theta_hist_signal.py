import numpy as np
import solar_position_calc as sc
from sys import argv as args
import mc
import os
import argparse
import gen_rescale_az_zen as gaz

##### SET UP COMMAND LINE ARGUMENTS #####
parser = parser = argparse.ArgumentParser()
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
ch             = args.ch
m              = args.m
nuType         = args.nt
nRun           = args.n
binning        = args.binning

ms  = [200,300,400,500,600,700,800,900,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000]
chans = [5,8,11]

SKIP = 25

np.random.seed(73)
ri  = np.random.randint(100000, size=(len(chans), len(ms)))
SEED_DICT = {chan:{} for chan in chans}
for i, chan in enumerate(chans):
    for j, mChi in enumerate(ms):
        SEED_DICT[chan][mChi] = ri[i][j]


if nuType=="nu":
    SEED = SEED_DICT[ch][m] + nRun
else:
    SEED = SEED_DICT[ch][m] - 1 - nRun

np.random.seed(SEED)

data_path    = "/data/user/jlazar/solar_WIMP/data/"
rSun        = 6.9e10  # radius of sun in cm
solarZenPdf = np.loadtxt("%s/solar_zenith_pdf.txt" % data_path)
z           = np.radians(np.linspace(0, 179, 1800))

deltaT = 60. * 30. # sec
jdStart = 2455349.5
jdEnd = 2457906.5
nDays = jdEnd - jdStart
nStep = deltaT / (3600. * 24.)
jds = np.linspace(jdStart, jdEnd, int(nDays / nStep) + 1)
azimuths = np.random.rand(len(jds)) * 2 * np.pi # solar azimuth positions
zeniths = np.random.choice(z, len(jds), p=solarZenPdf) # solar zenith positions

binsPerDecade = 10
highExp       = 6.5
lowExp        = 0.5
eBins         = np.logspace(lowExp, highExp, int(binsPerDecade * (highExp-lowExp) + 1))
if binning=="c":
    gammaBins = np.linspace(0, np.pi, 113) # 113 picked to best match binning on MEOWS wiki page
elif binning=="f":
    gammaBins = np.linspace(0, np.pi, 361) # bins of width one half degree



##### LOAD MC INFORMATION #####
mc   = np.load("/data/user/jlazar/solar_WIMP/data/mcRecarray.npy")
if nuType=="nu":
    indices = np.where(mc["i"]==14)[0]
else:
    indices = np.where(mc["i"]==-14)[0]
ow     = mc["oneWeight"][indices][nRun::SKIP]*1e-4
nu_e   = mc["nuE"][indices][nRun::SKIP]
nu_zen = mc["nuZen"][indices][nRun::SKIP]
nu_az  = mc["nuAz"][indices][nRun::SKIP]
reco_e = mc["recoE"][indices][nRun::SKIP]
if rescale_factor==1:
    print("Using MC reco angles")
    reco_zen = mc["recoZen"][indices][nRun::SKIP]
    reco_az  = mc["recoAz"][indices][nRun::SKIP]
else:
    print("Generating rescaled reco angles")
    reco_az, reco_zen = gaz.gen_new_zen_az(rescale_factor)
    reco_zen = reco_zen[indices][nRun::SKIP]
    reco_az  = reco_az[indices][nRun::SKIP]


dn_dz = np.load("%s/mc_dn_dz/ch%d_m%d_mc_dn_dz.npy" % (data_path, ch, m))[indices][nRun::SKIP]

#def loadFlux(nuType, ch, m):
#    mc_recarray = np.load("%s/mcRecarray.npy" % data_path)
#    if nuType == "nu":
#        i  = np.where(mc_recarray["i"]==14)[0]
#        #return np.load("%s/ch%d_m%d_nu_mu_flux.npy" % (dataPath, ch, m))
#    elif nuType == "nuBar":
#        #return np.load("%s/ch%d_m%d_nu_mu_bar_flux.npy" % (dataPath, ch, m))
#        i  = np.where(mc_recarray["i"]==-14)[0]
#    return np.load("%s/mc_dn_dz/ch%d_m%d_mc_dn_dz.npy" % (data_path, ch, m))[i]
#
#
#def loadMC(mc_rec_array_file, nuType):
#    return mc.MonteCarlo(mcFile, nuType)
#
#
#def truncateMC(monteCarlo, nRun):
#    monteCarlo.setNuZen(monteCarlo.nuZen[nRun::SKIP])
#    monteCarlo.setNuAz(monteCarlo.nuAz[nRun::SKIP])
#    monteCarlo.setNuE(monteCarlo.nuE[nRun::SKIP])
#    monteCarlo.setRecoZen(monteCarlo.recoZen[nRun::SKIP])
#    monteCarlo.setRecoAz(monteCarlo.recoAz[nRun::SKIP])
#    monteCarlo.setRecoE(monteCarlo.recoE[nRun::SKIP])
#    monteCarlo.setOneWeight(monteCarlo.oneWeight[nRun::SKIP])


def gammaCalc(dn_dz):
    numGammaTheta    = np.zeros((len(gammaBins)-1, len(eBins)-1))
    print("boing")
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
        #monteCarlo.setTrueGamma(zeniths[i], azimuths[i])
        #monteCarlo.setRecoGamma(zeniths[i], azimuths[i])
        n = np.where(nu_gamma <= gammaCut,
                     dn_dz[m] *                          \
                     ow[m] *                             \
                     (1. / solar_solid_angle) *       \
                     (1. / (4*np.pi*np.power(rad, 2))),
                     0
                    )
        
        
        hist = np.histogram2d(reco_gamma, reco_e[m], bins=[gammaBins, eBins], weights=n)
        numGammaTheta += hist[0]

    return numGammaTheta


def main():
    #dn_dz = loadFlux(nuType, ch, m)[nRun::SKIP]
    
    #monteCarlo = loadMC(mcFile, nuType)
    #truncateMC(monteCarlo, nRun)
    numGammaTheta = gammaCalc(dn_dz)
    np.save("%s/e_d_theta_hist/partial_hists/ch%d_m%d_f%f_%s_%s_%d_energy_delta_theta_hist_im_gonna_scream_test.npy" % (data_path, ch, m, rescale_factor, binning, nuType, nRun), numGammaTheta)
    print("bang")


main()
