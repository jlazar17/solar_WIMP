import numpy as np
import solar_position_calc as sc
from sys import argv as args
import mc
import dicts
import os

ms  = [200,300,400,500,600,700,800,900,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000]
chs = [5,8,11]

SKIP = 100

np.random.seed(74)
ri  = np.random.randint(100000, size=(len(chs), len(ms)))
SEED_DICT = {ch:{} for ch in chs}
for i, ch in enumerate(chs):
    for j, mChi in enumerate(ms):
        SEED_DICT[ch][mChi] = ri[i][j]

nRun = int(args[1])
nuType = args[2]
ch = int(args[3])
m = int(args[4])

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
azimuths = np.random.rand(len(jds)) * 2 * np.pi
zeniths = np.random.choice(z, len(jds), p=solarZenPdf)

binsPerDecade = 10
highExp       = 6.5
lowExp        = 0.5
eBins         = np.logspace(lowExp, highExp, int(binsPerDecade * (highExp-lowExp) + 1))
gammaBins     = np.linspace(0, np.pi, 113)# 113 picked to best match binning on MEOWS wiki page


####### LOAD IN MC INFORMATION ######
mc_arr  = np.load("%s/mcRecarray.npy" % data_path)
if nuType=="nu":
    indices = np.where(mc_arr["i"]==14)[0]
elif nuType=="nuBar":
    indices = np.where(mc_arr["i"]==-14)[0]
else:
    print("wrongo !")
    quit()
ow      = mc_arr["oneWeight"][indices][nRun::SKIP]
# true quantities
nu_e   = mc_arr["nuE"][indices][nRun::SKIP]
nu_zen = mc_arr["nuZen"][indices][nRun::SKIP]
nu_az  = mc_arr["nuAz"][indices][nRun::SKIP]
# reco quantities
reco_e = mc_arr["recoE"][indices][nRun::SKIP]
arr     = np.load("%s/rescaled_reco_angles.npy" % data_path)
reco_zen = arr["rescaled_reco_zen"][indices][nRun::SKIP]
reco_az  = arr["rescaled_reco_az"][indices][nRun::SKIP]

dn_dz = np.load("%s/mc_dn_dz/ch%d_m%d_mc_dn_dz.npy" % (data_path, ch, m))[indices][nRun::SKIP]

def opening_angle(zen1, az1, zen2, az2):
    return np.arccos(np.sin(zen1)*np.sin(zen2)*np.cos(az1-az2)+np.cos(zen1)*np.cos(zen2))

def loadFlux(nuType, ch, m):
    mc_recarray = np.load("%s/mcRecarray.npy" % data_path)
    if nuType == "nu":
        i  = np.where(mc_recarray["i"]==14)[0]
    elif nuType == "nuBar":
        i  = np.where(mc_recarray["i"]==-14)[0]
    return np.load("%s/mc_dn_dz/ch%d_m%d_mc_dn_dz.npy" % (data_path, ch, m))[i]


def loadMC(mc_rec_array_file, nuType):
    return mc.MonteCarlo(mcFile, nuType)


def truncateMC(monteCarlo, nRun):
    monteCarlo.setNuZen(monteCarlo.nuZen[nRun::SKIP])
    monteCarlo.setNuAz(monteCarlo.nuAz[nRun::SKIP])
    monteCarlo.setNuE(monteCarlo.nuE[nRun::SKIP])
    monteCarlo.setRecoZen(monteCarlo.recoZen[nRun::SKIP])
    monteCarlo.setRecoAz(monteCarlo.recoAz[nRun::SKIP])
    monteCarlo.setRecoE(monteCarlo.recoE[nRun::SKIP])
    monteCarlo.setOneWeight(monteCarlo.oneWeight[nRun::SKIP])


def gammaCalc(dn_dz):
    numGammaTheta    = np.zeros((len(gammaBins)-1, len(eBins)-1))
    print("boing")
    for i, jd in enumerate(jds):
        x = sc.nParameter(jd)
        obl = sc.solarObliquity(x)
        G = sc.g(x)
        rad = sc.solarR(G)
        solar_solid_angle = 2*np.pi*(1-np.cos(rSun/rad))

        gammaCut = np.arctan(rSun / rad)

        nu_gamma   = opening_angle(nu_zen, nu_az, zeniths[i], azimuths[i])
        reco_gamma = opening_angle(reco_zen, reco_az, zeniths[i], azimuths[i])
        #monteCarlo.setTrueGamma(zeniths[i], azimuths[i])
        #monteCarlo.setRecoGamma(zeniths[i], azimuths[i])
        #if i==0:
        #    print("max(dn_dz)==)" + str(max(dn_dz)))
        #    print("m=="+str(m))
        #    print("solar_solid_angle=="+str(solar_solid_angle))
        #    print("rad**2=="+str(rad**2))
        #    print("deltaT=="+str(deltaT))
        n = np.where(nu_gamma<= gammaCut,
                     dn_dz * ow * \
                     (1. / solar_solid_angle) * \
                     (1. / (4*np.pi*np.power(rad, 2))),
                     0
                    )
        
        

        hist = np.histogram2d(reco_gamma, reco_e, bins=[gammaBins, eBins], weights=n)
        numGammaTheta += hist[0]

    return numGammaTheta


def main():
#    nRun = int(args[1])
#    nuType = args[2]
#    ch = int(args[3])
#    m = int(args[4])
#    mcFile = args[5]
#    dn_dz = loadFlux(nuType, ch, m)[nRun::SKIP]
    
#    monteCarlo = loadMC(mcFile, nuType)
#    truncateMC(monteCarlo, nRun)
    numGammaTheta = gammaCalc(dn_dz)
#    if ~os.path.isdir("%s/e_d_theta_hist/partial_hists" % data_path):
#        os.mkdir("%s/e_d_theta_hist/partial_hists" % data_path)
    np.save("%s/e_d_theta_hist/partial_hists/ch%d_m%d_%s_%d_energy_delta_theta_hist_rescale_%d.npy" % (data_path, ch, m, nuType, nRun, plz), numGammaTheta)


main()
