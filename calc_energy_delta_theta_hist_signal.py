import numpy as np
import solar_position_calc as sc
from sys import argv as args
import mc
import dicts

ms  = [200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]
chs = [5,8,11]

SKIP = 100

np.random.seed(73)
ri  = np.random.randint(100000, size=(len(chs), len(ms)))
SEED_DICT = {ch:{} for ch in chs}
for i, ch in enumerate(chs):
    for j, mChi in enumerate(ms):
        SEED_DICT[ch][mChi] = ri[i][j]

nRun = int(args[1])
nuType = args[2]
ch = int(args[3])
m = int(args[4])
mcFile = args[5]

#BB_SEEDS  = {200 :31415, 300 :92653, 400 :27302, 500 :89793, 600 :93955, 700 :65695, 800 :59322, 900 :35163, 
#             1000:23846, 2000:55840, 3000:26433, 4000:96382, 5000:83279, 6000:, 10000:50288}
#WW_SEEDS  = {200 :7816,  300 :40628, 400 :54732, 500 :62089, 600 :60526, 700 :30218, 800 :65888, 900 :54407,
#             1000:98628, 2000:25176, 3000:3482 , 4000:49147, 5000:53421, 6000:, 10000:17067}
#TT_SEEDS  = {200 :95505, 300 :82231, 400 :47256, 500 :72535, 600 :18075, 700 :45910, 800 :30851, 900 :69113, 
#             1000 94081, 2000:97948, 3000:28481, 4000:65691, 5000:11745, 6000:, 10000:2841}
#SEED_DICT = {5:BB_SEEDS, 8:WW_SEEDS, 11:TT_SEEDS}

if nuType=="nu":
    SEED = SEED_DICT[ch][m] + nRun
else:
    SEED = SEED_DICT[ch][m] - 1 - nRun

np.random.seed(SEED)

def set_data_path():
    import re
    import os
    global data_path
    r = re.compile('cobalt.*.icecube.wisc.edu')
    if r.match(os.popen('hostname').readline().rstrip("\n")) is not None:
        data_path = "/data/user/jlazar/solar_WIMP/data/"
    elif os.popen('hostname').readline().rstrip("\n")=='MBP-FVFXC6EKHV2D.local':
        data_path = "/Users/jlazar/Documents/IceCube/solar_WIMP/data/"
    else:
        print("Machine not recognized")
        quit()

set_data_path()

#dataPath    = "/data/user/jlazar/solar_WIMP/data"
#dataPath    = "/Users/jlazar/Documents/IceCube/data"
rSun        = 6.9e10  # radius of sun in cm
solarZenPdf = np.loadtxt("%s/solar_zenith_pdf.txt" % data_path)
z           = np.radians(np.linspace(0, 179, 1800))
#sd_xs       = dicts.xs_limit_dict["sd"][ch][m]
#sd_ann_rate = dicts.sd_ann_rate_dict[m]
#rate        = sd_xs * sd_ann_rate
#rate        = sd_ann_rate

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


def loadFlux(nuType, ch, m):
    mc_recarray = np.load("%s/mcRecarray.npy" % data_path)
    if nuType == "nu":
        i  = np.where(mc_recarray["i"]==14)[0]
        #return np.load("%s/ch%d_m%d_nu_mu_flux.npy" % (dataPath, ch, m))
    elif nuType == "nuBar":
        #return np.load("%s/ch%d_m%d_nu_mu_bar_flux.npy" % (dataPath, ch, m))
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


def gammaCalc(dn_dz, monteCarlo):
    numGammaTheta    = np.zeros((len(gammaBins)-1, len(eBins)-1))
    for i, jd in enumerate(jds):
        x = sc.nParameter(jd)
        obl = sc.solarObliquity(x)
        G = sc.g(x)
        rad = sc.solarR(G)
        solar_solid_angle = 2*np.pi*(1-np.cos(rSun/rad))

        # gammaCut = np.arctan(rProductionRegion / rad)
        gammaCut = np.arctan(rSun / rad)

        monteCarlo.setTrueGamma(zeniths[i], azimuths[i])
        monteCarlo.setRecoGamma(zeniths[i], azimuths[i])
        #if i==0:
        #    print("max(dn_dz)==)" + str(max(dn_dz)))
        #    print("m=="+str(m))
        #    print("solar_solid_angle=="+str(solar_solid_angle))
        #    print("rad**2=="+str(rad**2))
        #    print("deltaT=="+str(deltaT))
        n = np.where(monteCarlo.trueGamma <= gammaCut,
                     dn_dz * monteCarlo.oneWeight,
                     #(1./m) * dn_dz * rate * (1./solar_solid_angle) * (1./(4*np.pi*rad**2)) * monteCarlo.oneWeight * deltaT,
                     0
                    )
        
        

        hist = np.histogram2d(monteCarlo.recoGamma, monteCarlo.recoE, bins=[gammaBins, eBins], weights=n)
        numGammaTheta += hist[0]

    return numGammaTheta


def main():
#    nRun = int(args[1])
#    nuType = args[2]
#    ch = int(args[3])
#    m = int(args[4])
#    mcFile = args[5]
    dn_dz = loadFlux(nuType, ch, m)[nRun::SKIP]
    
    monteCarlo = loadMC(mcFile, nuType)
    truncateMC(monteCarlo, nRun)
    numGammaTheta = gammaCalc(dn_dz, monteCarlo)
    np.save("%s/e_d_theta_hist/partial_hists/ch%d_m%d_%s_%d_energy_delta_theta_hist.npy" % (data_path, ch, m, nuType, nRun), numGammaTheta)


main()
