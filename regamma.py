"""
Read and MC file and write out the e_d_theta stuff fun. 
"""
import os

import numpy as np
import h5py

import gen_rescale_az_zen as gaz
import weight_MC as wmc
import solar_position_calc as sc
from controls import *

def strip_name(filepath): 
	return os.path.splitext(mcfile)[0].split("/")[-1].replace("_", "-")

class GammaCalcalater: 

	def __init__(self, mc_filepath, flux_filepath, lic_filepath, skip=1):

		self.skip = skip
		
		self.mc_filepath = mc_filepath
		self.flux_filepath = flux_filepath
		self.lic_filepath = lic_filepath

		self.mc_filename = strip_name(mc_filepath)
		self.flux_filename = strip_name(flux_filepath)

		self.check_comp()

		this.mc   = h5py.File(this.mcfilepath, "r")
		this.one_weights = self.calc_one_weights()
		self.dn_dz = np.load(args.fluxfile)[::SKIP]

	def calc_e_d_theta(self):
		slc = slice(:, :, self.skip)
		nu_e   = mc.root.NuEnergy[:]["value"][slc]
		nu_zen = mc.root.NuZenith[:]["value"][slc]
		nu_az  = mc.root.NuAzimuth[:]["value"][slc]
		reco_e = mc.root.MuExEnergy[:]["value"][slc]
		reco_zen = mc.root.MuExZenith[:]["value"][slc]
		reco_az  = mc.root.MuExAzimuth[:]["value"][slc]

		numGammaTheta = np.zeros((len(gammaBins)-1, len(eBins)-1))
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

		return numGammaTheta*self.skip

	def save(self, num_gamma): 
		np.save(
			f"/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/{self.flux_filename}_{self.mc_filename}_e_d_theta.npy",
			num_gamma.T[::-1]
		)

	def check_comp(self):
		assert(self.mc_filename==self.flux_filepath.split("/")[-1].split("_")[1][:-4])

	def calc_one_weights(self):
		try:
			ow = self.mc["oneweight"]["value"][::self.skip]*1e-4
		except:
			ow = wmc.weight_mc(self.lic_file)[::self.skip]*1.e-4
		return ow

	def __call__(self):
		SEED = abs(hash(flux_file)) % (2**32)
		np.random.seed(SEED)
		
		num_gamma =	self.calc_gamma(dn_dz, mc, fluxfile)
		self.save(num_gamma)

