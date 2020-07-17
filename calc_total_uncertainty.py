import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import os
import argparse
import emcee

import nuSQUIDSpy as nsq

import config
from physicsconstants import PhysicsConstants
from uncertainty_controls import ee, CC_central, CC_err_up, CC_err_down, qr_ch_dict 

def initialize_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ch",
                        type=int,
                        help="WIMPSim channel number bb:5, WW:8, tautau:11"
                       )
    parser.add_argument("-m",
                        type=int,
                        help="Dark matter mass"
                       )
    parser.add_argument('--ordering',
                        type=str,
                        help='normal (no) or inverted (io) ordering'
                       )
    parser.add_argument('-n',
                        type=int,
                        help='number of things to simulate. Default is 100',
                        default=100
                       )
    args = parser.parse_args()
    return args

def assymmetric_gaussian(x, mu, sigma0, sigma1):
    norm = 2 / (np.sqrt(np.pi)*(sigma0+sigma1))
    if x <= mu:
        val = norm * np.exp(-np.power(x-mu,2)/np.power(sigma0, 2))
    else:
        val = norm * np.exp(-np.power(x-mu,2)/np.power(sigma1, 2))
    return val

def get_emcee_params(ordering, ):
    if ordering=='no':
        mus    = np.asarray([33.82, 48.3, 8.61, 222, 7.39e-5, 2.523e-3, 1])
        sigmas = np.array([(0.76, 0.78), (1.9,1.1), (0.13,0.13), (28,38), (0.20e-5,0.21e-5), (0.030e-3,0.032e-3), (0.02, 0.038)]).T
    elif ordering=='io':
        mus    = np.array([33.82, 48.6, 8.65, 285, 7.39e-5, -2.509e-3, 1])
        sigmas = np.array([(0.76, 0.78), (1.5,1.1), (0.12,0.13), (26,24), (0.20e-5,0.21e-5), (0.030e-3,0.032e-3), (0.074,0.133)]).T
    else:
        print('invalid ordering option. Must be "no" or "io"')
        quit()
    return mus, sigmas

def lnprob(theta, mus, sigma0s, sigma1s, ):
    vals = np.zeros(len(theta), dtype=float)
    tups = [(theta[i], mus[i], sigma0s[i], sigma1s[i]) for i in range(len(theta))]
    for i, tup in enumerate(tups):
        vals[i] = np.log(assymmetric_gaussian(*tup))
    return np.sum(vals[np.isfinite(vals)])

def get_osc_params(mus, sigmas, nwalkers=1000, nsteps=500):
    ndim = len(mus)
    pos = [mus + 1e-3*mus*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mus, sigmas[0], sigmas[1]))
    sampler.run_mcmc(pos, nsteps)
    chain  = sampler.chain
    osc_params = chain[:,-1,:]
    return osc_params

def calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, xsec, e_min, nodes, param):
    e_max = m
    dn_dz = np.zeros((2, nodes))
    f = config.NuFlux("Pythia", e_min, e_max, nodes, qr_ch_dict[ch], m, param,
                           theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                           delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                           interactions=True, xsec=xsec, location='Earth', angle=0., latitude=0.)
    nu_mu_dn_dz     = np.asarray([tup[2] for tup in f]) * float(m)
    nu_mu_bar_dn_dz = np.asarray([tup[5] for tup in f]) * float(m)
    dn_dz[0][:] = nu_mu_dn_dz
    dn_dz[1][:] = nu_mu_bar_dn_dz
    return dn_dz

def main(ch, m, ordering, nwalkers, param):
    savedir = '/data/user/jlazar/solar_WIMP/data/total_uncertainties/'
    mus, sigmas = get_emcee_params(ordering)
    osc_params = get_osc_params(mus, sigmas, nwalkers=nwalkers)
    for osc_param in osc_params:
        print(tuple(osc_param))
        th12, th23, th13, delta, m12, m13, xs_scale = tuple(osc_param)
        xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5', ee, [xs_scale])
        savefname = '%d_%d_%s_%f_%f_%f_%f_%f_%f.npy' % (ch, m, ordering, th12, th23, th13, m12, m13, delta)
        print(savefname)
        flux = calc_flux(ch, m, th12, th23, th13, m12, m13, delta, xsec, 10, 200, param)
        np.save("%s/%s" % (savedir, savefname), flux)

if __name__=='__main__':
    param = PhysicsConstants()
    args = initialize_parser()
    ch        = args.ch
    m         = args.m
    ordering  = args.ordering
    nwalkers  = args.n
    seed = (hash(str(ch)+str(m)+ordering)+2) % 2**32
    np.random.seed(seed)
    main(ch, m, ordering, nwalkers, param)
