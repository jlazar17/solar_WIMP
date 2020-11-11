import numpy as np
import os
from os import sys
sys.path.append('/data/user/jlazar/solar_WIMP/')

from DM import DMAnnihilationJungmanSD

from scipy.special import gamma
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize
from scipy.stats import chi2
from scipy.integrate import quad

def parse_bool(thing):
    if thing in [True, '1', 1, 'yes', 'y']:
        return True
    elif thing in [False, '0', 0, 'no', 'n']:
        return False
    else:
        print('Input:%s not recognized' % thing)

def initialize_argparse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ch', type=int)
    parser.add_argument('-m', type=int)
    parser.add_argument('--opts', type=str, default='01')
    parser.add_argument('--savedir', type=str,
                        default='/data/user/jlazar/solar_WIMP/data/realizations/'
                       )
    args = parser.parse_args()
    return args

def factorial(x):
    x = np.asarray(x)
    return gamma(x+1)

class poisson():
    
    def __init__(self, mu):
        self.mu      = mu
        self.var     = mu
        self.std     = np.sqrt(mu)
        self.mean    = mu
        self.int_val = quad(self.val, 0, 10*self.std+mu)[0]
    
    def val(self, x):
        return np.exp(-self.mu) * np.power(self.mu, x) / factorial(x)
        
    def pdf(self, x):
        return self.val(x) / self.int_val
    
    def cdf(self, x):
        q = quad(self.pdf, 0, x)
        return q[0]

def log_likelihood(theta, mus, mub, data):
    logns, lognb = theta
    mu = np.exp(logns)*mus + np.exp(lognb)*mub
    llh = np.log(np.exp(-mu)*np.power(mu, data)/factorial(data))
    llh = llh[np.where(np.isfinite(llh))]
    return np.sum(llh)

class TS:
    
    def __init__(self, data, signal_dist, bg_dist, lxs_ini=0):
        
        self.data          = data
        self.signal_dist   = signal_dist
        self.bg_dist       = bg_dist
        self.lxs_ini       = lxs_ini
        self.null_bf_llh   = None
        self.signal_bf_llh = None
        self.bnds          = ((None, None), (None, None))
        
        self.bg_nll      = lambda x: -log_likelihood((-np.inf, x), 0, self.bg_dist, self.data)
        self.signal_nll  = lambda x: -log_likelihood(x, self.signal_dist, self.bg_dist, self.data)
        
    def set_null_bf_llh(self):
        self.bg_fit      = minimize(self.bg_nll, 0, method='BFGS', bounds=self.bnds)
        self.null_bf_llh = log_likelihood((np.log(0), self.bg_fit.x), 0, self.bg_dist, self.data)
        
    def set_signal_bf_llh(self):
        self.signal_fit    = minimize(self.signal_nll, (self.lxs_ini, 0), method='BFGS', bounds=self.bnds)
        self.signal_bf_llh = log_likelihood(self.signal_fit.x, self.signal_dist, self.bg_dist, self.data)
        
    def get_TS(self):
        if self.signal_bf_llh is None:
            self.set_signal_bf_llh()
        if self.null_bf_llh is None:
            self.set_null_bf_llh()
        return -2*(self.null_bf_llh-self.signal_bf_llh)

def get_distributions(ch, m, xs, scramble):
    if scramble:
        mu_bg    = np.load('/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/conv-numu_Nominal_e_d_theta_01.npy')
        mu_s_scr = DMAnnihilationJungmanSD(m,1e-39)/float(m)*np.load('/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/ch%d-m%d_Nominal_e_d_theta_01.npy' % (ch, m))
        mu_s     = DMAnnihilationJungmanSD(m,1e-39)/float(m)*np.load('/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/ch%d-m%d_Nominal_e_d_theta_00.npy' % (ch,m))

        bg_tmpl  = mu_bg+xs*mu_s_scr
        sig_tmpl = mu_s-mu_s_scr
    else:
        mu_bg    = np.load('/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/conv-numu_Nominal_e_d_theta_00.npy')
        mu_s     = DMAnnihilationJungmanSD(m,1e-39)/float(m)*np.load('/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/ch%d-m%d_Nominal_e_d_theta_00.npy' % (ch,m))

        bg_tmpl  = mu_bg
        sig_tmpl = mu_s
    return mu_bg, mu_s, bg_tmpl, sig_tmpl

def main(ch, m, savedir,  opts, n=10000):
    if not os.path.exists(savedir):
        os.makedir(savedir)
    scramble = parse_bool(opts[1])
    xss = np.logspace(-2, 2, 30)
    results = np.recarray((len(xss)*n),
                          dtype=[('bg_ts', float), ('sig_ts', float), ('inj_xs', float), ('fit_xs', float)]
                         )
    
    for i, xs in enumerate(xss):
        mu_bg, mu_s, bg_tmpl, sig_tmpl = get_distributions(ch, m, xs, scramble)

        slc       = slice(i*n, (i+1)*n)
        null_TS   = np.zeros(n)
        signal_TS = np.zeros(n)
        fit_xs    = np.zeros(n)
        for i in range(n):
            lxs_ini      = np.log(10)*(np.log10(xs)+0.5*(np.random.rand()-0.5))
            null_data    = np.random.poisson(mu_bg)
            ts           = TS(null_data, sig_tmpl, bg_tmpl, lxs_ini=lxs_ini)
            null_TS[i]   = ts.get_TS()
            inj_data     = np.random.poisson(mu_bg+xs*mu_s)
            ts           = TS(inj_data, sig_tmpl, bg_tmpl, lxs_ini=lxs_ini)
            signal_TS[i] = ts.get_TS()
            fit_xs[i]    = np.exp(ts.signal_fit.x[0])*1e-39
        results['bg_ts'][slc]  = null_TS
        results['sig_ts'][slc] = signal_TS
        results['inj_xs'][slc] = xs*1e-39
        results['fit_xs'][slc] = fit_xs
    np.save('%s/ch%d-m%d_%s_realiztions' % (savedir, ch, m, opts), results)

if __name__=='__main__':
    args = initialize_argparse()
    main(args.ch, args.m, args.savedir, args.opts, n=1000)
