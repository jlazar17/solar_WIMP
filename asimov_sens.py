
# coding: utf-8

# In[2]:


import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
plt.style.use("/data/user/jlazar/solar_WIMP/paper.mplstyle")

from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize
from scipy.optimize import fmin
from scipy.special import gamma
from scipy.optimize import curve_fit
from scipy.stats import poisson

import dicts

import interp_e_d_theta as edt

import sys

args = sys.argv

# In[3]:


sd_anns   = dicts.sd_ann_rate_dict # annihilation rate for xs=1pb
ms        = [300, 400, 500, 800, 1000, 2000, 3000, 4000, 5000, 8000, 10000]

nice_ch_strs = {5  : r"$b\bar{b}$",
                8  : r"$W^{+}W^{-}$",
                11 : r"$\tau^{+}\tau^{-}$"
               }

# In[4]:


m, r = zip(*sorted(sd_anns.items()))
m, r = np.asarray(m), np.asarray(r)
f = UnivariateSpline(np.log10(m), np.log10(r),)
def get_ann_rate(m):
    lm = np.log10(m)
    return np.power(10, f(lm))


# In[5]:


def round_sig_figs(x, sig_figs=2):
    return round(x, sig_figs-int(np.floor(np.log10(abs(x))))-1)


# In[6]:


def factorial(x):
    x = np.asarray(x)
    return gamma(x+1)


# In[7]:


def poisson_pdf(k, mu):
    return np.exp(-mu) * np.power(mu, k) / factorial(k)


# In[8]:


def neg_llh(mu_s, n_s, mu_b, n_b):
    """
    mu_s float: mean number of signal events
    n_s  int  : observed number of signal events
    mu_b float: mean number of background events
    n_b  int  : observed number of bacground events
    """
    mu_s = np.asarray(mu_s)
    mu_b = np.asarray(mu_b)
    n_b  = np.asarray(n_b)
    n_s  = np.asarray(n_s)
    mu = mu_s*n_s + mu_b*n_b
    lh = -np.log(np.exp(-mu) * np.power(mu, mu_b) / factorial(mu_b))
#     lh = lh[np.where(~np.isnan(lh))]
#     lh = lh[np.where(~np.isinf(lh))]
    return np.sum(lh)


# In[9]:


def model(mu_s, n_s, mu_b, n_b):
    mu_s = np.asarray(mu_s)
    mu_b = np.asarray(mu_b)
    n_b  = np.asarray(n_b)
    n_s  = np.asarray(n_s)
    mu = mu_s*n_s + mu_b*n_b
    return np.exp(-mu) * np.power(mu, mu_b) / factorial(mu_b)


# In[17]:


def MLERegression(params, ch, dp):
    """
    Find negative log-likelihood value for input params given input signal
    
    Parameters
    ----------
    params : tuple
        (WIMP mass in GeV, WIMP-nucleon cross section in picobarns, normalization on background)
    ch     : int
        WIMPSim channel number (bb:5, WW:8, tautau:11)
    dp     : np.ndarray
        Injected distribution of events for 7 years of data taking in E_{reco} and \Delta\theta_{reco}
    
    Returns
    -------
    sum_neg_llh : float
        neg-llh value for input params given injected event distribution
    """
    m, xs, nb = params[0], params[1], params[2]
    mu_s = make_signal(ch, m, xs)
    mu_b = np.load("/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/bg_tot_e_d_theta_hist.npy")
    yhat = (mu_s+nb*mu_b)
    neg_llh = -np.log(poisson_pdf(dp, yhat))
    sum_neg_llh = np.sum(neg_llh)
    return sum_neg_llh


# In[113]:


def make_signal(ch, m, xs):
    """
    Generate signal dist. in E_{reco} \Delta\theta_{reco} for given ann. channel, WIMP mass, and cross section

    Parameters
    ----------
    ch : int
        WIMPSim channel number (bb:5, WW:8, tautau:11)
    m  : int
        Mass of the WIMP in GeV
    xs : float
        WIMP-nucleon cross section in picobarns
    
    Returns
    -------
    mu_s : np.ndarray
        Distribution of expected number signal events for 7 years of data taking in E_{reco} and \Delta\theta_{reco}
    
    """
    delta_t = 30 * 60
    factor  = xs * (1. / m) * get_ann_rate(m) * 1e-3 * delta_t
    if m in ms:
        mu_s = factor*np.load("/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/ch%d_m%d_f1.000000_f_tot_e_d_theta_hist.npy" % (ch, m))
    else:
        m0   = ms[np.digitize(m, ms)-1]
        m1   = ms[np.digitize(m, ms)]
        x    = (m - m0) / (m1-m0)
        mu_s = factor*edt.interp_e_d_theta(ch, m0, m1, x)
    return mu_s


# In[124]:


def do_calc(ch, m_inj, xs_inj):
    sig_inj = make_signal(ch, m_inj, xs_inj)
    bg_inj  = np.load("/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/bg_tot_e_d_theta_hist.npy")
    data    = sig_inj + bg_inj
    
    xss = np.logspace(np.log10(xs_inj)-1.5, np.log10(xs_inj)+1.5, 51)
    #xss = np.logspace(np.log10(xs_inj)-1.5, np.log10(xs_inj)+1.5, 21)
    #mss = np.logspace(np.log10(550), 4, 21)
    mss = np.logspace(2.75, 4, 13)
    tups = [(m,x,1) for m in mss
                    for x in xss]
    outs = np.full(len(tups), 0.)
    for i, tup in enumerate(tups):
        outs[i] = MLERegression(tup, ch, data)
    return outs, mss, xss
    
    


# In[131]:
def plot(params, ch, m_inj, xs_inj):
    outs = params[0]
    mss  = params[1]
    xss  = params[2]

    plouts = outs.reshape((len(mss), len(xss)))[::-1].T
    plmin  = plouts.min()
    plouts = plouts - plmin + 0.000001

    X, Y   = np.log10(mss)[::-1], np.log10(xss)
    
    
    fig, ax = plt.subplots(1, 1)
    
    pcm = ax.pcolor(X, Y, plouts,
                     norm=colors.LogNorm(vmin=1e-2, vmax=plouts.max()),
                     cmap="PuBu",
                   )
    
    levels=[2.27886544, 2.9955]
    levels_dict = {2.27886544 : r"$1\sigma$",
                   2.9955     : "90\%"
                  }
    
    CS = ax.contourf(X, Y, plouts,
                     levels,
                     alpha=0
                    )
    CS2 = ax.contour(CS, levels=CS.levels, colors=("orange", "yellow"))
    CL  = ax.clabel(CS2, fmt=levels_dict, fontsize=25)
    
    ax.set_xlabel(r"$\log\left(\frac{m_{\chi}}{\rm{GeV}}\right)$")
    ax.set_ylabel(r"$\log\left(\frac{\sigma_{\chi N}}{\rm{fb}}\right)$")
    ax.set_title(nice_ch_strs[ch])
    
    ax.scatter(np.log10(m_inj), np.log10(xs_inj), color="r", marker="+", s=400)
    
    ax.legend(frameon=False)
    fig.colorbar(pcm, ax=ax, extend="min")
                    
    plt.savefig("/data/user/jlazar/solar_WIMP/plots/asimov_sens_ch%d_m%f_xs%f.png" % (ch, m_inj, xs_inj), dpi=400, bbox_inches="tight")

#def plot(params, ch, m_inj, xs_inj):
#    outs = params[0]
#    mss  = params[1]
#    xss  = params[2]
#    plouts = outs.reshape((len(mss), len(xss)))[::-1].T
#    plouts = plouts - plouts.min()+0.000001
#    
#    X, Y   = np.log10(mss)[::-1], np.log10(xss)
#    
#    fig, ax = plt.subplots(1, 1)
#    
#    pcm = ax.pcolor(X, Y, plouts,
#                     norm=colors.LogNorm(vmin=plouts.min(), vmax=plouts.max()),
#                     cmap="PuBu",
#                    )
#    
#    ax.set_xlabel(r"$\log\left(\frac{m_{\chi}}{GeV}\right)$")
#    ax.set_ylabel(r"$\log\left(\frac{\sigma_{\chi N}}{\rm{fb}}\right)$")
#    
#    ax.scatter(np.log10(m_inj), np.log10(xs_inj), color="r", marker="+", s=400)
#    fig.colorbar(pcm, ax=ax, extend='max')
#    
#    plt.savefig("plots/asimov_sens_ch%d_m%f_xs%f.png", dpi=400, bbox_inches="tight")
##    plt.show()

def main():
    ch     = int(args[1])
    m_inj  = float(args[2])
    xs_inj = float(args[3])
    plt_params = do_calc(ch, m_inj, xs_inj)
    plot(plt_params, ch, m_inj, xs_inj)

main()
