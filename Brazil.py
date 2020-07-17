#!/usr/bin/env python

import sys, time, os, glob, json, random, subprocess
import GolemFitPy as gf
import numpy as np
from optparse import OptionParser
from os.path import expandvars
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import pylab
from pylab import *
import scipy
from scipy import stats
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import scipy.stats
from scipy.stats import chisquare
import scipy.ndimage
import warnings

import scipy.interpolate as interpolate

try:
    import psutil
    def getMem():
        return 'Ram: '+str(np.round((psutil.Process(os.getpid()).memory_info().rss)/1.E9,4))+'GB'
except:
    def getMem():
        return 'Ram: psutil not installed'
