#!/usr/bin/env python
# 
# Vertical integral function for calculating column water vapor (CWV) aka
# IWV aka PW.
# 
# James Ruppert  
# jruppert@ou.edu  
# March 2025

from thermo_functions import *
import numpy as np

def vert_integral_hydro(data):
    p = data['p']
    try:
        sh = data['sh']
    except:
        sh = mixr2sh(data['mr'])
    var_int = -np.trapezoid(sh, p, axis=1)/9.81
    # Replace NaNs
    nans = np.where(np.isnan(p[:,5]))
    var_int[nans] = np.nan
    return var_int

def vert_integral(data, zbase=0, set_nans=True):
    p = np.ma.filled(data['p'], np.nan)
    t = np.ma.filled(data['tmpk'], np.nan)
    hght = data['hght']
    if hght.ndim == 1:
        iz = np.where(hght >= zbase)[0]
        hght = np.repeat(hght[np.newaxis,:], p.shape[0], axis=0)
    else:
        iz = np.where(hght[0,:] >= zbase)[0]
    try:
        sh = data['sh']
        # mr = sh2mixr(sh)
    except:
        mr = data['mr']
        sh = mixr2sh(mr)
    # rho = density_moist(t, sh, p) # kg/m3
    rd=287.04
    rv=461.5
    eps_r=rv/rd
    # The two forms of rho below are identical, just using r vs q
    # rho = p / ( rd * t * (1. + mr*eps_r)/(1.+mr) )
    # Below is the form of Schulz and Stevens (2018)
    rho = p / ( rd * t * (1 + sh*(eps_r - 1)) )
    cwv = np.trapezoid(sh[:,iz] * rho[:,iz], hght[:,iz], axis=1)
    # cwv = np.nansum((rho[:,iz] * sh[:,iz] * np.gradient(hght[:,iz], axis=1)), axis=1)
    # Replace NaNs
    if set_nans:
        nans = np.where(np.isnan(p[:,5]))
        cwv[nans] = np.nan
    return cwv
