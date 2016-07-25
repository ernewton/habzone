# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 01:01:13 2016

@author: ellie
"""

import numpy as np
import os
import stellar

def kep2dist(mstar, period, days=True):
    '''Semi-major axis from Keplerian orbital period'''
    if days:
        p = period/365.
    else:
        p = period
    a = ((p**2)*mstar)**(1./3.)
    return a

def kep2per(mstar, a):
    '''Keplerian orbital period from semi-major axis'''
    
    p = ((a**3)/mstar)**(1./2.)
    p *= 365.
    return p

def habzone(tstar, lstar, verbose=True, inner='moist', outer='maximum'):
    '''Calculate habitable zone insolation and semi-major axis'''

    # from web
    #hzfile = 'hz_coefficients.dat'
    #runaway = 1
    #moist = -1
    #maximum = 2
    
    # from paper erratum
    try:
        hzfile = os.path.join(os.path.abspath(os.path.dirname(__file__)))+'/hz_coefficients_updated.dat'
    except:
        hzfile='hz_coefficients_updated.dat'
    kop_coeff = np.genfromtxt(hzfile)
    runaway = 1
    moist = 2
    maximum = 3
    if inner=='moist':
        inn = 2
    elif inner=='runaway':
        inn = 1
    else:
        raise
    if outer=='maximum':
        out = 3
    else:
        raise
    
    tsun = 5780. # constant
    shp = len(kop_coeff[0])
    
    seffsun = kop_coeff[0]
    a = kop_coeff[1]
    b = kop_coeff[2]
    c = kop_coeff[3]
    d = kop_coeff[4]
    
    tdiff = tstar-tsun
    seff = np.zeros(shp)
    for i in np.arange(shp):
         seff[i] = seffsun[i] + a[i]*tdiff + b[i]*(tdiff**2) + c[i]*(tdiff**3) + d[i]*(tdiff**4)
    
    dist = np.sqrt(lstar/seff)
    
    if verbose:
        print "Runaway"
        print seff[runaway], dist[runaway]
        print "Maximum"
        print seff[maximum], dist[maximum]
        if moist > 0:
            print "Moist"
            print seff[moist], dist[moist]
    
    return seff[np.array([inn, out])], dist[np.array([inn, out])]


# supply stellar mass
# calculate habitable zone
# turn that into an orbital period
def habzone_m2p(mstar, verbose=False, standard='baraffe'):
    '''Calcualte habitable zone period from stellar mass'''
    tstar, lstar = stellar.from_mass(mstar, standard=standard)
    seff, dist = habzone(tstar, lstar, verbose=verbose)
    inner = kep2per(mstar, dist[0])
    outer = kep2per(mstar, dist[1])
    if verbose:
        print "Mstar = ", mstar
        print "-> Teff = ", tstar
        print "-> Lum = ", lstar
        print "Habitable zone insolation in Earth flux and distance in AU"
        print inner, outer
    # orbital period corresponding to the inner and outer edge of the habitable zone
    return inner, outer
