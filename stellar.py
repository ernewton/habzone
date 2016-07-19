# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 01:02:25 2016

@author: ellie
"""

import scipy.interpolate as interpolate
import numpy as np

def teff2mass(standard='baraffe'):
    '''Get mass from effective temperature by interpolating stellar models'''
    
    # teff, mass
    logt = False
    if standard == 'mamajek':
        teff, mass = np.genfromtxt('isochrones/mamajek_standards.dat', usecols=(1,15), unpack=True)
    elif standard == 'dartmouth':
        teff, mass = np.genfromtxt('isochrones/dartmouth_2Gyr_solar.dat', usecols=(2,1), unpack=True)
        logt = True
    elif standard == 'parsec':
        teff, mass = np.genfromtxt('isochrones/parsec_2Gyr_solar.dat',usecols=(5,3), unpack=True)
        logt = True
    elif standard == 'baraffe':
        teff, mass = np.genfromtxt('isochrones/baraffe15_2Gyr.dat',usecols=(1,0), unpack=True)
    else:
        return False

    if logt:
        teff = 10.**teff

    return interpolate.interp1d(teff, mass, bounds_error=False, kind='cubic')

def teff2lum(standard='baraffe'):
    '''Get lum from effective temperature by interpolating stellar models'''
    
    # teff, lum
    logt = False
    logl = False
    if standard == 'mamajek':
        teff, lum = np.genfromtxt('isochrones/mamajek_standards.dat', usecols=(1,5), unpack=True)
        logl = True
    elif standard == 'dartmouth':
        teff, lum = np.genfromtxt('isochrones/dartmouth_2Gyr_solar.dat', usecols=(2,4), unpack=True)
        logt = True
        logl = True
    elif standard == 'parsec':
        teff, lum = np.genfromtxt('isochrones/parsec_2Gyr_solar.dat',usecols=(5,4), unpack=True)
        logt = True
        logl = True
    elif standard == 'baraffe':
        teff, lum = np.genfromtxt('isochrones/baraffe15_2Gyr.dat',usecols=(1,2), unpack=True)
        logl = True
    else:
        return False

    if logl:
        lum = 10.**lum
    if logt:
        teff = 10.**teff
        
    return interpolate.interp1d(teff, lum, bounds_error=False, kind='cubic')


def from_mass(mymass, standard='baraffe'):    
    '''Get teff and lum from mass'''
    
    # teff, lum
    logt = False
    logl = False
    if standard == 'mamajek':
        teff, lum, mass = np.genfromtxt('isochrones/mamajek_standards.dat', 
                                        usecols=(1,5,15), unpack=True)
        logl = True
    elif standard == 'dartmouth':
        teff, lum, mass = np.genfromtxt('isochrones/dartmouth_2Gyr_solar.dat', 
                                  usecols=(2,4,1), unpack=True)
        logt = True
        logl = True
    elif standard == 'parsec':
        teff, lum, mass = np.genfromtxt('isochrones/parsec_2Gyr_solar.dat',
                                  usecols=(5,4,3), unpack=True)
        logt = True
        logl = True
    elif standard == 'baraffe':
        mass, teff, lum = np.genfromtxt('isochrones/baraffe15_2Gyr.dat',
                                        usecols=(0,1,2), unpack=True)
        logl = True
    else:
        return False

    if logl:
        lum = 10.**lum
    if logt:
        teff = 10.**teff
    
    pippo = interpolate.interp1d(mass, teff, bounds_error=False, kind='cubic')
    myteff = pippo(mymass)
    pippo = interpolate.interp1d(mass, lum, bounds_error=False, kind='cubic')
    mylum = pippo(mymass)
    
    return myteff, mylum
