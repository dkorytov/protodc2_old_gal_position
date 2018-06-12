#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.special import hyp2f1
from scipy.integrate import quad

from util2 import *

def M_of_r(r):
    return 1.0/(r+1.0)+np.log(r+1.0)

def M_of_rg(r,g):
    a = -1.0/((g-1.0)*(r+1))*r**(1-g)
    b = (r+1.0)*hyp2f1(1.0-g,1.0-g,2.0-g,-r)
    c = (r+1.0)*hyp2f1(1.0-g,3.0-g,2.0-g,-r)
    d = -2.0*(x+1.0)**g
    return a*(b+c+d)

def NFW(r):
    return r**2/(r*(1.0+r)**2)

def NFW_g(r,g):
    return r**2/(r**g*(1.0+r)**(3-g))

def get_Mstar_func():
    data = np.loadtxt("data/mstar.txt",skiprows=1)
    print data.shape
    f = interp1d(data[:,0],np.log10(data[:,1]))
    def _f(z):
        log10M = f(z)
        return 10**log10M
    return _f

def rnd_unit():
    s = np.random.normal(0,1.0,3)
    mag = np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)
    return s[0]/mag,s[1]/mag,s[2]/mag


def integrate_func(f,r,args=()):
    vals = np.zeros(r.size)
    err = np.zeros(r.size)
    for i in range(0,r.size-1):
        y,abserr = quad(f,r[i],r[i+1],args=args)
        vals[i+1]=y
        err[i]=abserr
    return np.cumsum(vals)
_g = 0.0
_rs = np.logspace(-3,2.5,200)
_Ms = M_of_r(_rs)-1.0
_f_r_of_Ms = interp1d(_Ms,_rs,kind='quadratic')
_f_Ms_of_r = interp1d(_rs,_Ms,kind='quadratic')

def set_gamma(g):
    global _g,_Mgs, _f_r_of_Mgs, _f_Mgs_of_r
    _g = g
    _Mgs = integrate_func(NFW_g,_rs,args=(g))
    _f_r_of_Mgs = interp1d(_Mgs,_rs,kind='quadratic')
    _f_Mgs_of_r = interp1d(_rs,_Mgs,kind='quadratic')

def get_gamma():
    return _g

def get_m200(fof_mass,z):
    return so_from_fof_z(fof_mass,z)

def get_r200(fof_mass,z):
    sod_mass = so_from_fof_z(fof_mass,z)
    r200 = rdelta_from_sod(sod_mass,200.0,z=z)
    return r200

def rnd_NFWg_radius(max_nR200,conc):
    #note max
    max_rs = max_nR200*conc
    max_rs_mass = _f_Mgs_of_r(max_rs)
    mass = max_rs_mass*np.random.rand()
    rs = _f_r_of_Mgs(mass)
    return rs/conc
    
def rnd_NFW_radius(max_nR200,conc):
    max_rs = max_nR200*conc
    max_rs_mass = _f_Ms_of_r(max_rs)
    mass = max_rs_mass*np.random.rand()
    rs = _f_r_of_Ms(mass)
    return rs/conc



if __name__ == "__main__":
    set_gamma(1.0)
    data_x = []
    data_y = []
    data_z = []
    data_r = []
    for i in range(0,1000000):
    #r = r_of_M();
        #r=1.0
        r=rnd_NFWg_radius(3.0,6)
        x,y,z = rnd_unit()
        data_x.append(x*r)
        data_y.append(y*r)
        data_z.append(z*r)
        data_r.append(np.sqrt(x*x*r*r + y*y*r*r + z*z*r*r))

    plt.figure()
    plt.plot(data_x,data_y,',',alpha=0.5)
    plt.axis("equal")
    plt.grid()
    

    plt.figure()
    rs = np.logspace(-3,np.log10(3),20)
    rs_avg = (rs[:-1]+rs[1:])/2.0
    rs_area = 4.0*np.pi*(rs[1:]**2 - rs[:-1]**2)
    rs_vol = 4.0/3.0*np.pi * (rs[1:]**3 - rs[:-1]**3)
    H, x_bins = np.histogram(data_r,bins=rs)
    plt.plot(rs_avg,H/rs_vol/np.sum(H/rs_vol),'-x',label='random sample')
    H2 = NFW_g(rs_avg*6.0,get_gamma())/rs_area
    plt.plot(rs_avg,H2/np.sum(H2),'-o',label='theory NFWg')
    H3 = NFW(rs_avg*6.0)/rs_area
    plt.plot(rs_avg,H3/np.sum(H3),'--s',label='theory NFW')
    plt.yscale('log')
    plt.xscale('log')
    r=_rs
    
    plt.figure()
    plt.title("radial density profile")
    plt.plot(r,NFW(r),'--')
    plt.plot(r,NFW_g(r,1.0))
    plt.plot(r,NFW_g(r,1.1))
    plt.plot(r,NFW_g(r,1.5))
    plt.plot(r,NFW_g(r,2.0))
    plt.plot(r,NFW_g(r,2.5))
    plt.plot(r,NFW_g(r,3.0))
    plt.plot(r,NFW_g(r,4.0))
    plt.yscale('log')
    plt.xscale('log')

    plt.figure()
    plt.title("radial count profile")
    plt.plot(r,integrate_func(NFW,r))
    plt.plot(r,M_of_r(r)-1,'--')
    plt.plot(r,integrate_func(NFW_g,r,args=(1.0)))
    plt.plot(r,integrate_func(NFW_g,r,args=(1.1)))
    plt.plot(r,integrate_func(NFW_g,r,args=(1.5)))
    plt.plot(r,integrate_func(NFW_g,r,args=(2.0)))
    plt.plot(r,integrate_func(NFW_g,r,args=(2.5)))
    plt.plot(r,integrate_func(NFW_g,r,args=(3.0)))
    plt.plot(r,integrate_func(NFW_g,r,args=(4.0)))
    plt.yscale('log')
    plt.xscale('log')
    
    plt.show()
