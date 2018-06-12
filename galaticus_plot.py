#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import h5py
import dtk

def logspace(data,num):
    return np.logspace(np.log10(np.min(data)),np.log10(np.max(data)),num)

hfile = h5py.File("test.hdf5",'r')


plt.figure()
H,x_bins,y_bins = np.histogram2d(hfile["gal_core/sdss_i"],hfile["gal_core/stellar_mass"],
                                 bins=(100,logspace(hfile["gal_core/stellar_mass"],100)))
plt.pcolormesh(x_bins,y_bins,H.T,cmap='BuPu',norm=clrs.LogNorm())
plt.yscale('log')
plt.figure()
plt.title("Core Satellite Galaxies")
slct = hfile["gal_core/parentIndex"][:]!=hfile["gal_core/nodeIndex"][:]
print np.sum(slct)
plt.plot(hfile["gal_core/sdss_i"][slct],hfile["gal_core/stellar_mass"][slct],'x',alpha=0.5)
plt.yscale('log')
plt.axhline(y=1e7,color='r')
plt.xlabel("SDSS i-band Mag")
plt.ylabel("Stellar Mass [Msun?]")
above = np.sum(hfile["gal_core/stellar_mass"][slct]>1e7)
below = hfile["gal_core/stellar_mass"][slct].size-above
plt.plot([],[],'r',label='1e7 stellar cut\nabove - %d\nbelow - %d'%(above,below))
plt.yscale('log')
plt.grid()
plt.legend(loc='best')


plt.figure()
plt.title("Halo Tag Particle Satellite Galaxies")
slct = hfile["gal_bhp/parentIndex"][:]!=hfile["gal_bhp/nodeIndex"][:]
plt.plot(hfile["gal_bhp/sdss_i"][slct],hfile["gal_bhp/stellar_mass"][slct],'x',alpha=0.5)
plt.axhline(y=1e7,color='r')
plt.xlabel("SDSS i-band Mag")
plt.ylabel("Stellar Mass [Msun?]")
above = np.sum(hfile["gal_bhp/stellar_mass"][slct]>1e7)
below = hfile["gal_bhp/stellar_mass"][slct].size-above
plt.plot([],[],'r',label='1e7 stellar cut\nabove - %d\nbelow - %d'%(above,below))
plt.yscale('log')
plt.grid()
plt.legend(loc='best')


plt.figure()
plt.title("Halo Particle Satellite Galaxies")
slct = hfile["gal_hp/parentIndex"][:]!=hfile["gal_hp/nodeIndex"][:]
plt.plot(hfile["gal_hp/sdss_i"][slct],hfile["gal_hp/stellar_mass"][slct],'x',alpha=0.5)
plt.axhline(y=1e7,color='r')
plt.xlabel("SDSS i-band Mag")
plt.ylabel("Stellar Mass [Msun?]")
above = np.sum(hfile["gal_hp/stellar_mass"][slct]>1e7)
below = hfile["gal_hp/stellar_mass"][slct].size-above
plt.plot([],[],'r',label='1e7 stellar cut\nabove - %d\nbelow - %d'%(above,below))
plt.yscale('log')
plt.grid()
plt.legend(loc='best')

plt.figure()
plt.title("Unmatched Satellite Galaxies")
slct = hfile["gal_blank/parentIndex"][:]!=hfile["gal_blank/nodeIndex"][:]
plt.plot(hfile["gal_blank/sdss_i"][slct],hfile["gal_blank/stellar_mass"][slct],'x',alpha=0.5)
plt.axhline(y=1e7,color='r')
plt.xlabel("SDSS i-band Mag")
plt.ylabel("Stellar Mass [Msun?]")
above = np.sum(hfile["gal_blank/stellar_mass"][slct]>1e7)
below = hfile["gal_blank/stellar_mass"][slct].size-above
plt.plot([],[],'r',label='1e7 stellar cut\nabove - %d\nbelow - %d'%(above,below))
plt.yscale('log')
plt.grid()
plt.legend(loc='best')
plt.figure()
plt.plot(hfile['gal_blank/host_halo_mass'],hfile['gal_blank/infall_halo_mass'])
plt.yscale('log')
plt.xscale('log')
print 


plt.show()
