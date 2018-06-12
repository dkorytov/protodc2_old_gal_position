#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import h5py 
import sys
import dtk

param = dtk.Param(sys.argv[1])

mt_loc = param.get_string("mt_loc")
sod_loc = param.get_string("sod_loc")

subfile_num = param.get_int('subfile_num')
step = param.get_int("step")
output = param.get_string("output").replace("${step}",str(step))

mt_cat = dtk.Catalog(mt_loc)
mt_cat.add_step(step)
for i in range(0,subfile_num):
    mt_cat.add_subfile(i)
mt_cat.add_var("/forestHalos/timestep",as_name='time_step')
mt_cat.add_var("/forestHalos/haloTag",as_name='htag')
mt_cat.read_hdf5()
mt_cat.apply_function('htag',dtk.frag_to_real)
slct = mt_cat[step]["time_step"] == step
mt_cat[step]["time_step"]=mt_cat[step]["time_step"][slct]
mt_cat[step]["htag"]=mt_cat[step]["htag"][slct]

sod_cat = dtk.Catalog(sod_loc)
sod_cat.add_step(step)
sod_cat.add_var("fof_halo_tag",as_name='htag')
sod_cat.add_var("sod_halo_mass")
sod_cat.add_var("sod_halo_radius")
sod_cat.add_var("sod_halo_cdelta")
sod_cat.read_gio()
sod_mass = sod_cat[step]["sod_halo_mass"]
sod_htag = sod_cat[step]["htag"]
srt = np.argsort(sod_htag)
sod_htag = sod_htag[srt]
sod_mass = sod_mass[srt]


#Plotting Rs as a funciton of M200
##########################################
#H,x_bins,y_bins = np.histogram2d(sod_cat[step]["sod_halo_mass"],sod_cat[step]["sod_halo_radius"]/sod_cat[step]["sod_halo_cdelta"],
#                                 bins = (np.logspace(11,16,100),np.logspace(-3,1,100)))

#plt.figure()
#plt.pcolormesh(x_bins,y_bins,H.T,cmap="PuBu",norm=clrs.LogNorm())
#plt.xscale('log')
#plt.xlabel("sod mass [Msun/h]")
#plt.yscale('log')
#plt.ylabel("Rs [Mpc/h]")
#plt.grid()
#plt.show()

halo_cat = dtk.Catalog()
halo_cat.join(mt_cat,sod_cat,join_on="htag")

def schechter_func(mr,m200,M,a,A,g):
    #TODO: Mass has no h in it?
    m200 = m200*3.0/0.7
    Nb = A*(m200/1e12)**g
    F = 10**(0.4*(M-mr)*(a+1.0))*np.exp(-10**(0.4*(M-mr)))
    return Nb*F

def global_lum_func(mr,m200):
    #Ting-Wen Lan: THe Galaxy luminosity function in groups and clusters: the faint-end upturn and the connection to the field luminosity function.
    Mb = -21.30
    ab = -1.03
    Ab = 0.081
    gb = 1.06
    Mf = -18
    af = -1.66
    Af = 0.26
    gf = 0.74
    return schechter_func(mr,m200,Mb,ab,Ab,gb)+schechter_func(mr,m200,Mf,af,Af,gf)


hfile = h5py.File(output,'r')
gal1_r = hfile["gal_core/sdss_r" ][:]
gal2_r = hfile["gal_bhp/sdss_r"  ][:]
gal3_r = hfile["gal_rnd/sdss_r"   ][:]
gal4_r = hfile["gal_blank/sdss_r"][:]
gal5_r = hfile["gal_nfw/sdss_r"][:]

gal1_hm = hfile["gal_core/host_halo_mass"] [:]
gal2_hm = hfile["gal_bhp/host_halo_mass"]  [:]
gal3_hm = hfile["gal_rnd/host_halo_mass"]   [:]
gal4_hm = hfile["gal_blank/host_halo_mass"][:]
gal5_hm = hfile["gal_nfw/host_halo_mass"][:]

gal1_ht = hfile["gal_core/host_halo_tag"] [:]
gal2_ht = hfile["gal_bhp/host_halo_tag"]  [:]
gal3_ht = hfile["gal_rnd/host_halo_tag"]   [:]
gal4_ht = hfile["gal_blank/host_halo_tag"][:]
gal5_ht = hfile["gal_nfw/host_halo_tag"][:]

print "cores: ",gal1_r.size
print "infall prtcl: ",gal2_r.size
print "rnd prtcl: ",gal3_r.size
print "blank: ",gal4_r.size
print "nfw: ",gal5_r.size


indx = np.searchsorted(sod_htag,gal1_ht)
slct1 = sod_htag[indx.clip(min=0,max=sod_htag.size-1)]==gal1_ht
indx1 = indx[slct1]

indx = np.searchsorted(sod_htag,gal2_ht)
slct2 = sod_htag[indx.clip(min=0,max=sod_htag.size-1)]==gal2_ht
indx2 = indx[slct2]

indx = np.searchsorted(sod_htag,gal3_ht)
slct3 = sod_htag[indx.clip(min=0,max=sod_htag.size-1)]==gal3_ht
indx3 = indx[slct3]

indx = np.searchsorted(sod_htag,gal4_ht)
slct4 = sod_htag[indx.clip(min=0,max=sod_htag.size-1)]==gal4_ht
indx4 = indx[slct4]

indx = np.searchsorted(sod_htag,gal5_ht)
slct5 = sod_htag[indx.clip(min=0,max=sod_htag.size-1)]==gal5_ht
indx5 = indx[slct5]

gal1_hm[slct1] = sod_mass[indx1]
gal2_hm[slct2] = sod_mass[indx2]
gal3_hm[slct3] = sod_mass[indx3]
gal4_hm[slct4] = sod_mass[indx4]
gal5_hm[slct5] = sod_mass[indx5]

print gal1_hm.size
print gal2_hm.size
plt.figure(figsize=(10,8))
x_bins = np.linspace(10,16,25)
plt.hist(np.log10(gal1_hm), bins = x_bins, histtype='step',label='core')
#plt.hist(np.log10(gal2_hm), bins = x_bins, histtype='step',label='bhp')
#plt.hist(np.log10(gal3_hm), bins = x_bins, histtype='step',label='hp')
#plt.hist(np.log10(gal4_hm), bins = x_bins, histtype='step',label='glact')
plt.hist(np.log10(gal5_hm), bins = x_bins, histtype='step',label='nfw')
plt.legend()


 
gal1_pi = hfile["gal_core/hostIndex"] [:]
gal2_pi = hfile["gal_bhp/hostIndex"]  [:]
gal3_pi = hfile["gal_rnd/hostIndex"]   [:]
gal4_pi = hfile["gal_blank/hostIndex"][:]
gal5_pi = hfile["gal_nfw/hostIndex"][:]

gal1_ni = hfile["gal_core/nodeIndex"] [:]
gal2_ni = hfile["gal_bhp/nodeIndex"]  [:]
gal3_ni = hfile["gal_rnd/nodeIndex"]   [:]
gal4_ni = hfile["gal_blank/nodeIndex"][:]
gal5_ni = hfile["gal_nfw/nodeIndex"][:]



slct1 = gal1_pi!=gal1_ni
slct2 = gal2_pi!=gal2_ni
slct3 = gal3_pi!=gal3_ni
slct4 = gal4_pi!=gal4_ni
slct5 = gal5_pi!=gal5_ni

print np.sum(slct1),np.sum(slct1==0)
print np.sum(slct2),np.sum(slct2==0)
print np.sum(slct3),np.sum(slct3==0)
print np.sum(slct4),np.sum(slct4==0)
print np.sum(slct5),np.sum(slct5==0)

#gal_r  = np.concatenate((gal1_r,  gal2_r,  gal3_r,  gal4_r))
#gal_hm  = np.concatenate((gal1_hm, gal2_hm, gal3_hm, gal4_hm))

gal_r  = np.concatenate((gal1_r[slct1],  gal2_r[slct2],  gal3_r[slct3],  gal4_r[slct4], gal5_r[slct5]))
gal_hm  = np.concatenate((gal1_hm[slct1], gal2_hm[slct2], gal3_hm[slct3], gal4_hm[slct4], gal5_hm[slct5]))

mass_bins_lin = np.linspace(11,15,10)
mass_bins     = 10**(mass_bins_lin)
mass_bins_avg = 10**((mass_bins_lin[1:]+mass_bins_lin[:-1])/2.0)

print mass_bins_avg

mr_bins = np.linspace(-26,-11,75)
mr_bins_avg = (mr_bins[1:]+mr_bins[:-1])/2.0
mr_bins_width = np.abs(mr_bins[1:]-mr_bins[:-1])

H,x_bins,y_bins = np.histogram2d(gal_r,gal_hm,bins=(mr_bins,mass_bins))
H.astype(dtype='float')

H1,x_bins,y_bins = np.histogram2d(gal1_r[slct1],gal1_hm[slct1],bins=(mr_bins,mass_bins))
H1.astype(dtype='float')

H2,x_bins,y_bins = np.histogram2d(gal2_r[slct2],gal2_hm[slct2],bins=(mr_bins,mass_bins))
H2.astype(dtype='float')

H3,x_bins,y_bins = np.histogram2d(gal3_r[slct3],gal3_hm[slct3],bins=(mr_bins,mass_bins))
H3.astype(dtype='float')

H4,x_bins,y_bins = np.histogram2d(gal4_r[slct4],gal4_hm[slct4],bins=(mr_bins,mass_bins))
H4.astype(dtype='float')

H5,x_bins,y_bins = np.histogram2d(gal5_r[slct5],gal5_hm[slct5],bins=(mr_bins,mass_bins))
H5.astype(dtype='float')


hm_counts,x_bins = np.histogram(halo_cat[step]["sod_halo_mass"],bins=mass_bins)
hm_counts.astype(dtype="float")
print hm_counts
f,axs = plt.subplots(3,3,sharex='col',sharey='all',figsize=(20,15))
for i in range(0,len(mass_bins_avg)):
    ax = axs[i/3,i%3]
    ax.set_title("%.2e<M200c<%.2e[Msun/h]\nN=%d"%(mass_bins[i],mass_bins[i+1],hm_counts[i]))
    m_axis = np.linspace(-11,-26,100)
    ax.plot(m_axis,global_lum_func(m_axis,mass_bins_avg[i]))
    ax.set_yscale('log')
    ax.grid()
    ax.set_ylim([1e-4,1e4])
    ax.set_ylabel("dN/dM [dex^-1]")
    ax.set_xlabel("Mr")

    #ax.scatter(mr_bins_avg,H[:,i]/mr_bins_width/hm_counts,marker='o',facecolor='none')
    if(hm_counts[i]>0):

        y = H[:,i]/mr_bins_width/hm_counts[i]
        err  = np.sqrt(H[:,i])/mr_bins_width/hm_counts[i]
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,1e-4))
        ax.errorbar(mr_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ok',mfc='none',mec='k')

        y = H1[:,i]/mr_bins_width/hm_counts[i]
        err  = np.sqrt(H1[:,i])/mr_bins_width/hm_counts[i]
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,1e-4))
        ax.errorbar(mr_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-or',mfc='none',mec='r')

        y = H2[:,i]/mr_bins_width/hm_counts[i]
        err  = np.sqrt(H2[:,i])/mr_bins_width/hm_counts[i]
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,1e-4))
        ax.errorbar(mr_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ob',mfc='none',mec='b')

        y = H3[:,i]/mr_bins_width/hm_counts[i]*1.1
        err  = np.sqrt(H3[:,i])/mr_bins_width/hm_counts[i]
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,1e-4))
        ax.errorbar(mr_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-sg',mfc='none',mec='g')

        y = H4[:,i]/mr_bins_width/hm_counts[i]
        err  = np.sqrt(H4[:,i])/mr_bins_width/hm_counts[i]
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,1e-4))
        ax.errorbar(mr_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ok',mfc='none',mec='k')
        
        y = H5[:,i]/mr_bins_width/hm_counts[i]
        err  = np.sqrt(H5[:,i])/mr_bins_width/hm_counts[i]
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,1e-4))
        ax.errorbar(mr_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-oc',mfc='none',mec='c')
  
    if(i==len(mass_bins_avg)-1):
        ax.errorbar([],[],yerr=[],fmt='-ok',mfc='none',mec='k',label='all')
        ax.errorbar([],[],yerr=[],fmt='-or',mfc='none',mec='r',label='core pos')
        ax.errorbar([],[],yerr=[],fmt='-ob',mfc='none',mec='b',label='infall prtcl')
        ax.errorbar([],[],yerr=[],fmt='-og',mfc='none',mec='g',label='rnd host prtcl')
        ax.errorbar([],[],yerr=[],fmt='-oc',mfc='none',mec='c',label='nfw profile')
        plt.legend(loc='best',framealpha=0.5)


dtk.save_figs(path='figs/'+param.file+'/gal_lum_func/',extension='.png')

plt.show()
