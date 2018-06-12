#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import h5py 
import sys
import dtk
from util import *


Mstar_func = get_Mstar_func()

def get_conc_from_sod_mass(sod_mass, z):
    mass = sod_mass/Mstar_func(z)
    a = 3.44
    b = 420.49
    c0 = 3.19
    m = -0.10
    return a*((mass/b)**m*(1.0+mass/b)**-m - 1.0)+c0

def get_sod_prop(halo_tag,halo_fof_mass,z):
    indx = np.searchsorted(sod_cat[step]["halo_tag"],halo_tag)
    if(indx < sod_cat[step]["halo_tag"].size):
        if(sod_cat[step]["halo_tag"][indx]==halo_tag):
            r200 = sod_cat[step]["sod_halo_radius"][indx]
            c200 = sod_cat[step]["sod_halo_cdelta"][indx]
            m200 = sod_cat[step]["sod_halo_mass"][indx]
            if(c200 < 0):
                c200 = 3.0
            return m200, r200, c200
    m200 = get_m200(halo_fof_mass,z)
    r200 = get_r200(halo_fof_mass,z)
    c200 = get_conc_from_sod_mass(m200,z)
    return m200, r200, c200
    


param = dtk.Param(sys.argv[1])

mt_loc = param.get_string("mt_loc")
sod_loc = param.get_string("sod_loc")
step = param.get_int('step')
output = param.get_string("output").replace("${step}",str(step))
subfile_num = param.get_int('subfile_num')
NFW_gamma = param.get_float("NFW_gamma")
NFW_maxR200 = param.get_float("NFW_maxR200")
set_gamma(NFW_gamma)
stepz = dtk.StepZ(200,0,500)

mt_cat1 = dtk.Catalog(mt_loc)
mt_cat1.add_step(step)
for i in range(0,subfile_num):
    mt_cat1.add_subfile(i)
mt_cat1.add_var("/forestHalos/timestep",as_name='time_step')
mt_cat1.add_var("/forestHalos/nodeIndex",as_name='nodeIndex')
mt_cat1.add_var("/forestHalos/nodeMass",as_name='halo_mass')
mt_cat1.add_var("/forestHalos/haloTag",as_name='halo_tag')
mt_cat1.add_var("/forestHalos/redshift",as_name='redshift')
mt_cat1.add_var("/forestHalos/timestep",as_name='timestep')
mt_cat1.add_var("/forestHalos/descendentIndex",as_name='descn_indx')
mt_cat1.add_var("/forestHalos/position",as_name='x',index=0);
mt_cat1.add_var("/forestHalos/position",as_name='y',index=1);
mt_cat1.add_var("/forestHalos/position",as_name='z',index=2);
mt_cat1.add_var("/forestHalos/velocity",as_name='vx',index=0);
mt_cat1.add_var("/forestHalos/velocity",as_name='vy',index=1);
mt_cat1.add_var("/forestHalos/velocity",as_name='vz',index=2);
mt_cat1.read_hdf5()
mt_cat1.apply_function('halo_tag',dtk.major_frag_to_real)
slct = mt_cat1[step]['timestep']==step
mt_cat = dtk.Catalog()
mt_cat.select(mt_cat1,slct)

sod_cat = dtk.Catalog(sod_loc)
sod_cat.add_step(step)
sod_cat.add_var("fof_halo_tag",as_name="halo_tag")
sod_cat.add_var("sod_halo_mass")
sod_cat.add_var("sod_halo_radius")
sod_cat.add_var("sod_halo_cdelta")
sod_cat.read_gio()

pos_types = ["gal_core",
             "gal_bhp",
             "gal_rnd",
             "gal_blank",
             "gal_nfw",
             "gal_central"]

def read_var(hfile,name):
    vals = []
    for i,pos_name in enumerate(pos_types):
        vals.append(hfile[pos_name+"/"+name][:])
        print i,vals[-1].size
    val = np.concatenate(vals)
    return val

def read_pos_type(hfile,name):
    vals = []
    for i,pos_name in enumerate(pos_types):
        vals.append(np.ones(hfile[pos_name+"/"+name][:].size,dtype=int)*i)
        print vals[-1][-10:]
    val = np.concatenate(vals)
    return val

def periodic_dist(x1,x2,box_len):
    dx = np.abs(x1-x2)
    dl = box_len/2.0
    slct = dx > dl
    dx[slct] = dx[slct]-box_len
    return dx

srt = np.argsort(mt_cat[step]["halo_tag"])

print output
hfile = h5py.File(output)
gal_mr   = read_var(hfile,"sdss_r")
gal_x   = read_var(hfile,"x")
gal_y   = read_var(hfile,"y")
gal_z   = read_var(hfile,"z")
gal_vx  = read_var(hfile,"vx")
gal_vy  = read_var(hfile,"vy")
gal_vz  = read_var(hfile,"vz")
gal_hmass= read_var(hfile,"host_halo_mass")
gal_htag = read_var(hfile,"host_halo_tag")
gal_id   = read_var(hfile,"nodeIndex")
gal_type = read_pos_type(hfile,"x")
gal_hindx = np.searchsorted(mt_cat[step]["halo_tag"],gal_htag,sorter=srt)
gal_hindx = srt[gal_hindx.clip(max=mt_cat[step]["halo_tag"].size-1)]
 



slct1 = gal_htag == mt_cat[step]["halo_tag"][gal_hindx]
slct2 = gal_mr < -16
slct =  slct1 & slct2
print "select size: ", np.sum(slct),gal_hindx.size
gal_hindx =gal_hindx[slct]
gal_mr = gal_mr[slct]
gal_x  = gal_x[slct]
gal_y  = gal_y[slct]
gal_z  = gal_z[slct]
gal_vx = gal_vx[slct]
gal_vy = gal_vy[slct]
gal_vz = gal_vz[slct]
gal_id = gal_id[slct]
gal_hmass=gal_hmass[slct]
gal_type=gal_type[slct]

gal_hr200 = np.zeros_like(gal_hmass)
print "getting halo r200s.."
for i in range(0,gal_hr200.size):
    if(i%10000 == 0):
        print "i=",i,'/',gal_hr200.size
    m200,r200,c200 = get_sod_prop(gal_htag[i],gal_hmass[i],stepz.get_z(step))
    gal_hr200[i]=r200

slct1 = gal_type==0 #core
slct2 = gal_type==1 #infall
slct3 = gal_type==2 #rnd 
slct4 = gal_type==3 #gltcs
slct5 = gal_type==4 #nfw
slct6 = gal_type==5 #central
print "B:", np.sum(slct1),np.sum(slct2),np.sum(slct3),np.sum(slct4),np.sum(slct5),np.sum(slct6)


gal_hx = mt_cat[step]["x"][gal_hindx]
gal_hy = mt_cat[step]["y"][gal_hindx]
gal_hz = mt_cat[step]["z"][gal_hindx]
gal_dx = periodic_dist(gal_x,gal_hx,256.0)
gal_dy = periodic_dist(gal_y,gal_hy,256.0)
gal_dz = periodic_dist(gal_z,gal_hz,256.0)
gal_dr3 = np.sqrt(gal_dx**2+gal_dy**2+gal_dz**2)/gal_hr200
gal_dr2 = np.sqrt(gal_dx**2+gal_dy**2)/gal_hr200
gal_dr3_mpc = np.sqrt(gal_dx**2+gal_dy**2+gal_dz**2)
gal_dr2_mpc = np.sqrt(gal_dx**2+gal_dy**2)


mass_bins_lin = np.linspace(11,15,10)
mass_bins     = 10**(mass_bins_lin)
mass_bins_avg = 10**((mass_bins_lin[1:]+mass_bins_lin[:-1])/2.0)
mass_bins_width = mass_bins_lin[1:]-mass_bins_lin[:-1]

r_bins = np.linspace(0,4.0,16)
r_bins_avg = (r_bins[1:]+r_bins[:-1])/2.0
r_bins_width = (r_bins[1:]-r_bins[:-1])
r_bins_volume = (r_bins[1:]**3-r_bins[:-1]**3)*4.0/3.0*np.pi


H, x_bins,y_binx = np.histogram2d(gal_hmass,gal_dr3,bins=(mass_bins,r_bins))
H1, x_bins,y_binx = np.histogram2d(gal_hmass[slct1],gal_dr3[slct1],bins=(mass_bins,r_bins))
H2, x_bins,y_binx = np.histogram2d(gal_hmass[slct2],gal_dr3[slct2],bins=(mass_bins,r_bins))
H3, x_bins,y_binx = np.histogram2d(gal_hmass[slct3],gal_dr3[slct3],bins=(mass_bins,r_bins))
H4, x_bins,y_binx = np.histogram2d(gal_hmass[slct4],gal_dr3[slct4],bins=(mass_bins,r_bins))
H5, x_bins,y_binx = np.histogram2d(gal_hmass[slct5],gal_dr3[slct5],bins=(mass_bins,r_bins))
H6, x_bins,y_binx = np.histogram2d(gal_hmass[slct6],gal_dr3[slct6],bins=(mass_bins,r_bins))

print gal_dr3[slct6]
print "starting the plot..."

hm_counts,x_bins = np.histogram(mt_cat[step]['halo_mass'],bins=mass_bins)
hm_counts.astype(dtype="float")
ymin = 1e-6
f,axs = plt.subplots(3,3,sharex='col',sharey='all')
unity = True
for i in range(0,len(mass_bins_avg)):
    ax = axs[i/3,i%3]
    ax.set_title("%.2e<Mfof<%.2e[Msun/h]\nN=%d"%(mass_bins[i],mass_bins[i+1],hm_counts[i]))
    m_axis = np.linspace(-11,-26,100)
    ax.set_yscale('log')
    ax.grid()
    if(unity):
        ax.set_ylabel("gal/volume")
        ax.set_ylim([ymin,1e2])
        ax.set_xlim([0,4.0])
    else:
        ax.set_ylabel("gal/volume [Mpc/h^-3]")
    ax.set_xlabel("r/R200")
    #ax.scatter(mr_bins_avg,H[:,i]/mr_bins_width/hm_counts,marker='o',facecolor='none')
    if(hm_counts[i]>0):
        norm = 1.0
        if(unity):
            norm = 1.0/(np.sum(H[i,:])/hm_counts[i])
            print norm
        y = H[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ok',mfc='none',mec='k')

        if(unity):
            norm = 1.0/np.sum(H1[i,:]/hm_counts[i])
        y = H1[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H1[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-or',mfc='none',mec='r')

        if(unity):
            norm = 1.0/np.sum(H2[i,:]/hm_counts[i])
        y = H2[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H2[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ob',mfc='none',mec='b')

        if(unity):
            norm = 1.0/np.sum(H3[i,:]/hm_counts[i])
        y = H3[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H3[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-og',mfc='none',mec='g')

        if(unity):
            norm = 1.0/np.sum(H4[i,:]/hm_counts[i])
        y = H4[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H4[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ok',mfc='none',mec='k')

        if(unity):
            norm = 1.0/np.sum(H5[i,:]/hm_counts[i])
        y = H5[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H5[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-sc',mfc='none',mec='c')

        if(unity):
            norm = 1.0/np.sum(H6[i,:]/hm_counts[i])
        y = H6[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H6[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-sm',mfc='none',mec='m')

    if(i==len(mass_bins_avg)-1):
        print "ping!!"
        ax.errorbar([],[],yerr=[],fmt='-ok',mfc='none',mec='k',label='all')
        ax.errorbar([],[],yerr=[],fmt='-or',mfc='none',mec='r',label='core pos')
        ax.errorbar([],[],yerr=[],fmt='-ob',mfc='none',mec='b',label='infall prtcl')
        ax.errorbar([],[],yerr=[],fmt='-og',mfc='none',mec='g',label='rnd host prtcl')
        ax.errorbar([],[],yerr=[],fmt='-sc',mfc='none',mec='c',label='NFW profile')
        ax.errorbar([],[],yerr=[],fmt='-sm',mfc='none',mec='m',label='Centrals')
        plt.legend(loc='best',framealpha=0.5)


#########################################
##### Mpc distance, not in r200s ########
#########################################
r_min = 1e-2
r_bins = np.linspace(0,6.0,16)
r_bins = np.logspace(np.log10(r_min),1,16)-r_min
r_bins_avg = (r_bins[1:]+r_bins[:-1])/2.0
r_bins_width = (r_bins[1:]-r_bins[:-1])
r_bins_volume = (r_bins[1:]**3-r_bins[:-1]**3)*4.0/3.0*np.pi


H, x_bins,y_binx = np.histogram2d(gal_hmass,gal_dr3_mpc,bins=(mass_bins,r_bins))
H1, x_bins,y_binx = np.histogram2d(gal_hmass[slct1],gal_dr3_mpc[slct1],bins=(mass_bins,r_bins))
H2, x_bins,y_binx = np.histogram2d(gal_hmass[slct2],gal_dr3_mpc[slct2],bins=(mass_bins,r_bins))
H3, x_bins,y_binx = np.histogram2d(gal_hmass[slct3],gal_dr3_mpc[slct3],bins=(mass_bins,r_bins))
H4, x_bins,y_binx = np.histogram2d(gal_hmass[slct4],gal_dr3_mpc[slct4],bins=(mass_bins,r_bins))
H5, x_bins,y_binx = np.histogram2d(gal_hmass[slct5],gal_dr3_mpc[slct5],bins=(mass_bins,r_bins))
H6, x_bins,y_binx = np.histogram2d(gal_hmass[slct6],gal_dr3_mpc[slct6],bins=(mass_bins,r_bins))

print gal_dr3[slct6]
print "starting the plot..."

hm_counts,x_bins = np.histogram(mt_cat[step]['halo_mass'],bins=mass_bins)
hm_counts.astype(dtype="float")
ymin = 1e-6
f,axs = plt.subplots(3,3,sharex='col')
unity = True
for i in range(0,len(mass_bins_avg)):
    ax = axs[i/3,i%3]
    ax.set_title("%.2e<Mfof<%.2e[Msun/h]\nN=%d"%(mass_bins[i],mass_bins[i+1],hm_counts[i]))
    m_axis = np.linspace(-11,-26,100)
    ax.set_yscale('log')
    ax.grid()
    if(unity):
        ax.set_ylabel("gal/volume")
        ax.set_ylim([ymin,1e6])
        ax.set_xlim([0,4.0])
    else:
        ax.set_ylabel("gal/volume [Mpc/h^-3]")
    ax.set_xlabel("r [Mpc/h]")

    #ax.scatter(mr_bins_avg,H[:,i]/mr_bins_width/hm_counts,marker='o',facecolor='none')
    if(hm_counts[i]>0):
        ax.set_xscale('log')
        norm = 1.0
        if(unity):
            norm = 1.0/(np.sum(H[i,:])/hm_counts[i])
            print norm
        y = H[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ok',mfc='none',mec='k')

        if(unity):
            norm = 1.0/np.sum(H1[i,:]/hm_counts[i])
        y = H1[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H1[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-or',mfc='none',mec='r')

        if(unity):
            norm = 1.0/np.sum(H2[i,:]/hm_counts[i])
        y = H2[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H2[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ob',mfc='none',mec='b')

        if(unity):
            norm = 1.0/np.sum(H3[i,:]/hm_counts[i])
        y = H3[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H3[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-og',mfc='none',mec='g')

        if(unity):
            norm = 1.0/np.sum(H4[i,:]/hm_counts[i])
        y = H4[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H4[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-ok',mfc='none',mec='k')

        if(unity):
            norm = 1.0/np.sum(H5[i,:]/hm_counts[i])
        y = H5[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H5[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-sc',mfc='none',mec='c')

        if(unity):
            norm = 1.0/np.sum(H6[i,:]/hm_counts[i])
        y = H6[i,:]/r_bins_volume/hm_counts[i]*norm
        err  = np.sqrt(H6[i,:])/r_bins_volume/hm_counts[i]*norm
        yerrtop = err 
        yerrbot = y-(np.maximum(y-err,ymin))
        ax.errorbar(r_bins_avg,y,yerr = [yerrbot,yerrtop],fmt='-sm',mfc='none',mec='m')

    if(i==len(mass_bins_avg)-1):
        ax.errorbar([],[],yerr=[],fmt='-ok',mfc='none',mec='k',label='all')
        ax.errorbar([],[],yerr=[],fmt='-or',mfc='none',mec='r',label='core pos')
        ax.errorbar([],[],yerr=[],fmt='-ob',mfc='none',mec='b',label='infall prtcl')
        ax.errorbar([],[],yerr=[],fmt='-og',mfc='none',mec='g',label='rnd host prtcl')
        ax.errorbar([],[],yerr=[],fmt='-sc',mfc='none',mec='c',label='NFW profile')
        ax.errorbar([],[],yerr=[],fmt='-sm',mfc='none',mec='m',label='Centrals')
        ax.legend(loc='best',framealpha=0.5)
dtk.save_figs(path='figs/'+param.file+"/gal_rad_prof/",extension=".png")
plt.show()
