#!/usr/bin/env python2.7
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import sys
import h5py
import dtk


param = dtk.Param(sys.argv[1])

fof_loc = param.get_string("fof_loc")
sod_loc = param.get_string("sod_loc")
bhp_loc = param.get_string("bhp_loc")
mt_loc  = param.get_string("mt_loc")
step   = param.get_int("step")
output = param.get_string("output").replace("${step}",str(step))

fof_cat = dtk.Catalog(fof_loc)
sod_cat = dtk.Catalog(sod_loc)
bhp_cat = dtk.Catalog(bhp_loc)
mt_cat  = dtk.Catalog(mt_loc)


print step

fof_cat.add_step(step)
sod_cat.add_step(step)
bhp_cat.add_step(step)
mt_cat.add_step(step)

fof_cat.add_var("fof_halo_tag")
fof_cat.add_var("fof_halo_mass")
fof_cat.add_var("fof_halo_center_x")
fof_cat.add_var("fof_halo_center_y")
fof_cat.add_var("fof_halo_center_z")

sod_cat.add_var("fof_halo_tag")
sod_cat.add_var("sod_halo_mass")
sod_cat.add_var("sod_halo_radius")

bhp_cat.add_var("fof_halo_tag")
bhp_cat.add_var("x")
bhp_cat.add_var("y")
bhp_cat.add_var("z")

mt_cat.add_var("forestHalos/haloTag",as_name='host_halo_tag')
mt_cat.add_var("forestHalos/timestep",as_name="time_step")

print "loading fof"
fof_cat.read_gio()

print "loading sod"
sod_cat.read_gio()

print "loading mt"
mt_cat.read_hdf5()

print "loading halo prtcls"
bhp_cat.read_gio()


print "merging fof & sod"
halo_cat = dtk.Catalog()
halo_cat.join(fof_cat,sod_cat,remove_matched=False) 

htags = np.unique(bhp_cat[step]["fof_halo_tag"])

arg_srt = np.argsort(halo_cat[step]["fof_halo_mass"])[::-1]

#htags = halo_cat[step]['fof_halo_tag'][arg_srt]

#htags = np.unique(halo_cat[step]['fof_halo_tag'])

hfile = h5py.File(output,'r')
gal_core_x    = hfile["gal_core/x"][:]
gal_core_y    = hfile["gal_core/y"][:]
gal_core_z    = hfile["gal_core/z"][:]
gal_core_htag = hfile["gal_core/host_halo_tag"][:]
gal_core_mr   = hfile["gal_core/sdss_r"][:]
gal_core_id   = hfile["gal_core/nodeIndex"][:]

gal_bhp_x    = hfile["gal_bhp/x"][:]
gal_bhp_y    = hfile["gal_bhp/y"][:]
gal_bhp_z    = hfile["gal_bhp/z"][:]
gal_bhp_htag = hfile["gal_bhp/host_halo_tag"][:]
gal_bhp_mr   = hfile["gal_bhp/sdss_r"][:]
gal_bhp_id   = hfile["gal_bhp/nodeIndex"][:]

gal_rnd_x    = hfile["gal_rnd/x"][:]
gal_rnd_y    = hfile["gal_rnd/y"][:]
gal_rnd_z    = hfile["gal_rnd/z"][:]
gal_rnd_htag = hfile["gal_rnd/host_halo_tag"][:]
gal_rnd_mr   = hfile["gal_rnd/sdss_r"][:]
gal_rnd_id   = hfile["gal_rnd/nodeIndex"][:]

gal_blank_x    = hfile["gal_blank/x"][:]
gal_blank_y    = hfile["gal_blank/y"][:]
gal_blank_z    = hfile["gal_blank/z"][:]
gal_blank_htag = hfile["gal_blank/host_halo_tag"][:]
gal_blank_mr   = hfile["gal_blank/sdss_r"][:]
gal_blank_id   = hfile["gal_blank/nodeIndex"][:]

gal_nfw_x    = hfile["gal_nfw/x"][:]
gal_nfw_y    = hfile["gal_nfw/y"][:]
gal_nfw_z    = hfile["gal_nfw/z"][:]
gal_nfw_htag = hfile["gal_nfw/host_halo_tag"][:]
gal_nfw_mr   = hfile["gal_nfw/sdss_r"][:]
gal_nfw_id   = hfile["gal_nfw/nodeIndex"][:]


gal_cent_x    = hfile['gal_central/x'][:]
gal_cent_y    = hfile['gal_central/y'][:]
gal_cent_z    = hfile['gal_central/z'][:]
gal_cent_htag = hfile["gal_central/host_halo_tag"][:]
gal_cent_mr   = hfile["gal_central/sdss_r"][:]
gal_cent_id   = hfile["gal_central/nodeIndex"][:]


slct = mt_cat[499]["time_step"]==499
print "step 499 halos", np.sum(slct)
htags = np.unique(mt_cat[499]["host_halo_tag"][slct])
print "unique tags: ",htags.size

print "starting loop"
i=0;
for htag in htags:
    i +=1
    #print htag,
    slct = halo_cat[step]["fof_halo_tag"]==htag
    if(np.sum(slct) == 0):
        #print "\tmoving on\n"
        continue
    outlier_limit = 2 #mpc
    outlier_count = 0
    max_outlier = 0
    fof_mass = halo_cat[step]["fof_halo_mass"][slct][0]
    sod_mass = halo_cat[step]["sod_halo_mass"][slct][0]
    fof_x = halo_cat[step]["fof_halo_center_x"][slct][0]
    fof_y = halo_cat[step]["fof_halo_center_y"][slct][0]
    fof_z = halo_cat[step]["fof_halo_center_z"][slct][0]
    sod_rad  = halo_cat[step]["sod_halo_radius"][slct][0]
    bhp_slct = bhp_cat[step]["fof_halo_tag"]==htag
    if(sod_mass < 1e13):
        continue
    #if(np.sum(bhp_slct)==0 ):
    #    continue

    x = bhp_cat[step]["x"][bhp_slct]-fof_x
    y = bhp_cat[step]["y"][bhp_slct]-fof_y
    plt.figure(figsize=(10,8))
    if(np.sum(bhp_slct)!=0):
        H,x_bins,y_bins = np.histogram2d(x,y,bins=(150,150))
        plt.pcolormesh(x_bins,y_bins,H.T+1.0,cmap='gray_r',norm = clrs.LogNorm())
    circle1 = plt.Circle((0,0),sod_rad,facecolor='None',edgecolor='k')
    plt.gca().add_artist(circle1)
    ax = plt.gca()
    ax.add_artist(circle1)
    plt.title(r"z=0    M$_{200}$=%.2e M$_{\odot}$/h"%(sod_mass))
    plt.xlabel('x [comv. Mpc/h]')
    plt.ylabel('y [comv. Mpc/h]')
    plt.grid()
    plt.axis('equal')
    plt.plot([],[],'k-',label=r'R$_{200}$=%.2f Mpc/h'%(sod_rad))



    
    mr_lim = -18
    slct1 = gal_core_htag == htag
    slct2 = gal_core_mr < mr_lim
    slct = slct1 & slct2
    #slct = slct1
    print np.sum(slct)
    if(np.sum(slct)>0):
        size = 10**(gal_core_mr[slct]/-2.5/3.0)/10
        plt.scatter(gal_core_x[slct]-fof_x, gal_core_y[slct]-fof_y,s=size,
                    facecolor='none',edgecolor='red',linewidths=2,label='core [%d]'%np.sum(slct))
        dx = gal_core_x[slct]-fof_x
        dy = gal_core_y[slct]-fof_y
        dz = gal_core_z[slct]-fof_z
        dr = np.sqrt(dx*dx + dy*dy + dz*dz)/sod_rad
        max_outlier = max(max_outlier,np.max(dr))
        outlier_count += np.sum(dr>outlier_limit)
        
    slct1 = gal_bhp_htag == htag
    slct2 = gal_bhp_mr < mr_lim
    slct = slct1 & slct2
    #slct = slct1
    if(np.sum(slct)>0):
        if(np.sum(slct)>0): #there is an issue with zero length s=size array
            size = 10**(gal_bhp_mr[slct]/-2.5/3.0)/10
        else:
            size = None
        plt.scatter(gal_bhp_x[slct]-fof_x, gal_bhp_y[slct]-fof_y,s=size,
                facecolor='none',edgecolor='b',linewidths=1,label='prtcl tag [%d]'%np.sum(slct))
        dx = gal_bhp_x[slct]-fof_x
        dy = gal_bhp_y[slct]-fof_y
        dz = gal_bhp_z[slct]-fof_z
        dr = np.sqrt(dx*dx + dy*dy + dz*dz)/sod_rad
        max_outlier = max(max_outlier,np.max(dr))
        outlier_count += np.sum(dr>outlier_limit)


    #plt.plot(gal_bhp_x[slct]-fof_x, gal_bhp_y[slct]-fof_y,'bx',label='gal bhp')
    slct1 = gal_rnd_htag == htag
    slct2 = gal_rnd_mr < mr_lim
    slct = slct1 & slct2
    print "rnd particles: ", np.sum(slct)
    if(np.sum(slct)>0):
    #slct = slct1
        size = 10**(gal_rnd_mr[slct]/-2.5/3.0)/10
        #if(np.sum(slct) == 0):
        #    size=None
        plt.scatter(gal_rnd_x[slct]-fof_x, gal_rnd_y[slct]-fof_y,s=size,marker='s',
                    facecolor='none',edgecolor='g',linewidths=1,label='rnd halo prtcl [%d]'%np.sum(slct))
        dx = gal_rnd_x[slct]-fof_x
        dy = gal_rnd_y[slct]-fof_y
        dz = gal_rnd_z[slct]-fof_z
        dr = np.sqrt(dx*dx + dy*dy + dz*dz)/sod_rad
        max_outlier = max(max_outlier,np.max(dr))
        outlier_count += np.sum(dr>outlier_limit)

        
    slct1 = gal_nfw_htag == htag
    slct2 = gal_nfw_mr < mr_lim
    slct = slct1 & slct2
    #slct = slct1
    print np.sum(slct)
    size = 10**(gal_nfw_mr[slct]/-2.5/3.0)/10
    if(np.sum(slct) > 0):
        plt.scatter(gal_nfw_x[slct]-fof_x, gal_nfw_y[slct]-fof_y,s=size,marker='s',
                    facecolor='none',edgecolor='c',linewidths=1,label='nfw profile [%d]'%np.sum(slct))


    slct1 = gal_cent_htag == htag
    slct2 = gal_cent_mr < mr_lim
    slct = slct1 & slct2
    print np.sum(slct)
    size = 10**(gal_cent_mr[slct]/-2.5/3.0)/10
    if(np.sum(slct) > 0):
        plt.scatter(gal_cent_x[slct]-fof_x, gal_cent_y[slct]-fof_y,s=size,marker='x',linewidth=3,
                    facecolor='none',edgecolor='m',linewidths=1,label='central [%d]'%np.sum(slct))
    
    
    
    #plt.plot(gal_rnd_x[slct]-fof_x, gal_rnd_y[slct]-fof_y,'gx',label='gal rnd')
    
    #size = 10**(gal_core_mi[slct]/-2.5/3.0)/10
    #plt.scatter(gal_core_x[slct],gal_core_y[slct],c='b')
    #size = 10**(core_sdss_i/-2.5/3.0)/10
    #sp = plt.scatter(core_x,core_y,s=size,c=(core_sdss_u-core_sdss_g),cmap='rainbow',vmin=1,vmax=1.7,lw=0,label='galaxy positions [%d]'%np.sum(core_slct))
    #cb = plt.colorbar(sp)
    #cb.set_label('u-g color')
    plt.legend(loc='best',framealpha=0.5)

    #if(i%100 ==99):
    #    dtk.save_figs(path='figs/'+param.file+'/plot/',extension='.png')
    #plt.show()
    #plt.show()
    x_min,x_max = plt.xlim()
    y_min,y_max = plt.ylim()
    x_min = min(-sod_rad*1.5,x_min)
    y_min = min(-sod_rad*1.5,y_min)
    x_max = max(sod_rad*1.5,x_max)
    print sod_rad*1.5, y_max
    y_max = max(sod_rad*1.5,y_max)
    plt.xlim(x_min,x_max)
    plt.ylim(y_min,y_max)
    if(outlier_count > -1):
       dtk.save_figs(path='figs/'+param.file+'/gal_halo_plot/',extension='.png')
       #plt.show()
    else:
        print "Nothing to show...",max_outlier,"mass: ",sod_mass
    
    plt.close('all')
