#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import dtk
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
h=0.702
param = dtk.Param(sys.argv[1])

sod_loc = param.get_string("sod_loc")
core_loc = param.get_string("core_loc")
bhp_loc = param.get_string("bhp_loc")
sbh_loc = param.get_string("sbh_loc")
steps   = param.get_int_list("steps")
step = steps[0]

sod_cat = dtk.Catalog(sod_loc)
core_cat= dtk.Catalog(core_loc)
bhp_cat  = dtk.Catalog(bhp_loc)
sbh_cat  = dtk.Catalog(sbh_loc)




sod_cat.add_steps(steps)
core_cat.add_steps(steps)
sbh_cat.add_steps(steps)
bhp_cat.add_steps(steps)



sod_cat.add_var("fof_halo_tag")
sod_cat.add_var("sod_halo_mass")
sod_cat.add_var("sod_halo_radius")
sod_cat.add_var("fof_halo_center_x")
sod_cat.add_var("fof_halo_center_y")
#sod_cat.add_var("fof_halo_center_z")

core_cat.add_var("fof_halo_tag")
core_cat.add_var("x")
core_cat.add_var("y")
#core_cat.add_var("z")

sbh_cat.add_var("fof_halo_tag")
sbh_cat.add_var("subhalo_mean_x")
sbh_cat.add_var("subhalo_mean_y")
sbh_cat.add_var("subhalo_count")


bhp_cat.add_var("fof_halo_tag")
bhp_cat.add_var("x")
bhp_cat.add_var("y")
#hp_cat.add_var("z")

print "reading sod cat"
sod_cat.read_gio()

print "reading core cat"
core_cat.read_gio()

print "reading subhalo cat"
sbh_cat.read_gio()
sbh_slct = sbh_cat[step]['subhalo_count']>99
for key in sbh_cat[step].keys():
    sbh_cat[step][key]=sbh_cat[step][key][sbh_slct]
print "reading bhp cat"
bhp_cat.read_gio()



sbh_srt = np.argsort(sbh_cat[step]['fof_halo_tag'])
core_srt = np.argsort(core_cat[step]['fof_halo_tag'])
bhp_srt = np.argsort(bhp_cat[step]['fof_halo_tag'])


for i in range(0,sod_cat[499]['fof_halo_tag'].size):
    htag = sod_cat[499]['fof_halo_tag'][i]
    x    =  sod_cat[499]['fof_halo_center_x'][i]/h
    y    =  sod_cat[499]['fof_halo_center_y'][i]/h
    rad  =  sod_cat[499]['sod_halo_radius'][i]/h
    mass =  sod_cat[499]['sod_halo_mass'][i]/h
    if(sod_cat[499]['sod_halo_mass'][i]<5e13):
        continue
    print "\n",i,htag,'%.2e'%sod_cat[499]['sod_halo_mass'][i]
    start,end = dtk.select_sorted(sbh_cat[step]['fof_halo_tag'],htag,sorter=sbh_srt)
    sbh_x = sbh_cat[step]['subhalo_mean_x'][sbh_srt[start:end]]/h-x
    sbh_y = sbh_cat[step]['subhalo_mean_y'][sbh_srt[start:end]]/h-y
    start,end = dtk.select_sorted(core_cat[step]['fof_halo_tag'],htag,sorter=core_srt)
    core_x = core_cat[step]['x'][core_srt[start:end]]/h-x
    core_y = core_cat[step]['y'][core_srt[start:end]]/h-y


    start,end = dtk.select_sorted(bhp_cat[step]['fof_halo_tag'],htag,sorter=bhp_srt)
    bhp_x = bhp_cat[step]['x'][bhp_srt[start:end]]/h-x
    bhp_y = bhp_cat[step]['y'][bhp_srt[start:end]]/h-y
    H,x_bins,y_bins = np.histogram2d(bhp_x,bhp_y,bins=(150,150))
    plt.figure()
    plt.pcolormesh(x_bins,y_bins,H.T,cmap='gray_r',norm = clrs.LogNorm())
    plt.plot(sbh_x,sbh_y,'x',mec='b',mfc='none',label=r'subhalos',mew=1.2,ms=8)
    plt.plot(core_x,core_y,'o',mec='r',mfc='none',label=r'cores',mew=1.0,ms=8)
    circle1 = plt.Circle((0,0),rad,facecolor='None',edgecolor='k')
    plt.gca().add_artist(circle1)
    plt.plot([],[],'k',label=r'R$_{200}$=%.2fMpc'%(rad))
    plt.axis('equal')
    plt.legend(loc='best')
    plt.xlabel(r'x [Mpc]')
    plt.ylabel(r'y [Mpc]')
    ax = plt.gca()
    ax.text(0.1,0.9,r'M$_{200c}$=%.2e M${_\odot}$'%mass,transform=ax.transAxes,fontsize=15)
    plt.grid()
    dtk.save_figs(path='figs/subhalos_cores/',extension='.png')
    plt.close()
    #plt.show()

