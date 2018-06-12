#!/usr/bin/env python2.7
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import sys
import dtk


param = dtk.Param(sys.argv[1])

fof_loc = param.get_string("fof_loc")
sod_loc = param.get_string("sod_loc")
core_loc = param.get_string("core_loc")
gal_loc  = param.get_string("gal_loc")
hp_loc = param.get_string("hp_loc")
mt_loc = param.get_string("mt_loc")
steps   = param.get_int_list("steps")
core_intacted_radius = param.get_float("core_intact_radius")

fof_cat = dtk.Catalog(fof_loc)
sod_cat = dtk.Catalog(sod_loc)
core_cat= dtk.Catalog(core_loc)
hp_cat  = dtk.Catalog(hp_loc)
gal_cat = dtk.Catalog(gal_loc)
mt_cat  = dtk.Catalog(mt_loc)
gal_mt_cat = dtk.Catalog()
step = steps[0]

fof_cat.add_steps(steps)
sod_cat.add_steps(steps)
core_cat.add_steps(steps)
hp_cat.add_steps(steps)
gal_cat.add_steps(steps)
mt_cat.add_steps(steps)

fof_cat.add_var("fof_halo_tag")
fof_cat.add_var("fof_halo_mass")
fof_cat.add_var("fof_halo_center_x")
fof_cat.add_var("fof_halo_center_y")
fof_cat.add_var("fof_halo_center_z")

sod_cat.add_var("fof_halo_tag")
sod_cat.add_var("sod_halo_mass")
sod_cat.add_var("sod_halo_radius")

core_cat.add_var("fof_halo_tag")
core_cat.add_var("x")
core_cat.add_var("y")
core_cat.add_var("z")
core_cat.add_var("radius")
core_cat.add_var("infall_tree_node_index",as_name="nodeIndex")

gal_cat.add_var("/Outputs/Output32/nodeData/nodeIndex",as_name='nodeIndex')
gal_cat.add_var("/Outputs/Output32/nodeData/spheroidLuminositiesStellar:SDSS_u:rest:z0.0000",as_name='sdss_su')
gal_cat.add_var("/Outputs/Output32/nodeData/spheroidLuminositiesStellar:SDSS_g:rest:z0.0000",as_name='sdss_sg')
gal_cat.add_var("/Outputs/Output32/nodeData/spheroidLuminositiesStellar:SDSS_r:rest:z0.0000",as_name='sdss_sr')
gal_cat.add_var("/Outputs/Output32/nodeData/spheroidLuminositiesStellar:SDSS_i:rest:z0.0000",as_name='sdss_si')
gal_cat.add_var("/Outputs/Output32/nodeData/spheroidLuminositiesStellar:SDSS_z:rest:z0.0000",as_name='sdss_sz')
gal_cat.add_var("/Outputs/Output32/nodeData/diskLuminositiesStellar:SDSS_u:rest:z0.0000",as_name='sdss_du')
gal_cat.add_var("/Outputs/Output32/nodeData/diskLuminositiesStellar:SDSS_g:rest:z0.0000",as_name='sdss_dg')
gal_cat.add_var("/Outputs/Output32/nodeData/diskLuminositiesStellar:SDSS_r:rest:z0.0000",as_name='sdss_dr')
gal_cat.add_var("/Outputs/Output32/nodeData/diskLuminositiesStellar:SDSS_i:rest:z0.0000",as_name='sdss_di')
gal_cat.add_var("/Outputs/Output32/nodeData/diskLuminositiesStellar:SDSS_z:rest:z0.0000",as_name='sdss_dz')
gal_cat.add_var("/Outputs/Output32/nodeData/parentIndex",as_name='parentIndex')

hp_cat.add_var("fof_halo_tag")
hp_cat.add_var("x")
hp_cat.add_var("y")
hp_cat.add_var("z")

mt_cat.add_var("/forestHalos/nodeIndex",as_name='parentIndex')
mt_cat.add_var("/forestHalos/nodeMass",as_name='host_halo_mass')


print "loading fof"
fof_cat.read_gio()

print "loading sod"
sod_cat.read_gio()

print "loading core"
core_cat.read_gio()

print "loading galaticus"
gal_cat.read_hdf5()

print "loading mtrees"
mt_cat.read_hdf5(verbose=True);

print "merging galaticus & mtrees"
gal_cat2 = dtk.Catalog()
slct = gal_cat[499]['parentIndex']==-1
slct_sat = gal_cat[499]['parentIndex']!=-1
gal_cat[499]['parentIndex'][slct]=gal_cat[499]['nodeIndex'][slct]
#gal_cat[499]['parentIndex'][slct_sat]=-1
gal_cat2.join(gal_cat,mt_cat,join_on='parentIndex',verbose=True,many_to_one=True)
print "merging cores & galaticus"
gal_core_cat = dtk.Catalog()
gal_core_cat.join(core_cat,gal_cat2,join_on='nodeIndex',verbose=True)

print "loading halo prtcls"
hp_cat.read_gio()


print "galaticus cores w/ mtree for host mass"


print "merging fof & sod"
halo_cat = dtk.Catalog()
halo_cat.join(fof_cat,sod_cat,remove_matched=False) 

gal_core_cat.apply_function('fof_halo_tag',dtk.frag_to_real)

htags = np.unique(hp_cat[step]["fof_halo_tag"])

arg_srt = np.argsort(halo_cat[step]["fof_halo_mass"])[::-1]

htags = halo_cat[step]['fof_halo_tag'][arg_srt]

htags = np.unique(gal_core_cat[step]["fof_halo_tag"])


plt.figure()
plt.title("Galaticus Galaxies w/o Cores")
mass = gal_cat2[499]['host_halo_mass']
sdss_i = -2.5*np.log10(gal_cat2[499]['sdss_di']+gal_cat2[499]['sdss_si'])
H1,x_bins,y_bins = np.histogram2d(mass, sdss_i,bins=(np.logspace(10,15,50),np.linspace(-25,-10,50)))
plt.pcolormesh(x_bins,y_bins,H1.T,cmap='PuBu',norm=clrs.LogNorm())
cb = plt.colorbar()
cb.set_label('density')
plt.xscale('log')
plt.xlabel('Host Halo FoF Mass [Msun/h]')
plt.ylabel('SDSS i-band Magnitude')
plt.tight_layout()
plt.grid()

plt.figure()
plt.title("Galatiucs Galaxies w/ Cores")
mass = gal_core_cat[499]['host_halo_mass']
sdss_i = -2.5*np.log10(gal_core_cat[499]['sdss_di']+gal_core_cat[499]['sdss_si'])
H2,x_bins,y_bins = np.histogram2d(mass, sdss_i,bins=(np.logspace(10,15,50),np.linspace(-25,-10,50)))
plt.pcolormesh(x_bins,y_bins,H2.T,cmap='PuBu',norm=clrs.LogNorm())
cb = plt.colorbar()
cb.set_label('density')
plt.xscale('log')
plt.xlabel('Host Halo FoF Mass [Msun/h]')
plt.ylabel('SDSS i-band Magnitude')
plt.tight_layout()
plt.grid()


plt.figure()
print H1.shape
print H2.shape
H3 = H2.astype(float)/(H1.astype(float)+H2.astype(float))
import numpy.ma as ma
H3T_ma = ma.masked_where(np.isnan(H3.T),H3.T)
#H_fract_w_core = (H1.astype(float)+H2.astype(float))
print H3
plt.title("Fraction of Galaxies w/ Cores")
plt.pcolor(x_bins,y_bins,H3T_ma,cmap='jet',vmax=np.nanmax(H3),vmin=np.nanmin(H3))
cb = plt.colorbar()
cb.set_label('Fraction with Cores')
plt.xscale('log')
plt.xlabel('Host Halo FoF Mass [Msun/h]')
plt.ylabel('SDSS i-band Magnitude')
plt.tight_layout()
plt.grid()



dtk.save_figs(path='figs/'+param.file+'/plot/',extension='.png')
plt.show()



print "starting loop"
i=0;
for htag in htags:
    #h = 709670038
    #h = 1060195455
    #if(htag != h):
    #    continue
    i +=1
    print htag
    slct = halo_cat[step]["fof_halo_tag"]==htag
    if(np.sum(slct) == 0):
        continue
    fof_mass = halo_cat[step]["fof_halo_mass"][slct][0]
    sod_mass = halo_cat[step]["sod_halo_mass"][slct][0]
    fof_x = halo_cat[step]["fof_halo_center_x"][slct][0]
    fof_y = halo_cat[step]["fof_halo_center_y"][slct][0]
    fof_z = halo_cat[step]["fof_halo_center_z"][slct][0]
    sod_rad  = halo_cat[step]["sod_halo_radius"][slct][0]
    hp_slct = hp_cat[step]["fof_halo_tag"]==htag
    core_slct = gal_core_cat[step]['fof_halo_tag']==htag
    print "core num: ", np.sum(core_slct)
    core_x = gal_core_cat[step]['x'][core_slct]-fof_x
    core_y = gal_core_cat[step]['y'][core_slct]-fof_y
    core_z = gal_core_cat[step]['z'][core_slct]-fof_z
    core_r = gal_core_cat[step]['radius'][core_slct]
    core_sdss_u = -2.5*np.log10(gal_core_cat[step]['sdss_du'][core_slct]+gal_core_cat[step]['sdss_su'][core_slct])
    core_sdss_g = -2.5*np.log10(gal_core_cat[step]['sdss_dg'][core_slct]+gal_core_cat[step]['sdss_sg'][core_slct])
    core_sdss_r = -2.5*np.log10(gal_core_cat[step]['sdss_dr'][core_slct]+gal_core_cat[step]['sdss_sr'][core_slct])
    core_sdss_i = -2.5*np.log10(gal_core_cat[step]['sdss_di'][core_slct]+gal_core_cat[step]['sdss_si'][core_slct])
    core_sdss_z = -2.5*np.log10(gal_core_cat[step]['sdss_dz'][core_slct]+gal_core_cat[step]['sdss_sz'][core_slct])
    core_intct = core_r <= core_intacted_radius
    core_dispt = core_r > core_intacted_radius
    print np.nonzero(slct),np.sum(slct)
    print "getting a plot ready"
    print "\t",htag,"%.4e"%(fof_mass),np.sum(hp_slct)

    if(np.sum(hp_slct)==0 ):
        continue

    x = hp_cat[step]["x"][hp_slct]-fof_x
    y = hp_cat[step]["y"][hp_slct]-fof_y
    lim_r = 2
    plt.figure()
    #    H,x_bins,y_bins = np.histogram2d(x,y,bins=(np.linspace(-lim_r,lim_r,150),np.linspace(-lim_r,lim_r,150)))
    H,x_bins,y_bins = np.histogram2d(x,y,bins=(150,150))
    circle1 = plt.Circle((0,0),sod_rad,facecolor='None',edgecolor='k')
    #circle1 = plt.Circle((0,0),.2/0.7,facecolor='None',edgecolor='r',linewidth=3.0)
    plt.gca().add_artist(circle1)
    ax = plt.gca()
    ax.add_artist(circle1)

    # plt.plot(x,y,'.',alpha=0.5)
    plt.title(r"z=0    M$_{200}$=%.2e M$_{\odot}$/h"%(sod_mass))
    plt.pcolormesh(x_bins,y_bins,H.T,cmap='gray_r',norm = clrs.LogNorm())
    plt.xlabel('x [comv. Mpc/h]')
    plt.ylabel('y [comv. Mpc/h]')
    plt.grid()
    #plt.ylim((-lim_r,lim_r))
    #plt.xlim((-lim_r,lim_r))
    plt.axis('equal')
    #plt.plot([],[],'r-',label=r'Radius=200kpc'%(sod_rad))
    plt.plot([],[],'k-',label=r'R$_{200}$=%.2f Mpc/h'%(sod_rad))
    #if(np.sum(core_dispt) > 0):    
    #    plt.plot(core_x[core_dispt]-fof_x,core_y[core_dispt]-fof_y,'s',mec='b',mfc='none',label='Disrupted Cores[%d]'%np.sum(core_dispt))
    #if(np.sum(core_intct) > 0):
    #    plt.plot(core_x[core_intct]-fof_x,core_y[core_intct]-fof_y,'o',mec='r',mfc='none',mew=1.0,label='galaxy positions [%d]'%np.sum(core_intct))
    size = 10**(core_sdss_i/-2.5/3.0)/10
    sp = plt.scatter(core_x,core_y,s=size,c=(core_sdss_u-core_sdss_g),cmap='rainbow',vmin=1,vmax=1.7,lw=0,label='galaxy positions [%d]'%np.sum(core_slct))
    print size
    cb = plt.colorbar(sp)
    cb.set_label('u-g color')
    plt.legend(loc='best',framealpha=0.5)

    #if(i%100 ==99):
    #    dtk.save_figs(path='figs/'+param.file+'/plot/',extension='.png')
    #plt.show()
    #plt.show()
    dtk.save_figs(path='figs/'+param.file+'/plot/',extension='.png')
    #plt.show()
    plt.close('all')
