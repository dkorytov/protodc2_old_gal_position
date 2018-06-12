#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import h5py 
import sys
import dtk

param = dtk.Param(sys.argv[1])

mt_loc = param.get_string("mt_loc")
fof_loc = param.get_string("fof_loc")
step = param.get_int('step')
output = param.get_string("output").replace("${step}",str(step))
subfile_num = param.get_int('subfile_num')


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
mt_cat1.add_var("/forestHalos/velocity",as_name='vz',index=0);
mt_cat1.add_var("/forestHalos/velocity",as_name='vy',index=1);
mt_cat1.add_var("/forestHalos/velocity",as_name='vx',index=2);
mt_cat1.read_hdf5()
mt_cat1.apply_function('halo_tag',dtk.major_frag_to_real)
slct = mt_cat1[step]['timestep']==step
mt_cat = dtk.Catalog()
mt_cat.select(mt_cat1,slct)

pos_types = ["gal_core",
             "gal_bhp",
             "gal_rnd",
             "gal_blank",
             "gal_nfw",
             "gal_central"]
pos_colors = ['r',
              'b',
              'g',
              'y',
              'c',
              'm']
pos_label = ['core',
             'infall prtcl',
             'rnd prtcl',
             'gltcs',
             'nfw',
             'central']

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

hfile = h5py.File(output)
gal_mr   = read_var(hfile,"sdss_r")
gal_x   = read_var(hfile,"x")
gal_y   = read_var(hfile,"y")
gal_z   = read_var(hfile,"z")
gal_vx  = read_var(hfile,"vz")
gal_vy  = read_var(hfile,"vx")
gal_vz  = read_var(hfile,"vx")
gal_v = np.sqrt(gal_vx**2 + gal_vy**2 + gal_vz**2)
gal_hmass= read_var(hfile,"host_halo_mass")
gal_htag = read_var(hfile,"host_halo_tag")
gal_id   = read_var(hfile,"nodeIndex")
gal_type = read_pos_type(hfile,"x")

srt = np.argsort(mt_cat[step]['halo_tag'])
gal_indx = srt[np.searchsorted(mt_cat[step]['halo_tag'],gal_htag,sorter=srt)]
match = gal_htag == mt_cat[step]['halo_tag'][gal_indx]
print "match: ",np.sum(match),'/',gal_indx.size
gal_hx = mt_cat[step]["x"][gal_indx]
gal_hy = mt_cat[step]["y"][gal_indx]
gal_hz = mt_cat[step]["z"][gal_indx]
gal_hvx = mt_cat[step]['vx'][gal_indx]
gal_hvy = mt_cat[step]['vy'][gal_indx]
gal_hvz = mt_cat[step]['vz'][gal_indx]

vel_bins = np.linspace(-5e3,5e3,200)
vel_bins_avg = (vel_bins[:-1]+vel_bins[1:])/2.0
mass_bins = np.logspace(11,15,10)
mass_bins_avg = (mass_bins[:-1]+mass_bins[1:])/2.0

# slct1 = np.abs(gal_vx-gal_hvx)>700 
# slct2 =gal_type==5
# slct3 = gal_hmass < 1e11
# slct = slct1 & slct2 & slct3
# print np.log10(np.min(gal_hmass[slct])),np.log10(np.max(gal_hmass[slct]))
# print np.min(np.abs(gal_vx[slct]-gal_hvx[slct])),np.max(np.abs(gal_vx[slct]-gal_hvx[slct]))
# for i in range(0,2):
#     print "\n\n\n\ti:",i,np.log10(gal_hmass[slct][i])
#     print "gal type:",pos_types[gal_type[slct][i]]
#     print "gal_hindx:", gal_id[slct][i]
#     print "halo_indx:", mt_cat[step]['nodeIndex'][gal_indx[slct][i]]
#     print "gal_step:", step
#     print "gal host htag:",gal_htag[slct][i]
#     print "halo htag    :",mt_cat[step]['halo_tag'][gal_indx[slct][i]]
#     print "halo_step:",mt_cat[step]['timestep'][gal_indx[slct][i]]
#     print "gal_pos:", gal_x[slct][i],  gal_y[slct][i],  gal_z[slct][i]
#     print "del pos:", gal_hx[slct][i] -gal_x[slct][i],  gal_hy[slct][i] - gal_y[slct][i], gal_hz[slct][i] -  gal_z[slct][i]
#     print "halo_pos:", gal_hx[slct][i],  gal_hy[slct][i],  gal_hz[slct][i]
#     print "gal_vel:", gal_vx[slct][i], gal_vy[slct][i], gal_vz[slct][i]
#     print "halo_vel:", gal_hvx[slct][i], gal_hvy[slct][i], gal_hvz[slct][i]
#     print "del vel:", gal_hvx[slct][i] -gal_vx[slct][i],  gal_hvy[slct][i] - gal_vy[slct][i], gal_hvz[slct][i] -  gal_vz[slct][i]
Hx_all,_,_ = np.histogram2d(gal_vx,gal_hmass,bins=(vel_bins,mass_bins))
Hcx_all,_,_ = np.histogram2d(gal_vx-gal_hvx,gal_hmass,bins=(vel_bins,mass_bins))
Hxs = []
Hcxs = []
for i,pos_type in enumerate(pos_types):
    slct = gal_type == i 
    print i, np.sum(slct)
    Hx,_,_ = np.histogram2d(gal_vx[slct],gal_hmass[slct],bins=(vel_bins,mass_bins))
    Hxs.append(Hx)
    Hcx,_,_= np.histogram2d(gal_vx[slct]-gal_hvx[slct],gal_hmass[slct],bins=(vel_bins,mass_bins))
    Hcxs.append(Hcx)
    
hm_counts,_ = np.histogram(mt_cat[step]['halo_mass'],bins=mass_bins)


f,axs = plt.subplots(3,3,sharex='col')#,sharey='all')
for i in range(0,len(mass_bins_avg)):
    ax = axs[i/3,i%3]
    ax.set_title("%.2e<Mfof<%.2e[Msun/h]\nN=%d"%(mass_bins[i],mass_bins[i+1],hm_counts[i]))

    ax.grid()
    ax.set_ylabel('num density')
    ax.set_xlabel('km/s')
    if(np.sum(Hx_all[:,i])>0):
        ax.plot(vel_bins_avg,Hx_all[:,i].T,'-k')
        ax.set_yscale('log')
        for j in range(0,len(pos_types)):
            if(j==5):
                a = 'x'
            else:
                a = ""
            ax.plot(vel_bins_avg,Hxs[j][:,i],'-'+a+pos_colors[j])
    if(i==0):
        for j in range(0,len(pos_types)):
            ax.plot([],[],'-'+pos_colors[j],label=pos_label[j])
            ax.legend(framealpha=0.5)


f,axs = plt.subplots(3,3,sharex='col')#,sharey='all')
for i in range(0,len(mass_bins_avg)):
    ax = axs[i/3,i%3]
    ax.set_title("%.2e<Mfof<%.2e[Msun/h]\nN=%d"%(mass_bins[i],mass_bins[i+1],hm_counts[i]))

    ax.grid()
    ax.set_ylabel('num density')
    ax.set_xlabel('km/s')
    if(np.sum(Hx_all[:,i])>0):
        ax.plot(vel_bins_avg,Hcx_all[:,i].T,'-k')
        ax.set_yscale('log')
        for j in range(0,len(pos_types)):
            if(j==5 or j==0):
                a = 'x'
            else:
                a = ""
            ax.plot(vel_bins_avg,Hcxs[j][:,i],'-'+a+pos_colors[j])
    if(i==0):
        for j in range(0,len(pos_types)):
            ax.plot([],[],'-'+pos_colors[j],label=pos_label[j])
            ax.legend(framealpha=0.5)

# plt.figure()
# plt.title('absolute velocity')
# plt.plot(vel_bins_avg,Hx_all,'-k',label='all')
# for i in range(0,len(Hxs)):
#     plt.plot(vel_bins_avg,Hxs[i],'-'+pos_colors[i],label=pos_label[i])
# plt.legend()
# plt.yscale('log')

# plt.figure()
# plt.title('relative velocity')
# plt.plot(vel_bins_avg,Hcx_all,'-k',label='all')
# for i in range(0,len(Hxs)):
#     plt.plot(vel_bins_avg,Hcxs[i],'-'+pos_colors[i],label=pos_label[i])
# plt.legend()
# plt.yscale('log')
dtk.save_figs(path='figs/'+param.file+"/gal_vel_prof/",extension=".png")
plt.show()
