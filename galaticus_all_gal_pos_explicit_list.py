#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import h5py 
import sys

import dtk
from util import *
def rand_unit_vect():
    vecs = np.random.normal(size )
    mags = np.linalg.norm(vecs,axis=-1)

def lum_to_mag(data1, data2):
    return -2.5*(np.log10(data1+data2))

def find_host(hostIndex,host_halo_tag,host_time_step,
              nodeIndex,node_halo_tag,node_time_step,
              nodeDescn,target_step):
    print "starting to find the host halo v1"
    srt = np.argsort(nodeIndex)
    for i in range(0,len(hostIndex)):
        if(i%10000 ==0):
            print float(i)/float(len(hostIndex))
        indx = i
        while(node_time_step[indx] < target_step):
            old_time = node_time_step[indx]
            descn = nodeDescn[indx]
            indx = srt[np.searchsorted(nodeIndex,descn,sorter=srt)]
            new_time = node_time_step[indx]
            if(new_time-old_time <= 0):
                print "**********"
            #print "step: ", node_time_step[indx], indx
        hostIndex[i] = indx
        host_halo_tag[i] = node_halo_tag[indx]
        host_time_step[i] = node_time_step[indx]

def find_host2(hostIndex,host_halo_tag,host_time_step,host_halo_mass,
               host_halo_x,host_halo_y, host_halo_z, host_halo_vx,host_halo_vy, host_halo_vz, 
               nodeIndex,node_halo_tag,node_time_step,node_halo_mass,
               node_halo_x,node_halo_y, node_halo_z, node_halo_vx,node_halo_vy, node_halo_vz, 
               nodeDescn,target_step):
    srt = np.argsort(nodeIndex)
    slct = host_time_step < target_step
    i=0;
    print "mt size: ",nodeIndex.size,node_halo_tag.size,node_time_step.size,nodeDescn.size
    print "gal size: ",hostIndex.size,host_halo_tag.size,host_time_step.size
    while(np.sum(slct)>0):
        i= i+1
        print "\n***",i,np.sum(slct)
        indx = srt[np.searchsorted(nodeIndex,hostIndex[slct],sorter=srt)]
        descn = nodeDescn[indx]
        indx_d = srt[np.searchsorted(nodeIndex,descn,sorter=srt)]
        hostIndex[slct]=nodeIndex[indx_d]
        host_time_step[slct]=node_time_step[indx_d]
        host_halo_tag[ slct] =node_halo_tag[indx_d]
        host_halo_mass[slct]=node_halo_mass[indx_d]
        host_halo_x[slct] = node_halo_x[indx_d]
        host_halo_y[slct] = node_halo_y[indx_d]
        host_halo_z[slct] = node_halo_z[indx_d]
        host_halo_vx[slct] = node_halo_vx[indx_d]
        host_halo_vy[slct] = node_halo_vy[indx_d]
        host_halo_vz[slct] = node_halo_vz[indx_d]
        slct = host_time_step < target_step

Mstar_func = get_Mstar_func()

def get_conc_from_sod_mass(sod_mass, z):
    mass = sod_mass/Mstar_func(z)
    a = 3.44
    b = 420.49
    c0 = 3.19
    m = -0.10
    return a*((mass/b)**m*(1.0+mass/b)**-m - 1.0)+c0

def get_sod_prop(halo_tag,halo_fof_mass,z):
    indx = np.searchsorted(sod_cat[step]["fof_halo_tag"],halo_tag)
    if(indx < sod_cat[step]["fof_halo_tag"].size):
        if(sod_cat[step]["fof_halo_tag"][indx]==halo_tag):
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
    
fof_sort = None
def get_fof_prop(halo_tag):
    global fof_sort
    if(fof_sort == None):
        fof_sort = np.argsort(fof_cat[step]["fof_halo_tag"])
    indx = dtk.search_sorted(fof_cat[step]["fof_halo_tag"],halo_tag,sorter=fof_sort)
    if(indx == -1):
        print "fof halo not found?!",halo_tag
        x = -1;
        y = -1;
        z = -1;
        vx = -1;
        vy = -1;
        vz = -1;
    else:
        x  = fof_cat[step]['x' ][indx]
        y  = fof_cat[step]['y' ][indx]
        z  = fof_cat[step]['z' ][indx]
        vx = fof_cat[step]['vx'][indx]
        vy = fof_cat[step]['vy'][indx]
        vz = fof_cat[step]['vz'][indx]
    return x,y,z,vx,vy,vz



stepz = dtk.StepZ(200,0,500)

param = dtk.Param(sys.argv[1])
fof_loc = param.get_string("fof_loc")
sod_loc = param.get_string("sod_loc")
core_loc = param.get_string("core_loc")
gal_loc  = param.get_string("gal_loc")
bhp_loc = param.get_string("bhp_loc")
mt_loc = param.get_string("mt_loc")
hp_loc = param.get_string("hp_loc")
step = param.get_int("step")
galaticus_step = param.get_int("galaticus_step")
galaticus_z    = param.get_string("galaticus_z")
subfile_num = param.get_int("subfile_num")
subfile_list = param.get_string_list('subfile_list')
output = param.get_string("output").replace("${step}",str(step))
rnd_halo_prtcl = param.get_bool("rnd_halo_prtcl")
NFW_profile = param.get_bool("NFW_profile")
NFW_gamma = param.get_float("NFW_gamma")
NFW_maxR200 = param.get_float("NFW_maxR200")
fof_cat = dtk.Catalog(fof_loc)
fof_cat.add_step(step)

sod_cat = dtk.Catalog(sod_loc)
sod_cat.add_step(step)

gal_cat = dtk.Catalog(gal_loc)
gal_cat.add_step(step,in_file_step=galaticus_step)

core_cat = dtk.Catalog(core_loc)
core_cat.add_step(step)

mt_cat = dtk.Catalog(mt_loc)
mt_cat.add_step(step)

for i in range(0,subfile_num):
    mt_cat.add_subfile(i)

gal_cat.set_explicit_files(subfile_list)


core_cat.add_var("fof_halo_tag")
core_cat.add_var("x",as_name='core_x')
core_cat.add_var("y",as_name='core_y')
core_cat.add_var("z",as_name='core_z')
core_cat.add_var("vx",as_name='core_vx')
core_cat.add_var("vy",as_name='core_vy')
core_cat.add_var("vz",as_name='core_vz')
core_cat.add_var("radius")
core_cat.add_var("infall_tree_node_index",as_name="nodeIndex")
core_cat.read_gio()

gal_cat.add_var("/Outputs/Output${step}/nodeData/nodeIndex",as_name='nodeIndex')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_u:rest:z"+galaticus_z,as_name='sdss_u')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_g:rest:z"+galaticus_z,as_name='sdss_g')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_r:rest:z"+galaticus_z,as_name='sdss_r')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_i:rest:z"+galaticus_z,as_name='sdss_i')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_z:rest:z"+galaticus_z,as_name='sdss_z')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidMassStellar",as_name='stellar_mass')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_u:rest:z"+galaticus_z,as_name='sdss_du')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_g:rest:z"+galaticus_z,as_name='sdss_dg')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_r:rest:z"+galaticus_z,as_name='sdss_dr')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_i:rest:z"+galaticus_z,as_name='sdss_di')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_z:rest:z"+galaticus_z,as_name='sdss_dz')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskMassStellar",as_name='stellar_mass_disk')

gal_cat.add_var("/Outputs/Output${step}/nodeData/positionPositionX",as_name="gltcs_x")
gal_cat.add_var("/Outputs/Output${step}/nodeData/positionPositionY",as_name="gltcs_y")
gal_cat.add_var("/Outputs/Output${step}/nodeData/positionPositionZ",as_name="gltcs_z")
gal_cat.add_var("/Outputs/Output${step}/nodeData/positionVelocityX",as_name="gltcs_vx")
gal_cat.add_var("/Outputs/Output${step}/nodeData/positionVelocityY",as_name="gltcs_vy")
gal_cat.add_var("/Outputs/Output${step}/nodeData/positionVelocityZ",as_name="gltcs_vz")
gal_cat.add_var("/Outputs/Output${step}/nodeData/nodeIndex",as_name='hostIndex')
gal_cat.add_var("/Outputs/Output${step}/nodeData/nodeIsIsolated",as_name="isIsolated")
print "reading gal catalog "
gal_cat.read_hdf5()
print "done"

gal_cat[step]['sdss_u'] = lum_to_mag(gal_cat[step]['sdss_u'],gal_cat[step]['sdss_du'])
gal_cat[step]['sdss_g'] = lum_to_mag(gal_cat[step]['sdss_g'],gal_cat[step]['sdss_dg'])
gal_cat[step]['sdss_r'] = lum_to_mag(gal_cat[step]['sdss_r'],gal_cat[step]['sdss_dr'])
gal_cat[step]['sdss_i'] = lum_to_mag(gal_cat[step]['sdss_i'],gal_cat[step]['sdss_di'])
gal_cat[step]['sdss_z'] = lum_to_mag(gal_cat[step]['sdss_z'],gal_cat[step]['sdss_dz'])
gal_cat[step]['stellar_mass'] = gal_cat[step]['stellar_mass'] + gal_cat[step]['stellar_mass_disk']


slct = gal_cat[step]["hostIndex"]==-1
gal_cat[step]["hostIndex"][slct]=gal_cat[step]["nodeIndex"][slct]

mt_cat.add_var("/forestHalos/nodeIndex",as_name='nodeIndex')
mt_cat.add_var("/forestHalos/nodeMass",as_name='infall_halo_mass')
mt_cat.add_var("/forestHalos/haloTag",as_name='infall_halo_tag')
mt_cat.add_var("/forestHalos/redshift",as_name='infall_redshift')
mt_cat.add_var("/forestHalos/timestep",as_name='infall_timestep')
mt_cat.add_var("/forestHalos/descendentIndex",as_name='descn_indx')
mt_cat.add_var("/forestHalos/position",as_name='infall_halo_x',index=0);
mt_cat.add_var("/forestHalos/position",as_name='infall_halo_y',index=1);
mt_cat.add_var("/forestHalos/position",as_name='infall_halo_z',index=2);
mt_cat.add_var("/forestHalos/velocity",as_name='infall_halo_vx',index=0);
mt_cat.add_var("/forestHalos/velocity",as_name='infall_halo_vy',index=1);
mt_cat.add_var("/forestHalos/velocity",as_name='infall_halo_vz',index=2);

mt_cat.read_hdf5()

mt_cat.apply_function('infall_halo_tag',dtk.major_frag_to_real)

sod_cat.add_var("fof_halo_tag")
sod_cat.add_var("sod_halo_mass")
sod_cat.add_var("sod_halo_radius")
sod_cat.add_var("sod_halo_cdelta")
sod_cat.add_var("fof_halo_center_x", as_name ='x')
sod_cat.add_var("fof_halo_center_y", as_name ='y')
sod_cat.add_var("fof_halo_center_z", as_name ='z')
sod_cat.read_gio()

fof_cat.add_var("fof_halo_tag")
fof_cat.add_var("fof_halo_mass")
fof_cat.add_var("fof_halo_center_x",as_name ='x')
fof_cat.add_var("fof_halo_center_y",as_name ='y')
fof_cat.add_var("fof_halo_center_z",as_name ='z')
fof_cat.add_var("fof_halo_mean_vx",as_name ='vx')
fof_cat.add_var("fof_halo_mean_vy",as_name ='vy')
fof_cat.add_var("fof_halo_mean_vz",as_name ='vz')
fof_cat.read_gio()


gal_cat2 = dtk.Catalog()
print "\n\nJoining galaticus & merger trees for nodeIndex"
#gal_cat2.join(gal_cat,mt_cat,join_on='nodeIndex',verbose=True,random_to_one=True,remove_matched=False)
gal_cat2.quick_join(gal_cat,mt_cat,'nodeIndex')

gal_cat3 = dtk.Catalog()
print "\n\nJoining galaticus & merger trees for hostIndex"
mt_cat.rename("nodeIndex",       "hostIndex")
mt_cat.rename("infall_halo_tag", "host_halo_tag")
mt_cat.rename("infall_halo_mass","host_halo_mass")
mt_cat.rename("infall_redshift", "host_redshift")
mt_cat.rename("infall_timestep", "host_timestep")
mt_cat.rename("infall_halo_x","host_halo_x")
mt_cat.rename("infall_halo_y","host_halo_y")
mt_cat.rename("infall_halo_z","host_halo_z")
mt_cat.rename("infall_halo_vx","host_halo_vx")
mt_cat.rename("infall_halo_vy","host_halo_vy")
mt_cat.rename("infall_halo_vz","host_halo_vz")
#gal_cat3.join(gal_cat2,mt_cat,join_on='hostIndex',verbose=True,many_to_one=True,remove_matched=False)
gal_cat3.quick_join(gal_cat2,mt_cat,'hostIndex')
find_host2(gal_cat3[step]["hostIndex"],gal_cat3[step]["host_halo_tag"],gal_cat3[step]["host_timestep"],gal_cat3[step]["host_halo_mass"],
           gal_cat3[step]["host_halo_x"], gal_cat3[step]["host_halo_y"], gal_cat3[step]["host_halo_z"],gal_cat3[step]["host_halo_vx"], gal_cat3[step]["host_halo_vy"], gal_cat3[step]["host_halo_vz"],
           mt_cat[step] ["hostIndex"],  mt_cat[step]["host_halo_tag"],  mt_cat[step]["host_timestep"],mt_cat[step]["host_halo_mass"],
           mt_cat[step]["host_halo_x"],mt_cat[step]["host_halo_y"],mt_cat[step]["host_halo_z"],mt_cat[step]["host_halo_vx"],mt_cat[step]["host_halo_vy"],mt_cat[step]["host_halo_vz"],
           mt_cat[step]["descn_indx"],step)

# ni = 2247296214057881197
# slct = mt_cat[step]['hostIndex']==ni
# mt_htag = mt_cat[step]['host_halo_tag'][slct]
# fof_slct = fof_cat[step]['fof_halo_tag']==mt_htag
# print "mtree fof tag via nodeIndex: ",mt_cat[step]['host_halo_tag'][slct]
# print "mtree pos: ",mt_cat[step]['host_halo_x'][slct],mt_cat[step]['host_halo_y'][slct],mt_cat[step]['host_halo_z'][slct]
# print "mtree vel: ",mt_cat[step]['host_halo_vx'][slct],mt_cat[step]['host_halo_vy'][slct],mt_cat[step]['host_halo_vz'][slct]
# print "fof cat pos: ",fof_cat[step]['x'][fof_slct],fof_cat[step]['y'][fof_slct],fof_cat[step]['z'][fof_slct]
# print "fof cat vel: ",fof_cat[step]['vx'][fof_slct],fof_cat[step]['vy'][fof_slct],fof_cat[step]['vz'][fof_slct]
# print "======================================="

# #copy the fof host halo position & vel into galaticus pos & velocity
# h=0.71*stepz.get_a(step)*1.07*1.0054
# diff_x = (0.1<np.abs(gal_cat3[step]['gltcs_x'] - gal_cat3[step]['host_halo_x']/h)) & (gal_cat3[step]['isIsolated']==1)
# diff_y = 0.1<(gal_cat3[step]['gltcs_y'] - gal_cat3[step]['host_halo_y']/h)
# diff_z = 0.1<(gal_cat3[step]['gltcs_z'] - gal_cat3[step]['host_halo_z']/h)
# diff_vx = gal_cat3[step]['gltcs_vx'] != gal_cat3[step]['host_halo_vx']
# diff_vy = gal_cat3[step]['gltcs_vy'] != gal_cat3[step]['host_halo_vy']
# diff_vz = gal_cat3[step]['gltcs_vz'] != gal_cat3[step]['host_halo_vz']


# print "this is number of different Xs: ",np.sum(diff_x)
# print gal_cat3[step]['gltcs_x'][diff_x], gal_cat3[step]['host_halo_x'][diff_x]/h, np.average(gal_cat3[step]['gltcs_x'][diff_x]/ gal_cat3[step]['host_halo_x'][diff_x]*h)
# print "this is number of different Ys: ",np.sum(diff_y)
# print gal_cat3[step]['gltcs_y'][diff_y], gal_cat3[step]['host_halo_y'][diff_y]/h, np.average(gal_cat3[step]['gltcs_y'][diff_y]/ gal_cat3[step]['host_halo_y'][diff_y]*h)
# print "this is number of different Zs: ",np.sum(diff_z)
# print gal_cat3[step]['gltcs_z'][diff_z], gal_cat3[step]['host_halo_z'][diff_z]/h, np.average(gal_cat3[step]['gltcs_z'][diff_z]/ gal_cat3[step]['host_halo_z'][diff_z]*h)
# print "this is number of different VXs: ",np.sum(diff_vx)
# print np.average(gal_cat3[step]['gltcs_vx'] /gal_cat3[step]['host_halo_vx'])
# print "this is number of different VYs: ",np.sum(diff_vy)
# print np.average(gal_cat3[step]['gltcs_vy'] /gal_cat3[step]['host_halo_vy'])
# print "this is number of different VZs: ",np.sum(diff_vz)
# print np.average(gal_cat3[step]['gltcs_vz'] /gal_cat3[step]['host_halo_vz'])

# #gal_indx = 32150
# #fof_indx = 3693
# #mt_indx  = 514211
# gal_indx = np.where(diff_x)[0][0]
# htag = gal_cat3[step]['host_halo_tag'][gal_indx]
# fof_indx = np.where(fof_cat[step]['fof_halo_tag']==htag)[0]
# mt_indx  = np.where( (mt_cat[step]['host_halo_tag']==htag) & (mt_cat[step]['host_timestep']==step))[0]

# print "gal_cat index=",np.where(diff_x)
# htag = gal_cat3[step]['host_halo_tag'][gal_indx]
# print "\n halo tag: ",htag
# print "isIsolate: ",gal_cat3[step]['isIsolated'][gal_indx]
# slct = fof_cat[step]['fof_halo_tag']==htag
# print "fof indx:",np.where(slct)
# print "fof mass: %.1e"%fof_cat[step]['fof_halo_mass'][fof_indx]
# print "fof pos: ",fof_cat[step]['x'][fof_indx],fof_cat[step]['y'][fof_indx],fof_cat[step]['z'][fof_indx]
# slct = (mt_cat[step]['host_halo_tag']==htag) & (mt_cat[step]['host_timestep']==step)
# print "mt_cat indx:",np.where(slct)
# print "mtree htag: ",mt_cat[step]['host_halo_tag'][mt_indx]
# print "mtree pos: ",mt_cat[step]['host_halo_x'][mt_indx],mt_cat[step]['host_halo_y'][mt_indx],mt_cat[step]['host_halo_z'][mt_indx]
# print "gltcs pos: ",gal_cat3[step]['gltcs_x'][gal_indx]*h,gal_cat3[step]['gltcs_y'][gal_indx]*h,gal_cat3[step]['gltcs_z'][gal_indx]*h
# print "host pos:  ",gal_cat3[step]['host_halo_x'][gal_indx],gal_cat3[step]['host_halo_y'][gal_indx],gal_cat3[step]['host_halo_z'][gal_indx]
# print ""
# print "gltcs vel: ",gal_cat3[step]['gltcs_vx'][gal_indx],gal_cat3[step]['gltcs_vy'][gal_indx],gal_cat3[step]['gltcs_vz'][gal_indx]
# print "mtree vel: ",mt_cat[step]['host_halo_vx'][mt_indx],mt_cat[step]['host_halo_vy'][mt_indx],mt_cat[step]['host_halo_vz'][mt_indx]
# print "fof vel:   ",fof_cat[step]['vx'][fof_indx],fof_cat[step]['vy'][fof_indx],fof_cat[step]['vz'][fof_indx]
# print "hostIndex: ",gal_cat3[step]['hostIndex'][gal_indx]
# diff_slct = gal_cat3[step]['isIsolated'] == 1
# slct = gal_cat3[step]['hostIndex'][diff_slct]==gal_cat3[step]['nodeIndex'][diff_slct]
# print "host/infall indx for centrals: ",np.sum(slct)/slct.size,np.sum(slct),slct.size

# diff = gal_cat3[step]['infall_halo_x'][diff_slct]==gal_cat3[step]['host_halo_x'][diff_slct]
# print "diff infall/host central x pos: ",np.sum(diff)/diff.size, np.sum(diff), diff.size
# diff = gal_cat3[step]['infall_halo_y'][diff_slct]==gal_cat3[step]['host_halo_y'][diff_slct]
# print "diff infall/host central y pos: ",np.sum(diff)/diff.size, np.sum(diff), diff.size
# diff = gal_cat3[step]['infall_halo_z'][diff_slct]==gal_cat3[step]['host_halo_z'][diff_slct]
# print "diff infall/host central z pos: ",np.sum(diff)/diff.size, np.sum(diff), diff.size

# diff = gal_cat3[step]['infall_halo_x'][diff_slct==0]==gal_cat3[step]['host_halo_x'][diff_slct==0]
# print "diff infall/host sat x pos: ",np.sum(diff)/diff.size, np.sum(diff), diff.size

# slct = gal_cat3[step]['host_timestep'][diff_slct]==step
# print "same host step central: ",np.sum(slct)/slct.size, np.sum(slct),slct.size
# slct = gal_cat3[step]['infall_timestep'][diff_slct]==step
# print "same infall step central: ",np.sum(slct)/slct.size, np.sum(slct),slct.size

# slct = gal_cat3[step]['host_timestep'][diff_slct==0]==step
# print "same host step sat: ",np.sum(slct)/slct.size, np.sum(slct),slct.size
# slct = gal_cat3[step]['infall_timestep'][diff_slct==0]==step
# print "same infall step sat: ",np.sum(slct)/slct.size, np.sum(slct),slct.size


# diff_x_val = np.abs(gal_cat3[step]['gltcs_x'][diff_slct]*h - gal_cat3[step]['infall_halo_x'][diff_slct])
# diff_y_val = np.abs(gal_cat3[step]['gltcs_y'][diff_slct]*h - gal_cat3[step]['infall_halo_y'][diff_slct])
# diff_z_val = np.abs(gal_cat3[step]['gltcs_z'][diff_slct]*h - gal_cat3[step]['infall_halo_z'][diff_slct])
# diff_vx_val = np.abs(gal_cat3[step]['gltcs_vx'][diff_slct]*h - gal_cat3[step]['infall_halo_vx'][diff_slct])
# diff_vy_val = np.abs(gal_cat3[step]['gltcs_vy'][diff_slct]*h - gal_cat3[step]['infall_halo_vy'][diff_slct])
# diff_vz_val = np.abs(gal_cat3[step]['gltcs_vz'][diff_slct]*h - gal_cat3[step]['infall_halo_vz'][diff_slct])

# diff_r_val = np.sqrt(diff_x_val**2 + diff_y_val**2 + diff_z_val**2)

# wrong_slct = diff_r_val > 5e-3
# print "====================================="
# print "fract wrong: ", np.sum(wrong_slct)/wrong_slct.size, np.sum(wrong_slct), wrong_slct.size
# print "fract again: ", float(np.sum(wrong_slct))/float(wrong_slct.size)
# print "====================================="
# plt.figure()
# plt.plot(gal_cat3[step]['host_halo_mass'][diff_slct],diff_x_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel('diff x')
# plt.xlabel('Mfof')
# plt.grid()

# plt.figure()
# plt.title("Centrals Only")
# plt.plot(gal_cat3[step]['host_halo_z'][diff_slct],diff_z_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.ylabel('diff z [Mpc/h]')
# plt.xlabel('fof z position [Mpc/h]')
# plt.grid()

# plt.figure()
# plt.title("Centrals Only")
# plt.plot(gal_cat3[step]['host_halo_mass'][diff_slct],diff_r_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel('diff r [Mpc/h]')
# plt.xlabel('Mfof Host[Msun/h]')
# plt.grid()


# plt.figure()
# plt.plot(diff_x_val,diff_y_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel('diff y')
# plt.xlabel('diff x')
# plt.grid()

# plt.figure()
# plt.plot(diff_x_val,diff_z_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel('diff z')
# plt.xlabel('diff x')
# plt.grid()

# plt.figure()
# plt.plot(diff_y_val,diff_z_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel('diff z')
# plt.xlabel('diff y')
# plt.grid()

# plt.figure()
# plt.plot(gal_cat3[step]['host_halo_mass'][diff_slct],diff_vx_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel('delta vx [km/s]')
# plt.xlabel('Mfof Host [Msun/h]')

# plt.figure()
# plt.plot(diff_x_val,diff_vx_val,'.',alpha=0.2)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel('delta vx [km/s]')
# plt.xlabel('Mfof Host [Msun/h]')


# plt.show()

# exit()

##################################################
# This is the end of the testing phase of things #
##################################################


gal_cat3[step]['gltcs_x']=np.copy(gal_cat3[step]['host_halo_x'])
gal_cat3[step]['gltcs_y']=np.copy(gal_cat3[step]['host_halo_y'])
gal_cat3[step]['gltcs_z']=np.copy(gal_cat3[step]['host_halo_z'])
gal_cat3[step]['gltcs_vx']=np.copy(gal_cat3[step]['host_halo_vx'])
gal_cat3[step]['gltcs_vy']=np.copy(gal_cat3[step]['host_halo_vy'])
gal_cat3[step]['gltcs_vz']=np.copy(gal_cat3[step]['host_halo_vz'])

gal_central = dtk.Catalog()
print "all count = ",gal_cat3[step]['host_halo_tag'].size
slct = gal_cat3[step]["isIsolated"]==1
print "central count =",np.sum(slct)
print "non-central count = ",np.sum(slct==0)
gal_central.select(gal_cat3,slct)
print "gal_central count = ",gal_central[step]['host_halo_tag'].size
print "gal_cat3 count = ",gal_cat3[step]['host_halo_tag'].size
#this catalog will have all the galaxies that have cores
gal_cores = dtk.Catalog()
print "\n\nJoining galaticus & core catalog on nodeIndex"
#gal_cores.join(gal_cat3,core_cat,join_on='nodeIndex')
gal_cores.quick_join(gal_cat3,core_cat,'nodeIndex',one_to_random=True)

print "these are the keys here..",gal_cat3.step_data.keys()

bhp_cat = dtk.Catalog(bhp_loc)
bhp_cat.add_step(step)
bhp_cat.add_var("id",as_name="infall_halo_tag")
bhp_cat.add_var("x",as_name="bhp_x")
bhp_cat.add_var("y",as_name="bhp_y")
bhp_cat.add_var("z",as_name="bhp_z")
bhp_cat.add_var("vx",as_name="bhp_vx")
bhp_cat.add_var("vy",as_name="bhp_vy")
bhp_cat.add_var("vz",as_name="bhp_vz")
bhp_cat.add_var("fof_halo_tag",as_name="host_halo_tag")
print "reading big halo particle file..."
bhp_cat.read_gio(verbose=True)
#bhp_cat.read_none()

gal_bhp = dtk.Catalog()
print "Joining non-core galatiucs & big halo particles"
gal_bhp.quick_join(gal_cat3,bhp_cat,'infall_halo_tag',req_also='host_halo_tag')
#gal_bhp.join(gal_cat3,bhp_cat,join_on='infall_halo_tag',req_also='host_halo_tag') #just to use with read none


gal_rnd = dtk.Catalog()
if(rnd_halo_prtcl):
    gal_rnd.quick_join(gal_cat3,bhp_cat,'host_halo_tag',one_to_random=True)
    #gal_rnd.join(gal_cat3,bhp_cat,join_on='host_halo_tag') #just to use with read none

if(NFW_profile):
    set_gamma(NFW_gamma)
    srt = np.argsort(sod_cat[step]["fof_halo_tag"])
    sod_cat[step]["fof_halo_tag"] = sod_cat[step]["fof_halo_tag"][srt]
    sod_cat[step]["sod_halo_mass"] = sod_cat[step]["sod_halo_mass"][srt]
    sod_cat[step]["sod_halo_radius"] = sod_cat[step]["sod_halo_radius"][srt]
    sod_cat[step]["sod_halo_cdelta"] = sod_cat[step]["sod_halo_cdelta"][srt]
    sod_cat[step]['x']=sod_cat[step]['x'][srt]
    sod_cat[step]['y']=sod_cat[step]['y'][srt]
    sod_cat[step]['z']=sod_cat[step]['z'][srt]
    srt = np.argsort(fof_cat[step]['fof_halo_tag'])
    fof_cat[step]['fof_halo_tag'] = fof_cat[step]['fof_halo_tag'][srt]
    fof_cat[step]['fof_halo_mass']= fof_cat[step]['fof_halo_mass'][srt]
    fof_cat[step]['x'] = fof_cat[step]['x'][srt]
    fof_cat[step]['y'] = fof_cat[step]['y'][srt]
    fof_cat[step]['z'] = fof_cat[step]['z'][srt]
    fof_cat[step]['vx'] = fof_cat[step]['vx'][srt]
    fof_cat[step]['vy'] = fof_cat[step]['vy'][srt]
    fof_cat[step]['vz'] = fof_cat[step]['vz'][srt]

    for i in range(0,gal_cat3[step]["host_halo_tag"].size):
        if(i%10000 == 0):
            print "i: ",i
        sod_m200, sod_r200, sod_c = get_sod_prop(gal_cat3[step]["host_halo_tag"][i],gal_cat3[step]["host_halo_mass"][i],stepz.get_z(step))
        #x,y,z,vx,vy,vz  = get_fof_prop(gal_cat3[step]["host_halo_tag"][i])
        x =gal_cat3[step]['host_halo_x'][i]
        y =gal_cat3[step]['host_halo_y'][i]
        z =gal_cat3[step]['host_halo_z'][i]
        vx =gal_cat3[step]['host_halo_vx'][i]
        vy =gal_cat3[step]['host_halo_vy'][i]
        vz =gal_cat3[step]['host_halo_vz'][i]

        r1 = rnd_NFWg_radius(NFW_maxR200,sod_c)
        r = r1*sod_r200
        ux,uy,uz = rnd_unit()  #random unit vector in 3D
        x = x + ux*r
        y = y + uy*r
        z = z + uz*r
        #Evrard et al 2007
        #https://arxiv.org/pdf/astro-ph/0702241.pdf
        a = stepz.get_a(step)
        omega_m = 0.223
        omega_l = 1-omega_m
        hz = 0.71*np.sqrt(omega_m/a**3 + omega_l)
        dispersion = 1082.9 * (hz*sod_m200/1e15)**0.3361
        vel = np.random.normal(scale=dispersion,size=(3))
        gal_cat3[step]['gltcs_x' ][i] =  x
        gal_cat3[step]['gltcs_y' ][i] =  y
        gal_cat3[step]['gltcs_z' ][i] =  z
        gal_cat3[step]['gltcs_vx'][i] = vx+vel[0]
        gal_cat3[step]['gltcs_vy'][i] = vy+vel[1]
        gal_cat3[step]['gltcs_vz'][i] = vz+vel[2]
    

hfile = h5py.File(output,"w")
hfile.create_dataset("gal_central/x", data=gal_central[step]["gltcs_x"])
hfile.create_dataset("gal_central/y", data=gal_central[step]["gltcs_y"])
hfile.create_dataset("gal_central/z", data=gal_central[step]["gltcs_z"])
hfile.create_dataset("gal_central/vx",data=gal_central[step]["gltcs_vx"])
hfile.create_dataset("gal_central/vy",data=gal_central[step]["gltcs_vy"])
hfile.create_dataset("gal_central/vz",data=gal_central[step]["gltcs_vz"])
hfile.create_dataset("gal_central/sdss_u",data=gal_central[step]["sdss_u"])
hfile.create_dataset("gal_central/sdss_g",data=gal_central[step]["sdss_g"])
hfile.create_dataset("gal_central/sdss_r",data=gal_central[step]["sdss_r"])
hfile.create_dataset("gal_central/sdss_i",data=gal_central[step]["sdss_i"])
hfile.create_dataset("gal_central/sdss_z",data=gal_central[step]["sdss_z"])
hfile.create_dataset("gal_central/stellar_mass",data=gal_central[step]["stellar_mass"])
hfile.create_dataset("gal_central/hostIndex",data=gal_central[step]["hostIndex"])
hfile.create_dataset("gal_central/nodeIndex",data=gal_central[step]["nodeIndex"])
hfile.create_dataset("gal_central/host_halo_tag",data=gal_central[step]["host_halo_tag"])
hfile.create_dataset("gal_central/host_halo_mass",data=gal_central[step]["host_halo_mass"])
hfile.create_dataset("gal_central/host_halo_x",data=gal_central[step]["host_halo_x"])
hfile.create_dataset("gal_central/host_halo_y",data=gal_central[step]["host_halo_y"])
hfile.create_dataset("gal_central/host_halo_z",data=gal_central[step]["host_halo_z"])
hfile.create_dataset("gal_central/infall_halo_tag",data=gal_central[step]["infall_halo_tag"])
hfile.create_dataset("gal_central/infall_halo_mass",data=gal_central[step]["infall_halo_mass"])


hfile.create_dataset("gal_core/x",data=gal_cores[step]["core_x"])
hfile.create_dataset("gal_core/y",data=gal_cores[step]["core_y"])
hfile.create_dataset("gal_core/z",data=gal_cores[step]["core_z"])
hfile.create_dataset("gal_core/vx",data=gal_cores[step]["core_vx"])
hfile.create_dataset("gal_core/vy",data=gal_cores[step]["core_vy"])
hfile.create_dataset("gal_core/vz",data=gal_cores[step]["core_vz"])
hfile.create_dataset("gal_core/sdss_u",data=gal_cores[step]["sdss_u"])
hfile.create_dataset("gal_core/sdss_g",data=gal_cores[step]["sdss_g"])
hfile.create_dataset("gal_core/sdss_r",data=gal_cores[step]["sdss_r"])
hfile.create_dataset("gal_core/sdss_i",data=gal_cores[step]["sdss_i"])
hfile.create_dataset("gal_core/sdss_z",data=gal_cores[step]["sdss_z"])
hfile.create_dataset("gal_core/stellar_mass",data=gal_cores[step]["stellar_mass"])
hfile.create_dataset("gal_core/hostIndex",data=gal_cores[step]["hostIndex"])
hfile.create_dataset("gal_core/nodeIndex",data=gal_cores[step]["nodeIndex"])
hfile.create_dataset("gal_core/host_halo_tag",data=gal_cores[step]["host_halo_tag"])
hfile.create_dataset("gal_core/host_halo_mass",data=gal_cores[step]["host_halo_mass"])
hfile.create_dataset("gal_core/host_halo_x",data=gal_cores[step]["host_halo_x"])
hfile.create_dataset("gal_core/host_halo_y",data=gal_cores[step]["host_halo_y"])
hfile.create_dataset("gal_core/host_halo_z",data=gal_cores[step]["host_halo_z"])
hfile.create_dataset("gal_core/infall_halo_tag",data=gal_cores[step]["infall_halo_tag"])
hfile.create_dataset("gal_core/infall_halo_mass",data=gal_cores[step]["infall_halo_mass"])

hfile.create_dataset("gal_bhp/x",data=gal_bhp[step]["bhp_x"])
hfile.create_dataset("gal_bhp/y",data=gal_bhp[step]["bhp_y"])
hfile.create_dataset("gal_bhp/z",data=gal_bhp[step]["bhp_z"])
hfile.create_dataset("gal_bhp/vx",data=gal_bhp[step]["bhp_vx"])
hfile.create_dataset("gal_bhp/vy",data=gal_bhp[step]["bhp_vy"])
hfile.create_dataset("gal_bhp/vz",data=gal_bhp[step]["bhp_vz"])
hfile.create_dataset("gal_bhp/sdss_u",data=gal_bhp[step]["sdss_u"])
hfile.create_dataset("gal_bhp/sdss_g",data=gal_bhp[step]["sdss_g"])
hfile.create_dataset("gal_bhp/sdss_r",data=gal_bhp[step]["sdss_r"])
hfile.create_dataset("gal_bhp/sdss_i",data=gal_bhp[step]["sdss_i"])
hfile.create_dataset("gal_bhp/sdss_z",data=gal_bhp[step]["sdss_z"])
hfile.create_dataset("gal_bhp/stellar_mass",data=gal_bhp[step]["stellar_mass"])
hfile.create_dataset("gal_bhp/hostIndex",data=gal_bhp[step]["hostIndex"])
hfile.create_dataset("gal_bhp/nodeIndex",data=gal_bhp[step]["nodeIndex"])
hfile.create_dataset("gal_bhp/host_halo_tag",data=gal_bhp[step]["host_halo_tag"])
hfile.create_dataset("gal_bhp/host_halo_mass",data=gal_bhp[step]["host_halo_mass"])
hfile.create_dataset("gal_bhp/host_halo_x",data=gal_bhp[step]["host_halo_x"])
hfile.create_dataset("gal_bhp/host_halo_y",data=gal_bhp[step]["host_halo_y"])
hfile.create_dataset("gal_bhp/host_halo_z",data=gal_bhp[step]["host_halo_z"])
hfile.create_dataset("gal_bhp/infall_halo_tag",data=gal_bhp[step]["infall_halo_tag"])
hfile.create_dataset("gal_bhp/infall_halo_mass",data=gal_bhp[step]["infall_halo_mass"])

hfile.create_dataset("gal_rnd/x",data=gal_rnd[step]["bhp_x"])
hfile.create_dataset("gal_rnd/y",data=gal_rnd[step]["bhp_y"])
hfile.create_dataset("gal_rnd/z",data=gal_rnd[step]["bhp_z"])
hfile.create_dataset("gal_rnd/vx",data=gal_rnd[step]["bhp_vx"])
hfile.create_dataset("gal_rnd/vy",data=gal_rnd[step]["bhp_vy"])
hfile.create_dataset("gal_rnd/vz",data=gal_rnd[step]["bhp_vz"])
hfile.create_dataset("gal_rnd/sdss_u",data=gal_rnd[step]["sdss_u"])
hfile.create_dataset("gal_rnd/sdss_g",data=gal_rnd[step]["sdss_g"])
hfile.create_dataset("gal_rnd/sdss_r",data=gal_rnd[step]["sdss_r"])
hfile.create_dataset("gal_rnd/sdss_i",data=gal_rnd[step]["sdss_i"])
hfile.create_dataset("gal_rnd/sdss_z",data=gal_rnd[step]["sdss_z"])
hfile.create_dataset("gal_rnd/stellar_mass",data=gal_rnd[step]["stellar_mass"])
hfile.create_dataset("gal_rnd/hostIndex",data=gal_rnd[step]["hostIndex"])
hfile.create_dataset("gal_rnd/nodeIndex",data=gal_rnd[step]["nodeIndex"])
hfile.create_dataset("gal_rnd/host_halo_tag",data=gal_rnd[step]["host_halo_tag"])
hfile.create_dataset("gal_rnd/host_halo_mass",data=gal_rnd[step]["host_halo_mass"])
hfile.create_dataset("gal_rnd/host_halo_x",data=gal_rnd[step]["host_halo_x"])
hfile.create_dataset("gal_rnd/host_halo_y",data=gal_rnd[step]["host_halo_y"])
hfile.create_dataset("gal_rnd/host_halo_z",data=gal_rnd[step]["host_halo_z"])
hfile.create_dataset("gal_rnd/infall_halo_tag",data=gal_rnd[step]["infall_halo_tag"])
hfile.create_dataset("gal_rnd/infall_halo_mass",data=gal_rnd[step]["infall_halo_mass"])

if(not NFW_profile):
    hfile.create_dataset("gal_nfw/x",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/y",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/z",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/vx",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/vy",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/vz",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/sdss_u",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/sdss_g",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/sdss_r",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/sdss_i",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/sdss_z",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/stellar_mass",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/hostIndex",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_nfw/nodeIndex",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_nfw/host_halo_tag",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_nfw/host_halo_mass",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/host_halo_x",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/host_halo_y",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_nfw/host_halo_z",data=np.zeros(0,dtype='f4;'))
    hfile.create_dataset("gal_nfw/infall_halo_tag",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_nfw/infall_halo_mass",data=np.zeros(0,dtype='f4'))


    hfile.create_dataset("gal_blank/x",data=gal_cat3[step]["gltcs_x"])
    hfile.create_dataset("gal_blank/y",data=gal_cat3[step]["gltcs_y"])
    hfile.create_dataset("gal_blank/z",data=gal_cat3[step]["gltcs_z"])
    hfile.create_dataset("gal_blank/vx",data=gal_cat3[step]["gltcs_vx"])
    hfile.create_dataset("gal_blank/vy",data=gal_cat3[step]["gltcs_vy"])
    hfile.create_dataset("gal_blank/vz",data=gal_cat3[step]["gltcs_vz"])
    hfile.create_dataset("gal_blank/sdss_u",data=gal_cat3[step]["sdss_u"])
    hfile.create_dataset("gal_blank/sdss_g",data=gal_cat3[step]["sdss_g"])
    hfile.create_dataset("gal_blank/sdss_r",data=gal_cat3[step]["sdss_r"])
    hfile.create_dataset("gal_blank/sdss_i",data=gal_cat3[step]["sdss_i"])
    hfile.create_dataset("gal_blank/sdss_z",data=gal_cat3[step]["sdss_z"])
    hfile.create_dataset("gal_blank/stellar_mass",data=gal_cat3[step]["stellar_mass"])
    hfile.create_dataset("gal_blank/hostIndex",data=gal_cat3[step]["hostIndex"])
    hfile.create_dataset("gal_blank/nodeIndex",data=gal_cat3[step]["nodeIndex"])
    hfile.create_dataset("gal_blank/host_halo_tag",data=gal_cat3[step]["host_halo_tag"])
    hfile.create_dataset("gal_blank/host_halo_mass",data=gal_cat3[step]["host_halo_mass"])
    hfile.create_dataset("gal_blank/host_halo_x",data=gal_cat3[step]["host_halo_x"])
    hfile.create_dataset("gal_blank/host_halo_y",data=gal_cat3[step]["host_halo_y"])
    hfile.create_dataset("gal_blank/host_halo_z",data=gal_cat3[step]["host_halo_z"])
    hfile.create_dataset("gal_blank/infall_halo_tag",data=gal_cat3[step]["infall_halo_tag"])
    hfile.create_dataset("gal_blank/infall_halo_mass",data=gal_cat3[step]["infall_halo_mass"])

else:
    hfile.create_dataset("gal_nfw/x",data=gal_cat3[step]["gltcs_x"])#note these gltcs_x positions
    hfile.create_dataset("gal_nfw/y",data=gal_cat3[step]["gltcs_y"])#are not from galaticus but from
    hfile.create_dataset("gal_nfw/z",data=gal_cat3[step]["gltcs_z"])#an NFW profile
    hfile.create_dataset("gal_nfw/vx",data=gal_cat3[step]["gltcs_vx"])
    hfile.create_dataset("gal_nfw/vy",data=gal_cat3[step]["gltcs_vy"])
    hfile.create_dataset("gal_nfw/vz",data=gal_cat3[step]["gltcs_vz"])
    hfile.create_dataset("gal_nfw/sdss_u",data=gal_cat3[step]["sdss_u"])
    hfile.create_dataset("gal_nfw/sdss_g",data=gal_cat3[step]["sdss_g"])
    hfile.create_dataset("gal_nfw/sdss_r",data=gal_cat3[step]["sdss_r"])
    hfile.create_dataset("gal_nfw/sdss_i",data=gal_cat3[step]["sdss_i"])
    hfile.create_dataset("gal_nfw/sdss_z",data=gal_cat3[step]["sdss_z"])
    hfile.create_dataset("gal_nfw/stellar_mass",data=gal_cat3[step]["stellar_mass"])
    hfile.create_dataset("gal_nfw/hostIndex",data=gal_cat3[step]["hostIndex"])
    hfile.create_dataset("gal_nfw/nodeIndex",data=gal_cat3[step]["nodeIndex"])
    hfile.create_dataset("gal_nfw/host_halo_tag",data=gal_cat3[step]["host_halo_tag"])
    hfile.create_dataset("gal_nfw/host_halo_mass",data=gal_cat3[step]["host_halo_mass"])
    hfile.create_dataset("gal_nfw/host_halo_x",data=gal_cat3[step]["host_halo_x"])
    hfile.create_dataset("gal_nfw/host_halo_y",data=gal_cat3[step]["host_halo_y"])
    hfile.create_dataset("gal_nfw/host_halo_z",data=gal_cat3[step]["host_halo_z"])
    hfile.create_dataset("gal_nfw/infall_halo_tag",data=gal_cat3[step]["infall_halo_tag"])
    hfile.create_dataset("gal_nfw/infall_halo_mass",data=gal_cat3[step]["infall_halo_mass"])

    hfile.create_dataset("gal_blank/x",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/y",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/z",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/vx",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/vy",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/vz",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/sdss_u",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/sdss_g",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/sdss_r",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/sdss_i",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/sdss_z",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/stellar_mass",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/hostIndex",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_blank/nodeIndex",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_blank/host_halo_tag",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_blank/host_halo_mass",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/host_halo_x",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/host_halo_y",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/host_halo_z",data=np.zeros(0,dtype='f4'))
    hfile.create_dataset("gal_blank/infall_halo_tag",data=np.zeros(0,dtype='i8'))
    hfile.create_dataset("gal_blank/infall_halo_mass",data=np.zeros(0,dtype='f4'))


x = np.concatenate((                gal_central[step]['gltcs_x'],         gal_cores[step]['core_x'],           gal_bhp[step]['bhp_x'],gal_cat3[step]['gltcs_x']))
y = np.concatenate((                gal_central[step]['gltcs_y'],         gal_cores[step]['core_y'],           gal_bhp[step]['bhp_y'],gal_cat3[step]['gltcs_y']))
z = np.concatenate((                gal_central[step]['gltcs_z'],         gal_cores[step]['core_z'],           gal_bhp[step]['bhp_z'],gal_cat3[step]['gltcs_z']))
vx = np.concatenate((               gal_central[step]['gltcs_vx'],        gal_cores[step]['core_vx'],          gal_bhp[step]['host_halo_vx'],gal_cat3[step]['gltcs_vx']))
vy = np.concatenate((               gal_central[step]['gltcs_vy'],        gal_cores[step]['core_vy'],          gal_bhp[step]['host_halo_vy'],gal_cat3[step]['gltcs_vy']))
vz = np.concatenate((               gal_central[step]['gltcs_vz'],        gal_cores[step]['core_vz'],          gal_bhp[step]['host_halo_vz'],gal_cat3[step]['gltcs_vz']))
sdss_u = np.concatenate((           gal_central[step]['sdss_u'],          gal_cores[step]['sdss_u'],           gal_bhp[step]['sdss_u'],gal_cat3[step]['sdss_u']))
sdss_g = np.concatenate((           gal_central[step]['sdss_g'],          gal_cores[step]['sdss_g'],           gal_bhp[step]['sdss_g'],gal_cat3[step]['sdss_g']))
sdss_r = np.concatenate((           gal_central[step]['sdss_r'],          gal_cores[step]['sdss_r'],           gal_bhp[step]['sdss_r'],gal_cat3[step]['sdss_r']))
sdss_i = np.concatenate((           gal_central[step]['sdss_i'],          gal_cores[step]['sdss_i'],           gal_bhp[step]['sdss_i'],gal_cat3[step]['sdss_i']))
sdss_z = np.concatenate((           gal_central[step]['sdss_z'],          gal_cores[step]['sdss_z'],           gal_bhp[step]['sdss_z'],gal_cat3[step]['sdss_z']))
stellar_mass = np.concatenate((     gal_central[step]['stellar_mass'],    gal_cores[step]['stellar_mass'],     gal_bhp[step]['stellar_mass'],gal_cat3[step]['stellar_mass']))
hostIndex = np.concatenate((        gal_central[step]['hostIndex'],       gal_cores[step]['hostIndex'],        gal_bhp[step]['hostIndex'],gal_cat3[step]['hostIndex']))
nodeIndex = np.concatenate((        gal_central[step]['nodeIndex'],       gal_cores[step]['nodeIndex'],        gal_bhp[step]['nodeIndex'],gal_cat3[step]['nodeIndex']))
host_halo_tag = np.concatenate((    gal_central[step]['host_halo_tag'],   gal_cores[step]['host_halo_tag'],    gal_bhp[step]['host_halo_tag'],gal_cat3[step]['host_halo_tag']))
host_halo_mass = np.concatenate((   gal_central[step]['host_halo_mass'],  gal_cores[step]['host_halo_mass'],   gal_bhp[step]['host_halo_mass'],gal_cat3[step]['host_halo_mass']))
host_halo_x = np.concatenate((      gal_central[step]['host_halo_x'],     gal_cores[step]['host_halo_x'],      gal_bhp[step]['host_halo_x'],gal_cat3[step]['host_halo_x']))
host_halo_y = np.concatenate((      gal_central[step]['host_halo_y'],     gal_cores[step]['host_halo_y'],      gal_bhp[step]['host_halo_y'],gal_cat3[step]['host_halo_y']))
host_halo_z = np.concatenate((      gal_central[step]['host_halo_z'],     gal_cores[step]['host_halo_z'],      gal_bhp[step]['host_halo_z'],gal_cat3[step]['host_halo_z']))
host_halo_vx = np.concatenate((     gal_central[step]['host_halo_vx'],    gal_cores[step]['host_halo_vx'],     gal_bhp[step]['host_halo_vx'],gal_cat3[step]['host_halo_vx']))
host_halo_vy = np.concatenate((     gal_central[step]['host_halo_vy'],    gal_cores[step]['host_halo_vy'],     gal_bhp[step]['host_halo_vy'],gal_cat3[step]['host_halo_vy']))
host_halo_vz = np.concatenate((     gal_central[step]['host_halo_vz'],    gal_cores[step]['host_halo_vz'],     gal_bhp[step]['host_halo_vz'],gal_cat3[step]['host_halo_vz']))
infall_halo_mass = np.concatenate(( gal_central[step]['infall_halo_mass'],gal_cores[step]['infall_halo_mass'], gal_bhp[step]['infall_halo_mass'],gal_cat3[step]['infall_halo_mass']))
infall_halo_tag = np.concatenate((  gal_central[step]['infall_halo_mass'],gal_cores[step]['infall_halo_tag'],  gal_bhp[step]['infall_halo_tag'],gal_cat3[step]['infall_halo_tag']))
gal_pos_type = np.concatenate((np.zeros_like(gal_central[step]['infall_halo_tag'],dtype=int)+-1,
                               np.zeros_like(gal_cores[step]['infall_halo_tag'],dtype=int)+0,
                               np.zeros_like(gal_bhp[step]['infall_halo_tag'],dtype=int)+1,
                               np.zeros_like(gal_cat3[step]['infall_halo_tag'],dtype=int)+3))



hfile.create_dataset("gal_all/x",data=x)
hfile.create_dataset("gal_all/y",data=y)
hfile.create_dataset("gal_all/z",data=z)
hfile.create_dataset("gal_all/vx",data=vx)
hfile.create_dataset("gal_all/vy",data=vy)
hfile.create_dataset("gal_all/vz",data=vz)
hfile.create_dataset("gal_all/sdss_u",data=sdss_u)
hfile.create_dataset("gal_all/sdss_g",data=sdss_g)
hfile.create_dataset("gal_all/sdss_r",data=sdss_r)
hfile.create_dataset("gal_all/sdss_i",data=sdss_i)
hfile.create_dataset("gal_all/sdss_z",data=sdss_z)
hfile.create_dataset("gal_all/stellar_mass",data=stellar_mass)
hfile.create_dataset("gal_all/hostIndex",data=hostIndex)
hfile.create_dataset("gal_all/nodeIndex",data=nodeIndex)
hfile.create_dataset("gal_all/host_halo_tag",data=host_halo_tag)
hfile.create_dataset("gal_all/host_halo_mass",data=host_halo_mass)
hfile.create_dataset("gal_all/host_halo_x",data=host_halo_x)
hfile.create_dataset("gal_all/host_halo_y",data=host_halo_y)
hfile.create_dataset("gal_all/host_halo_z",data=host_halo_z)
hfile.create_dataset("gal_all/host_halo_vx",data=host_halo_vx)
hfile.create_dataset("gal_all/host_halo_vy",data=host_halo_vy)
hfile.create_dataset("gal_all/host_halo_vz",data=host_halo_vz)
hfile.create_dataset("gal_all/infall_halo_tag",data=infall_halo_tag)
hfile.create_dataset("gal_all/infall_halo_mass",data=infall_halo_mass)
hfile.create_dataset("gal_all/pos_type",data=gal_pos_type)



hfile.close()
#slct_bhp = gal_bhp[step]['hostIndex']==-1
#slct_gal = gal_cat2[step]['hostIndex']==-1

#print "centrals with halo particles found: ", np.sum(slct_bhp)
#print "centrals with halo partilces not found", np.sum(slct_gal)

#print "sat with halo particles found: ", np.sum(slct_bhp==0)
#print "sat with halo particles not found: ", np.sum(slct_gal==0)


