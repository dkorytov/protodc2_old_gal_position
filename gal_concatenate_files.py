#!/usr/bin/env python2.7

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import dtk
import h5py
import sys


input_locs = []
output_loc = "tmp.hdf5"


groups = ['gal_bhp','gal_blank','gal_central','gal_core','gal_nfw','gal_rnd']
columns = ['nodeIndex','x','y','z','vx','vy','vz']


output_hfile = h5py.File(output_loc,'w')
for group in groups:
    print("working on group ",group)
    output_group = output_hfile[group]
    for column in columns:
        print("\tcolumn")
        data_list = []
        for input_loc in input_locs:
            print("\t\tinput_loc")
            hgroup = h5py.File(input_loc,'r')[group]
            data_list.append(hgroup[column].value)
        data = np.concatenate(data_list)
        output_group[column] = data
        output_hfile.flush()
        
        
