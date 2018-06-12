#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import h5py 
import sys

import dtk
from util import *


#param = dtk.Param(sys.argv[1])
output_gio_pos = "output/gal_499_pos.gio"
output_gio = "output/gal_499.gio"

gal_id = dtk.gio_read(output_gio,"nodeIndex") 
gal_x  = dtk.gio_read(output_gio,"x") 
gal_y  = dtk.gio_read(output_gio,"y") 
gal_z  = dtk.gio_read(output_gio,"z") 
mask  = 0x000FFFFF00000000
mask2 = 0x0FF0000000000000
print mask,mask2,gal_id.dtype
gal_id2 = (gal_id & mask  )>>32
gal_id3 = (gal_id & mask2  )>>52
print gal_id2[:10]
print gal_id3[:10]
target = 2247296544770364691
slct = gal_id == target
print gal_id[:10]
print gal_id.size
print np.sum(slct)
