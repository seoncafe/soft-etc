#!/usr/bin/env python
import numpy as np
from fbm_lib import *

#seed   = None
seed   = 9876540
nx     = 1024
ny, nz = nx,nx
mach   = 2

kcut_list = [0,1,2,3,4,5]

for kcut in kcut_list:
   fname = 'M%03d_%04d_kcut%d.fits.gz' % (mach*10, nx, kcut)
   print('making...',fname)
   a = fbm3d_ISM(nx,ny,nz,mach=mach,k_cut_index=kcut)
   a.writeto(fname)
