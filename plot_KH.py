import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.colors import hsv_to_rgb

N = 100
dname='PPMP_CTU'

for i in range(0,N+1):

  f = h5py.File(dname+'/hdf5/'+str(i)+'.h5', 'r')
  head = f.attrs
  gamma = head['gamma'][0]
  nx = head['dims'][0]
  ny = head['dims'][1]
  d  = np.array(f['density'])
  mx = np.array(f['momentum_x'])
  my = np.array(f['momentum_y'])
  mz = np.array(f['momentum_z'])
  E  = np.array(f['Energy'])
  vx = mx/d
  vy = my/d
  vz = mz/d
  p  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0)
  e  = p/d/(gamma - 1.0)

  min_d = 0.01
  max_d = 2.0
  dpo = np.clip(d, min_d, max_d)
  dpo  = (dpo-min_d)/(max_d-min_d)

  min_m = -2.0
  max_m = 2.0
  mpo = np.clip(mx, min_m, max_m)
  mpo = (mpo-min_m)/(max_m-min_m)

  cmin = 120./360
  cmax = 300./360

  H = (cmax-cmin)*mpo + cmin
  S = (1.0 - 0.75*dpo)
  V = dpo

  HSV = np.dstack((H,S,V))
  RGB = hsv_to_rgb(HSV)

  rotate = ndimage.rotate(RGB, 90)

  fig = plt.figure(figsize=(nx/100., ny/100.), dpi=100, frameon=False)
  a0 = plt.axes([0,0,1,1])
  a0.text(20, 240, 'PPMP with CTU', color='white')
  plt.imshow(rotate, origin='lower')
  plt.axis('off')
  plt.savefig(dname+'/png/'+str(i)+'.png', dpi=100)
  plt.close(fig)


