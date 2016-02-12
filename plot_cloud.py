import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.colors import hsv_to_rgb
from matplotlib.patches import Rectangle

m_h = 1.672622e-24 # in grams
d_c = m_h
l_c = 3.08567758e18   # cm in 1 pc
t_c = 3.15569e10      # seconds in a kyr
v_c = l_c / t_c
p_c = d_c * v_c * v_c
kb = 1.380658e-16     # k_boltzmann in ergs/K
m_c = d_c * l_c * l_c * l_c / 1.9891e33 # characteristic mass in solar masses
istart=0
iend=250
dname='n1/'

for i in range(istart,iend+1):

  print(i)
  f = h5py.File(dname+'hdf5/'+str(i)+'.h5', 'r')
  head = f.attrs
  gamma = head['gamma'][0]
  nx = head['dims'][0]
  ny = head['dims'][1]
  nz = head['dims'][2]
  d  = np.array(f['density'])
  #mx = np.array(f['momentum_x'])
  #my = np.array(f['momentum_y'])
  #mz = np.array(f['momentum_z'])
  #E  = np.array(f['Energy'])
  GE = np.array(f['GasEnergy'])
  #vx = mx/d
  #vy = my/d
  #vz = mz/d
  #p  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0)
  p = GE * (gamma-1.0)
  #e  = p/d/(gamma - 1.0)
  e  = GE/d

  n = d*d_c/m_h # number density
  T =  (p*p_c)/(n*kb) # temperature
  T = np.clip(T, 1, 1e9)
  log_T = np.log10(T)

  log_d = np.log10(d)

  d_cloud = np.mean(d[d>0.001])
  d_med = np.median(d[d>0.001])
  m_cloud = (np.sum(d[d>0.001])*(20./128)**3)*m_c
  print("cloud density: %8.5f" %d_cloud)
  print("cloud mass:    %8.5f" %m_cloud)
  print("median density:%8.5f" %d_med)
  print("d  range     = %8.5f %11.5f" % (np.min(d),np.max(d)))
  print("T  range     = %e %e" % (np.min(T),np.max(T)))
  print("log d range  = %8.5f %11.5f" % (np.min(log_d),np.max(log_d)))
  print("log T range  = %e %e" % (np.min(log_T),np.max(log_T)))



  #make color image
  cmin=0/360.
  cmax=300/360.
  #set the temperature scale
  min_T = 2.0 
  max_T = 8.0
  cT = (np.clip(log_T, min_T, max_T) - min_T + 0.01) / (max_T-min_T)
  H = (cmin-cmax)*cT + cmax
  #print np.min(H), np.max(H)

  #set the density scale
  d_min = 0.0001
  d_max = 150 

  #produce baseline density plot
  rpd = (log_d - np.log10(d_min))/(np.log10(d_max) - np.log10(d_min))
  V = rpd

  #eliminate accidental negative values
  V = np.clip(V, 0.01, 0.99)
  S = 1.0 - V

  HSV = np.dstack((H,S,V))
  RGB = hsv_to_rgb(HSV)

  rotate = ndimage.rotate(RGB, 90)

  fig = plt.figure(figsize=(0.960*10, 0.256*10), dpi=100, frameon=False)
  #fig = plt.figure(figsize=(1.920*10, 0.256*10), dpi=100, frameon=False)
  #fig = plt.figure(figsize=(1.920*10, 0.768*10), dpi=100, frameon=False)
  #fig = plt.figure(figsize=(0.640*10, 0.320*10), dpi=100, frameon=False)
  a0 = plt.axes([0,0,1,1])
  plt.imshow(rotate, origin='lower')
  plt.axis('off')
  a0.text(50, 220, str(10*i)+' kyr', color='white')
  plt.savefig(dname+'png/'+str(i)+'.png', dpi=100)
  plt.close(fig)


