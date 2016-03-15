import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 

m_h = 1.672622e-24 # in grams
d_c = m_h
l_c = 3.08567758e18   # cm in 1 pc
t_c = 3.15569e10      # seconds in a kyr
v_c = l_c / t_c
p_c = d_c * v_c * v_c
kb = 1.380658e-16     # k_boltzmann in ergs/K
m_c = d_c * l_c * l_c * l_c / 1.9891e33 # characteristic mass in solar masses
P_init = 2.744870e-12

f = open('heat_eq_cool.txt', 'r')
lines = f.readlines()
f.close()

neq = []
Teq = []

for line in lines:
  p = line.split()
  neq.append(float(p[0]))
  Teq.append(float(p[1]))

neq = np.array(neq)
Teq = np.array(Teq)
T_a = P_init / (pow(10,neq)*kb)
T_a = np.log10(T_a)

dname='sphere_wind/n01/'
istart=0
iend=100
scale=5
for i in range(istart,iend+1):

  f = h5py.File(dname+'hdf5/'+str(i)+'.h5', 'r')
  head = f.attrs
  gamma = head['gamma'][0]
  nx = head['dims'][0]
  ny = head['dims'][1]
  nz = head['dims'][2]
  d  = np.array(f['density']).astype(float)
  mx = np.array(f['momentum_x']).astype(float)
  my = np.array(f['momentum_y']).astype(float)
  mz = np.array(f['momentum_z']).astype(float)
  E  = np.array(f['Energy']).astype(float)
  GE = np.array(f['GasEnergy']).astype(float)
  vx = mx/d
  vy = my/d
  vz = mz/d
  #p  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0)
  p  = GE * (gamma - 1.0)
  #e  = p/d/(gamma - 1.0)
  e  = GE/d

  n = d*d_c/m_h # number density
  log_n = np.log10(n)
  T = (p*p_c)/(n*kb) # temperature
  log_T = np.log10(T)

  n_hist = np.reshape(log_n, nx*ny*nz)
  v_hist = np.reshape(vx*v_c/1e5, nx*ny*nz)
  T_hist = np.reshape(log_T, nx*ny*nz)

  fig1 = plt.figure()
  plt.hist2d(n_hist, T_hist, bins=100, range=[[-5,3],[1,8]], norm=LogNorm())
  plt.colorbar()
  plt.plot(neq, Teq, color='black', linestyle='--')
  plt.plot(neq, T_a, color='black')
  plt.xlabel('$log_{10}(n)$ [$cm^{-3}$]')
  plt.ylabel('$log_{10}(T)$ [K]')
  plt.text(1,7.5,str(scale*i)+' kyr')
  fig1.savefig(dname+'phase_diagrams/'+'nT_'+str(i)+'.png')
  plt.close(fig1)

  fig2 = plt.figure()
  plt.hist2d(n_hist, v_hist, bins=100, range=[[-5,3],[0,1500]], norm=LogNorm())
  plt.colorbar()
  plt.xlabel('$log_{10}(n)$ [$cm^{-3}$]')
  plt.ylabel('$v$ [km $s^{-1}$]')
  plt.text(1,1300,str(scale*i)+' kyr')
  fig2.savefig(dname+'phase_diagrams/'+'nv_'+str(i)+'.png')
  plt.close(fig2)

  fig3 = plt.figure()
  plt.plot(neq, Teq, color='black', linestyle='--')
  plt.plot(neq, T_a, color='black')
  plt.axis([-5,3,1,8])
  h1, xe1, ye1 = np.histogram2d(n_hist, T_hist, bins=100, weights=v_hist, range=[[-5,3],[1,8]])
  h2, xe2, ye2 = np.histogram2d(n_hist, T_hist, bins=100, range=[[-5,3],[1,8]])
  h = h1/h2
  h = h.T
  plt.imshow(h, origin='lower', extent=[-5,3,1,8], interpolation='nearest')
  plt.colorbar(label='v [km $s^{-1}$]')
  plt.xlabel('$log_{10}(n)$ [$cm^{-3}$]')
  plt.ylabel('$log_{10}(T)$ [K]')
  plt.text(1,7.5,str(scale*i)+' kyr')
  fig3.savefig(dname+'phase_diagrams/'+'nvT_'+str(i)+'.png')
  plt.close(fig3)

