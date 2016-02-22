import h5py
import numpy as np

ns = 0 
ne = 400
n_procs = 8 

for n in range(ns, ne+1):

  dname = './n10_big/hdf5/'
  f = h5py.File(dname+str(n)+'.h5', 'w')

  file = h5py.File(dname+'raw/'+str(n)+'.h5.0', 'r')
  head = file.attrs
  nx = head['dims'][0]
  ny = head['dims'][1]
  nz = head['dims'][2]
  nxl = head['dims_local'][0]
  nyl = head['dims_local'][1]
  nzl = head['dims_local'][2]
  xs = head['offset'][0]
  ys = head['offset'][1]
  zs = head['offset'][2]
  f.attrs['dims'] = [nx, ny, nz]
  f.attrs['gamma'] = [head['gamma'][0]]
  f.attrs['t'] = [head['t'][0]]
  f.attrs['dt'] = [head['dt'][0]]
  f.attrs['n_step'] = [head['n_step'][0]]

  d  = f.create_dataset("density", (nx, ny))
  mx = f.create_dataset("momentum_x", (nx, ny))
  my = f.create_dataset("momentum_y", (nx, ny))
  mz = f.create_dataset("momentum_z", (nx, ny))
  E  = f.create_dataset("Energy", (nx, ny))
  GE  = f.create_dataset("GasEnergy", (nx, ny))

  d[xs:xs+nxl,ys:ys+nyl]  = file['density']
  mx[xs:xs+nxl,ys:ys+nyl] = file['momentum_x']
  my[xs:xs+nxl,ys:ys+nyl] = file['momentum_y']
  mz[xs:xs+nxl,ys:ys+nyl] = file['momentum_z']
  E[xs:xs+nxl,ys:ys+nyl]  = file['Energy']
  GE[xs:xs+nxl,ys:ys+nyl] = file['GasEnergy']

  for i in range(1, n_procs):

    file = h5py.File(dname+'raw/'+str(n)+'.h5.'+str(i), 'r')
    head = file.attrs
    nxl = head['dims_local'][0]
    nyl = head['dims_local'][1]
    nzl = head['dims_local'][2]
    xs = head['offset'][0]
    ys = head['offset'][1]
    zs = head['offset'][2]

    d[xs:xs+nxl,ys:ys+nyl]  = file['density']
    mx[xs:xs+nxl,ys:ys+nyl] = file['momentum_x']
    my[xs:xs+nxl,ys:ys+nyl] = file['momentum_y']
    mz[xs:xs+nxl,ys:ys+nyl] = file['momentum_z']
    E[xs:xs+nxl,ys:ys+nyl]  = file['Energy']
    GE[xs:xs+nxl,ys:ys+nyl] = file['GasEnergy']

