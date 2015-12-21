import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

f = h5py.File('../output/Noh_3D/1.h5', 'r')
head = f.attrs
gamma = head['gamma'][0]
d  = f['density']
mx = f['momentum_x']
my = f['momentum_y']
mz = f['momentum_z']
E  = f['Energy']
#e  = f['GasEnergy']
d  = np.array(d)
mx = np.array(mx)
my = np.array(my)
mz = np.array(mz)
E  = np.array(E)
#e  = np.array(e)
vx = mx/d
vy = my/d
vz = mz/d
p  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0)
e  = p/d/(gamma - 1.0)
#e  = e/d

levels = np.linspace(4.0, 64.0, 31) 
x = np.linspace(0.0, 1.0, 64, endpoint=False)
y = np.linspace(0.0, 1.0, 64, endpoint=False)
dmin, dmax = np.min(d), np.max(d)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
plt.pcolor(x, y, d[:,:,0], cmap=cm.jet, vmin=dmin, vmax=dmax)
plt.colorbar(label='Density')
plt.axis([0.0, 1.0, 0.0, 1.0])
plt.contour(x, y, d[:,:,0], levels, colors='black')
plt.show()
#fig.savefig('implosion.png')
