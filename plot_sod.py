import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('../output/Sod_3D/1.h5', 'r')
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

fig = plt.figure()
ax1 = fig.add_subplot(2,2,1)
ax1.plot(d[0,0,:], 'o')
plt.axis([0, 100, 0, 1.1])
plt.ylabel('density')
ax2 = fig.add_subplot(2,2,2)
plt.axis([0, 100, -0.1, 1.1])
ax2.plot(vz[0,0,:], 'o')
plt.ylabel('velocity')
ax3 = fig.add_subplot(2,2,3)
ax3.plot(p[0,0,:], 'o')
plt.ylabel('pressure')
plt.axis([0, 100, 0, 1.1])
ax4 = fig.add_subplot(2,2,4)
ax4.plot(e[0,0,:], 'o')
plt.ylabel('internal energy')
plt.axis([0, 100, 1.5, 3.7])

plt.show()
