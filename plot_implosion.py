import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

f = h5py.File('../output/implosion_2D/1.h5', 'r')
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

levels = np.linspace(0.125, 1.0, 36) 
x = np.linspace(0.0, 0.3, 400, endpoint=False)
y = np.linspace(0.0, 0.3, 400, endpoint=False)
pmin, pmax = np.min(p), np.max(p)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
plt.pcolor(x, y, p, cmap=cm.jet, vmin=pmin, vmax=pmax)
plt.colorbar(label='Pressure')
plt.axis([0.0, 0.3, 0.0, 0.3])
plt.contour(x, y, d, levels, colors='black')
plt.show()
fig.savefig('implosion.png')
