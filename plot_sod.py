import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('./hdf5/1.h5', 'r')
head = f.attrs
gamma = head['gamma'][0]
d  = np.array(f['density'])
mx = np.array(f['momentum_x'])
my = np.array(f['momentum_y'])
mz = np.array(f['momentum_z'])
E  = np.array(f['Energy'])
#e  = np.array(f['GasEnergy'])
vx = mx/d
vy = my/d
vz = mz/d
p  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0)
e  = p/d/(gamma - 1.0)
#e  = e/d

f2 = open('exact.txt', 'r')
f2.readline()
lines=f2.readlines()
f2.close()

x = []
de = []
ve = []
pe = []
ie = []

for line in lines:
  p = line.split()
  x.append(float(p[0]))
  de.append(float(p[1]))
  ve.append(float(p[2]))
  pe.append(float(p[3]))
  ie.append(float(p[4]))

x = np.array(x)
de = np.array(de)
ve = np.array(ve)
pe = np.array(pe)
ie = np.array(ie)

fig = plt.figure()
ax1 = plt.axes([0.1, 0.6, 0.35, 0.35])
plt.axis([0, 100, 0, 1.1])
ax1.plot(d, 'o', markersize=3, color='black')
ax1.plot(x*100, de, color='black')
plt.ylabel('density')
ax2 = plt.axes([0.6, 0.6, 0.35, 0.35])
plt.axis([0, 100, -0.1, 1.1])
ax2.plot(vx, 'o', markersize=3, color='black')
ax2.plot(x*100, ve, color='black')
plt.ylabel('velocity')
ax3 = plt.axes([0.1, 0.1, 0.35, 0.35])
plt.axis([0, 100, 0, 1.1])
ax3.plot(p, 'o', markersize=3, color='black')
ax3.plot(x*100, pe, color='black')
plt.ylabel('pressure')
ax4 = plt.axes([0.6, 0.1, 0.35, 0.35])
plt.axis([0, 100, 0.5, 2.5])
ax4.plot(e, 'o', markersize=3, color='black')
ax4.plot(x*100, ie, color='black')
plt.ylabel('internal energy')

plt.savefig("sod.png", dpi=300);
plt.close(fig)
