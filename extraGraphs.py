import matplotlib.pyplot as plt
import numpy as np

sref = [0, 54.79166667, 92.12934524, 78.83928571, 180.5952381, 210.5059524, 343.5714286]
mref = [0, 93.56666404, 345.8340111, 769.9130264, 1214.661427, 2060.937142, 3012.751853]
tref = [0, 2.312783988, 6.13595752, -215.4845505, -207.5078057, -499.9967706, -486.356065]

wingtipdistance = [0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 18, 18, 20]

srefLastGrad = (sref[-1]-sref[-2])/(wingtipdistance[-4]-wingtipdistance[-5])
srefEnd6 = sref[-1] + (wingtipdistance[-3] - wingtipdistance[-4])*srefLastGrad
sref.append(srefEnd6)
sref.append(0)
sref.append(0)

mrefLastGrad = (mref[-1]-mref[-2])/(wingtipdistance[-4]-wingtipdistance[-5])
mrefEnd6 = mref[-1] + (wingtipdistance[-3] - wingtipdistance[-4])*mrefLastGrad
mref.append(mrefEnd6)
mref.append(mrefEnd6)
mref.append(mrefEnd6)

trefLastGrad = (tref[-1]-tref[-2])/(wingtipdistance[-4]-wingtipdistance[-5])
trefEnd6 = tref[-1] + (wingtipdistance[-3] - wingtipdistance[-4])*trefLastGrad
tref.append(trefEnd6)
tref.append(trefEnd6)
tref.append(trefEnd6)

# Plot Bending Moment Curve

plt.figure(1)
plt.plot(wingtipdistance, sref, 'k-', label='Shear Force - S$ref$')
plt.plot([18,18], [0, 435], 'k--', label='Fuselage Boundary')
plt.xlabel('Distance from Wing-Tip (m)')
plt.xlim([0, wingtipdistance[-1]])
plt.xticks(np.arange(0, 22, 3))
plt.ylim([0, max(sref)+30])
plt.ylabel('Shear Force (kN)')
plt.legend(loc='upper left', shadow=True)
plt.grid(True)
plt.title("Cumulative Shear Force over Loaded Wing at 2.5g")
plt.savefig('sref1.png')

plt.figure(2)
plt.plot(wingtipdistance, mref, 'k-', label='Bending Moment Force - M$ref$')
plt.plot([18,18], [0, 3640], 'k--', label='Fuselage Boundary')
plt.xlabel('Distance from Wing-Tip (m)')
plt.xlim([0, wingtipdistance[-1]])
plt.xticks(np.arange(0, 22, 3))
plt.ylim([0, max(mref)+160])
plt.ylabel('Bending Moment (kNm)')
plt.legend(loc='upper left', shadow=True)
plt.grid(True)
plt.title("Bending Moment over Loaded Wing at 2.5g")
plt.savefig('mref1.png')

plt.figure(3)
plt.plot(wingtipdistance, tref, 'k-', label='Torque - T$ref$')
plt.plot([18,18], [0, -520], 'k--', label='Fuselage Boundary')
plt.xlim([0, wingtipdistance[-1]])
plt.xlabel('Distance from Wing-Tip (m)')
plt.xticks(np.arange(0, 22, 3))
plt.ylabel('Torque (kNm)')
plt.legend(loc='upper right', shadow=True)
plt.grid(True)
plt.title("Torque through Loaded Wing at 2.5g")
plt.savefig('tref1.png')

plt.show()
