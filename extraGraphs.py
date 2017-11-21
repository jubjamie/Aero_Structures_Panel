import matplotlib.pyplot as plt
import numpy as np

sref = [0, 54.79166667, 92.12934524, 78.83928571, 180.5952381, 210.5059524, 343.5714286, 0]
mref = [0, 93.56666404, 345.8340111, 769.9130264, 1214.661427, 2060.937142, 3012.751853]
tref = [0, 2.312783988, 6.13595752, -215.4845505, -207.5078057, -499.9967706, -486.356065]
wingtipdistance = [0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 20]

# Plot Bending Moment Curve

plt.figure(1)
plt.plot(wingtipdistance, sref, 'b-', label='Shear Force - S$ref$')
plt.xlabel('Distance from Wing-Tip (m)')
plt.ylabel('Shear Force (kN)', color='b')
plt.tick_params('y', colors='b')
plt.legend(loc='upper left', shadow=True)
plt.grid(True)
plt.title("Cumulative Shear Force over Loaded Wing at 2.5g")

plt.show()
