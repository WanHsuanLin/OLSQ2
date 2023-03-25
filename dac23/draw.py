import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams.update({'font.size': 14})
# width = 23
# height = 12
# linewidth = 8 #default 1.5
# markersize = 20

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
x_label = ['', r"5$\times$5", r'6$\times$6', r'7$\times$7', r'8$\times$8', r'9$\times$9']
x = np.arange(5, 10)
y = np.arange(15, 39, 3)
X, Y = np.meshgrid(x, y)

# Z = np.array([[6.61, 9.35, 12.62, 12.36, 17.16, 19.72, 20.95, 104.13],
            #  [11.04, 13.96, 19.54, 21.1, 27.08, 37.84, 22.87, 29.84],
            #  [16.85, 22.21, 26.45, 16.76, 61.93, 29.66, 40.15, 29.86],
            #  [21.2, 29.68, 20.48, 22.14, 34.6, 248.14, 221.73, 213.48],
            #  [19.83, 23.67, 30.68, 138.04, 75.82, 227.32, 558.38, 508.34]])

Z = np.array([[51.61, 324.46, 731.07, 904.28, 2073.51, 1977.53, 6584.44, 13609.38],
[244.18, 875.58, 1019.2, 1361.78, 6213.32, 9844.04, 12606.27, 27478.22],
[623.01, 1562.69, 2720.55, 4352.07, 13367.83, 13004.06, 27073.6, 69496.54],
[1497.4, 3616.37, 6483.33, 11879.47, 25830.72, 84290.19, 2500*60, 2500*60],
[2059.87, 12161.26, 21120.41, 43288.48, 2500*60, 2500*60, 2500*60, 2500*60]])

Z = Z/60
Z = Z.T

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, linewidth=2, cmap='cividis', vmin = 0, vmax = 2550)
# surf = ax.plot_surface(X, Y, Z, linewidth=0, antialiased=False, cmap='cividis')

# Customize the z axis.
ax.set_zlim(0, 2550)
ax.set_xlim(4, 10)
ax.set_xticklabels(x_label)
# ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.0f}')

plt.xlabel("Grid Size")
plt.ylabel("Gate Number")
ax.set_zlabel("Z3 Solving Time (min)")

plt.show()
# plt.savefig("runtime_new_olsq.pdf", dpi=300, bbox_inches='tight')