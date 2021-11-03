import sys
import matplotlib.pyplot as plt
import matplotlib.colors as col
import numpy as np


output_file = str(sys.argv[1])
if len(sys.argv) == 3:
    input_gamma = sys.argv[2]
else:
    input_gamma = 1.0;

# Stability function
S = lambda z, gamma=1.0: (1.0 + z * (1.0 - 2.0 * gamma) + z**2 * (gamma**2 - 2.0 * gamma + 0.5)) / (1.0 - gamma * z)**2

absS = lambda z: np.abs(S(z, gamma=input_gamma))

# Compute F(x) on a meshgrid
grid1D = np.linspace(-7.0, 7.0, 180, endpoint=True)
X, Y = np.meshgrid(grid1D, grid1D)
Z = absS(X + 1.0j * Y)

# Contour plot distinguishes absS < 1 and absS > 1
fig = plt.figure()
cmap = col.ListedColormap(['lime','w'])
bounds = [0.0, 1.0, 2.0]
norm = col.BoundaryNorm(bounds, cmap.N)
cs1 = plt.contourf(X, Y, Z, cmap=cmap, norm=norm, levels=bounds, extend='both')
linewidth = 0.2
plt.contour(X, Y, Z, colors='k', levels=[0.0, 1.0], linewidths=linewidth)
plt.plot([-6.0, 6.0], [0.0, 0.0], color='k', linewidth=linewidth)
plt.plot([0.0, 0.0], [-6.0, 6.0], color='k', linewidth=linewidth)
plt.xlabel('Re')
plt.ylabel('Im')
plt.axis('square')
#fig.colorbar(cs1)
plt.savefig(output_file, bbox_inches='tight')
