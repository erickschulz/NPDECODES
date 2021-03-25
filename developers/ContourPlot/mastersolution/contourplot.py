import sys
import matplotlib.pyplot as plt
import numpy as np

# Read .csv file
input_file = str(sys.argv[1])
output_file = str(sys.argv[2])
isoline1 = np.genfromtxt(input_file, delimiter=',', skip_footer=2)
isoline2 = np.genfromtxt(input_file, delimiter=',', skip_header=2)

# Compute F(x) on a meshgrid
grid1D = np.linspace(-2.0, 4.0, 100)
X, Y = np.meshgrid(grid1D, grid1D)
F = lambda x, y: (x**2 + y**2)**2 - x**3 - y**3
Z = F(X, Y)

# Plot isolines on top of F(x)
fig, ax = plt.subplots()
cs = ax.contourf(X, Y, Z, 50, cmap=plt.cm.ocean_r)
cbar = fig.colorbar(cs)
plt.plot(isoline1[0], isoline1[1], '--r', label=r'Using $\mathbf{grad}\,F,\ \mathbf{y}_0 = (1, 0)$')
plt.plot(isoline2[0], isoline2[1], '-g', label=r'Using $F,\ \mathbf{y}_0 = (2, 0)$')
plt.legend(loc='upper right', framealpha=1.0)
plt.savefig(output_file, bbox_inches='tight')
