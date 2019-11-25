from matplotlib.pyplot import figure, savefig
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt, linspace, meshgrid
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

tR = genfromtxt(input_file, delimiter=',')
t = tR[:,0]
R = tR[:,1:]

m, n = R.shape
x = linspace(0, 1, n)
X, T = meshgrid(x, t)

fig = figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, T, R, cmap='cool', edgecolor='black', alpha=0.5)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u(x,t)')
savefig(output_file)

print('Generated ' + output_file)
