import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

tR = np.genfromtxt("solution.csv", delimiter=',')
t = tR[:,0]
R = tR[:,1:]

m, n = R.shape
x = np.linspace(0, 1, n)
X, T = np.meshgrid(x, t)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, T, R, cmap='cool', edgecolor='black', alpha=0.5)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u(x,t)')
plt.savefig("viswave.png")

print("The plot has been written to viswave.png.")
