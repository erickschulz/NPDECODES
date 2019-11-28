from matplotlib.pyplot import figure, savefig, cm
from mesh_reader import mesh_reader
from mpl_toolkits.mplot3d import Axes3D
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

x, y, z, triangles = mesh_reader(input_file)


fig = figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
ax.plot_trisurf(x, y, triangles, z, linewidth=0.2, antialiased=True, cmap=cm.Spectral)
fig.savefig(output_file, bbox_inches='tight')
