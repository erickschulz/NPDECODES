from matplotlib.pyplot import figure, title, triplot, savefig, xlabel, ylabel
from mesh_reader import mesh_reader
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

x, y, _, triangles = mesh_reader(input_file)

fig = figure(0, figsize=(10, 10))
title(input_file.split('/')[-1])
triplot(x, y, triangles, linewidth=0.5)
xlabel('x')
ylabel('y')
fig.savefig(output_file, bbox_inches='tight')
