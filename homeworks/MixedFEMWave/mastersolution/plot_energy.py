from matplotlib.pyplot import figure, plot, savefig, title, xlabel, ylabel
from numpy import genfromtxt, linspace
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
t = data[0]
energy = data[1]

fig = figure()
title('Evolution of the energy')
plot(t, energy, '-')
xlabel('t')
ylabel('E')
savefig(output_file)

print('Generated ' + output_file)
