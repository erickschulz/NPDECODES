from matplotlib.pyplot import figure, legend, plot, savefig, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
t = data[0]
energy = data[1]

fig = figure()
plot(t, energy, '-')#, label='T = 0.3')
xlabel('t')
ylabel('H(u(t))')
#legend()
savefig(output_file)

print('Generated ' + output_file)
