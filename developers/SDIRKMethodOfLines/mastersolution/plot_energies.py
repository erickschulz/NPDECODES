from matplotlib.pyplot import figure, plot, savefig, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
t = data[0]
energies = data[1]

figure()
plot(t, energies)
xlabel(r'time $t$')
ylabel(r'thermal energy $E_h(t)$')
savefig(output_file)

print('Generated ' + output_file)
