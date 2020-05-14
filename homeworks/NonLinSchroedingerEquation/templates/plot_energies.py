from matplotlib.pyplot import figure, legend, plot, savefig, title, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
t = data[0]
E_kin = data[1]
E_int = data[2]

fig = figure()
title('Energies Along the Approximate Solution')
plot(t, E_kin + E_int, '-', label='total energy')
plot(t, E_kin, '-', label='kinetic energy')
plot(t, E_int, '-', label='interaction energy')
xlabel('t')
ylabel('H(u(t))')
legend()
savefig(output_file)

print('Generated ' + output_file)
