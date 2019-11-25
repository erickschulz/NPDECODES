from csv import reader
from matplotlib.pyplot import figure, legend, plot, savefig, xlabel, ylabel
from numpy import array
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

with open(input_file, 'r') as inputFile:
	rows = list(reader(inputFile, delimiter=','))
	t = array([float(ri) for ri in rows[0]])
	E_pot = array([float(ri) for ri in rows[1]])
	E_kin = array([float(ri) for ri in rows[2]])
inputFile.close()

# converts an array of size m+1 to size m
Average = lambda x: array([(x[i] + x[i + 1]) * 0.5 for i in range(len(x) - 1)])

t_averaged = Average(t)
E_pot_averaged = Average(E_pot)
E_tot = E_pot_averaged + E_kin

fig = figure()
plot(t, E_pot, label='potential energy')
plot(t_averaged, E_kin, '--', label='kinetic energy')
plot(t_averaged, E_tot, label='total energy')
xlabel('time t')
ylabel('energy')
legend()
savefig(output_file)

print('Generated ' + output_file)
