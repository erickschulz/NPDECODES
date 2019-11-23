from numpy import array
import matplotlib.pyplot as plt
import csv
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

with open(input_file, 'r') as inputFile:
	reader = csv.reader(inputFile, delimiter=',')
	rows = list(reader)
	t = array([float(ri) for ri in rows[0]])
	E_pot = array([float(ri) for ri in rows[1]])
	E_kin = array([float(ri) for ri in rows[2]])
inputFile.close()

# converts an array of size m+1 to size m
Reduce = lambda x: array([(x[i] + x[i + 1]) * 0.5 for i in range(len(x) - 1)])

t_reduced = Reduce(t)
E_pot_reduced = Reduce(E_pot)
E_tot = E_pot_reduced + E_kin

fig = plt.figure()
plt.plot(t, E_pot, label='potential energy')
plt.plot(t_reduced, E_kin, '--', label='kinetic energy')
plt.plot(t_reduced, E_tot, label='total energy')
plt.xlabel('time t')
plt.ylabel('energy')
plt.legend()
plt.savefig(output_file)

print("Generated ", output_file)
