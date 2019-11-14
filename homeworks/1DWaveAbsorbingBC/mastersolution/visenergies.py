import numpy as np
import matplotlib.pyplot as plt
import csv

with open('energies.csv', 'r') as inputFile:
	reader = csv.reader(inputFile, delimiter=',')
	rows = list(reader)
	t = np.array([float(ri) for ri in rows[0]])
	E_pot = np.array([float(ri) for ri in rows[1]])
	E_kin = np.array([float(ri) for ri in rows[2]])
inputFile.close()

# converts an array of size m+1 to size m
Reduce = lambda x: np.array([(x[i] + x[i + 1]) * 0.5 for i in np.arange(len(x) - 1)])

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
plt.savefig("visenergies.png")

print("The plot has been written to visenergies.png.")
