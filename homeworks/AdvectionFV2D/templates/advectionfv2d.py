import numpy as np
import matplotlib.pyplot as plt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = np.genfromtxt(input_file, delimiter=',')

num_cells = data[:,0]
l2_errors = data[:,1]

fig = plt.figure()
plt.plot(num_cells, l2_errors, 'o-')
plt.grid()
plt.xlabel("Number of Cells")
plt.xscale('log', basex=2)
plt.ylabel("L2 Error")
plt.yscale('log', basey=2)
plt.tight_layout()

plt.savefig(output_file)
print('Generated ' + output_file)
