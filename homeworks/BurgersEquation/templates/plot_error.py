import numpy as np
import matplotlib.pyplot as plt
from sys import argv

input_filename = str(argv[1])
output_filename = str(argv[2])

data = np.genfromtxt(input_filename, delimiter=',')
h = data[0]
error_short = data[1]
error_long = data[2]

x = [0.01, 0.1]

fig = plt.figure()
plt.loglog(h, error_short, 'o-', label='T = 0.3')
plt.loglog(h, error_long, 'o-', label='T = 3.0')
plt.loglog(x, x, '--', label='slope 1')
plt.xlabel('mesh-width h')
plt.ylabel('error')
plt.legend()
plt.savefig(output_filename)

print('Generated ' + output_filename)
