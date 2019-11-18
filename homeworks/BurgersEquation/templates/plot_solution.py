import numpy as np
import matplotlib.pyplot as plt
from sys import argv

input_filename = str(argv[1])
output_filename = str(argv[2])

data = np.genfromtxt(input_filename, delimiter=',')
x = data[0]
solution_short = data[1]
solution_long = data[2]

fig = plt.figure()
plt.plot(x, solution_short, '-', label='T = 0.3')
plt.plot(x, solution_long, '--', label='T = 3.0')
plt.xlabel('x')
plt.ylabel('u(x, T)')
plt.legend()
plt.savefig(output_filename)

print('Generated ' + output_filename)
