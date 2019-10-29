import os
import sys
import numpy as np
import matplotlib.pyplot as plt

input_file = sys.argv[1]
output_file = os.path.splitext(input_file)[0] + '.eps'

data = np.genfromtxt(input_file, delimiter=',')
x = data[0]
ufinal = data[1]

plt.figure()
plt.plot(x, ufinal, '-')
plt.xlabel('x')
plt.ylabel('u(x, T)')
plt.savefig(output_file)

print("Generated " + output_file)
