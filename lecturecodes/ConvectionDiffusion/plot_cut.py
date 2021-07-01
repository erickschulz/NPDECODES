import sys
import os
import numpy as np
import matplotlib.pyplot as plt

input_file = os.path.abspath(sys.argv[1])
output_folder = os.path.dirname(input_file)

data = np.genfromtxt(input_file, delimiter=',')
T = data[:,0]
Values = data[:,1]
plt.plot(T,Values)
plt.show()