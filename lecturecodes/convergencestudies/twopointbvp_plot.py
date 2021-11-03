import sys
import os
import numpy as np
import matplotlib.pyplot as plt

input_file = os.path.abspath(sys.argv[1])
output_folder = os.path.dirname(input_file)

data = np.genfromtxt(input_file, delimiter=',')
N = data[:,0]
norm_max = data[:,1]
norm_H1 = data[:,2]
norm_L2 = data[:,3]

plt.plot(N, norm_max, 'o', markersize=1, label='maximum norm')
plt.plot(N, norm_L2, 'o', markersize=1, label='$L^2$ norm')
plt.plot(N, norm_H1, 'o', markersize=1, label='$H^1$ semi-norm')
plt.legend(loc='upper right')
plt.grid()
plt.xlabel('M')
plt.ylabel('error norms')
plt.savefig(os.path.join(output_folder, '1DFD2pBVPerrlin.eps'))

plt.gcf().clear()

plt.loglog(N, norm_max, 'o', markersize=1, label='maximum norm')
plt.loglog(N, norm_L2, 'o', markersize=1, label='$L^2$ norm')
plt.loglog(N, norm_H1, 'o', markersize=1, label='$H^1$ semi-norm')
plt.legend(loc='lower left')
plt.grid()
plt.xlabel('M')
plt.ylabel('error norms')
plt.savefig(os.path.join(output_folder, '1DFD2pBVPerrlog.eps'))
