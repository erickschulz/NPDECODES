#!/usr/bin/env python
#-*- codin:utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

#input/output
input_folder = os.path.abspath(sys.argv[1])
input_file = os.path.join(input_folder,"convergence_potential.csv")
output_folder = os.path.join(input_folder, "plots")
os.makedirs(output_folder, exist_ok=True)

#read data
data = np.genfromtxt(input_file, delimiter=',', skip_header=1)
h = data[:,0]
error_potential = data[:,1]
h2 = np.square(h)
h2 = 2.0*error_potential[-1]/h2[-1]*h2
#Plot errors
plt.loglog(h,error_potential, 'o-', markersize=5, label="Error in u(x) (Potential Method)")
plt.loglog(h, h2,'--', label="O(h^2)" )
#Label plot
plt.legend()
plt.xlabel('h')
plt.ylabel('Error ')
plt.grid()

plt.savefig(os.path.join(output_folder, 'convergence_potential.eps'))