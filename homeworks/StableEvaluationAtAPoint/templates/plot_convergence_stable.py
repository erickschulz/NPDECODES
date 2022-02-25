#!/usr/bin/env python
#-*- codin:utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt


input_folder = os.path.abspath(sys.argv[1])
input_file = os.path.join(input_folder,"convergence_stable.csv")
output_folder = os.path.join(input_folder, "plots")
os.makedirs(output_folder, exist_ok=True)

#read data
data = np.genfromtxt(input_file, delimiter=',', skip_header=1)
h = data[:,0]
error_direct = data[:,1]
error_stable = data[:,2]
h2 = np.square(h)
h2_stable = 2.0*error_stable[1]/h2[1]*h2
h2_direct = 3.0*error_direct[-1]/h2[-1]*h2
#Plot errors

plt.loglog(h,error_direct, 'o-', markersize=5, label="Error in u(x) (Direct)")
plt.loglog(h,error_stable, 'o-', markersize=5, label="Error in u(x) (Stable)")
plt.loglog(h,h2_stable, '--', label="O(h^2)")
plt.loglog(h,h2_direct, '--', label="O(h^2)")

#Label plot
plt.legend()
plt.xlabel('h')
plt.ylabel('Error')
plt.grid()

plt.savefig(os.path.join(output_folder, 'convergence_stable.eps'))