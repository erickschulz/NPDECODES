#!/usr/bin/env python
#-*- codin:utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpltools import annotation


input_folder = os.path.abspath(sys.argv[1])
input_file = os.path.join(input_folder,"results_errors.txt")
output_folder = os.path.join(input_folder, "plots")
os.makedirs(output_folder, exist_ok=True)


#read header
with open(input_file) as f:
    error_names = f.readline()[:-1].split(",")[1:]

#read data
data = np.genfromtxt(input_file, delimiter=',', skip_header=1)
h = data[:,0]
errors = data[:,1:]

#Plot errors
for i in range(len(error_names)):
    plt.loglog(h,errors[:,i], 'o-', markersize=5, label=error_names[i])

#Add slope triangles
min_error = np.min(errors)
annotation.slope_marker((h[1],min_error), 1.0, invert=False)
annotation.slope_marker((h[3],min_error), 2.0, invert=False)

#Label plot
plt.legend()
plt.xlabel('h')
plt.ylabel('Error Norms')
plt.grid()

plt.savefig(os.path.join(output_folder, 'convergence.eps'))
