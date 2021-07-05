#!/usr/bin/env python
#-*- codin:utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpltools import annotation


input_file = os.path.abspath(sys.argv[1])
output_folder = os.path.dirname(input_file)

#read header
with open(input_file) as f:
    error_names = f.readline()[:-1].split(",")[1:]

#exctract data
data = np.genfromtxt(input_file, delimiter=',', skip_header=1)
h = data[:,0]
errors = data[:,1:]

#Plot errors
for i in range(len(error_names)):
    plt.loglog(h,errors[:,i], 'o-', markersize=5, label=error_names[i])


min_error = np.min(errors)
annotation.slope_marker((h[1],min_error), 1.0, invert=False)
annotation.slope_marker((h[3],min_error), 2.0, invert=False)

#Label plot
plt.legend()
plt.xlabel('h')
plt.ylabel('Error Norms')
plt.grid()

plt.savefig(os.path.join(output_folder, 'convergence.eps'))
