#!/usr/bin/env python
#-*- codin:utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

input_folder = os.path.abspath(sys.argv[1])
output_folder = os.path.join(input_folder, "plots")
os.makedirs(output_folder, exist_ok=True)

input_filenames = ["results_upwind.txt", "results_supg.txt", "results_standard_FEM.txt"]
output_filenames = ["upwind.eps", "supg.eps", "standardfem.eps"]
titles = ["Upwind Quadrature Solution", "Streamline-Diffusion (SUPG) Solution", "Standard FEM Solution"]


for i in range(len(input_filenames)):
    file = os.path.join(input_folder,input_filenames[i])
    data = np.genfromtxt(file,delimiter=',')
    T = data[:,0]
    values = data[:,3]
    plt.plot(T,values, linewidth=0.8)
    plt.xlabel("$t$")
    plt.ylabel("$u(\gamma(t))$")
    plt.title(titles[i])
    plt.savefig(os.path.join(output_folder,output_filenames[i]), bbox_inches="tight")
    plt.gcf().clear()

