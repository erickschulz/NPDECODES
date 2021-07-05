import sys
import os
import numpy as np
import matplotlib.pyplot as plt

input_folder = os.path.abspath(sys.argv[1])
output_folder = input_folder

input_filenames = ["results_upwind.txt", "results_supg.txt", "results_standard_FEM.txt"]
output_filenames = ["upwindcut.eps", "supgcut.eps", "standardfemcut.eps"]
titles = ["Upwind Quadrature Solution", "Streamline-Diffusion (SUPG) Solution", "Standard FEM Solution"]
for i in range(len(input_filenames)):
    file = os.path.join(input_folder,input_filenames[i])
    data = np.genfromtxt(file,delimiter=',')
    T = data[:,0]
    values = data[:,1]
    plt.plot(T,values)
    plt.xlabel("$t$")
    plt.ylabel("$u(\gamma(t))$")
    plt.suptitle(titles[i])
    plt.title("along $\gamma(t) = (t,1-t)$")

    plt.savefig(os.path.join(output_folder,output_filenames[i]))
    plt.gcf().clear()

