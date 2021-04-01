# Small python plotting script for CSV data
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

# Read command line arguments
if len(argv) < 3:
    print("Usage: " + str(argv[0]) + " <infile> <outputdir>")
input_file = str(argv[1])
output_dir = str(argv[2])
print("Reading data from ", input_file, ", writing plot to ", output_dir)


# Read CSV data using dedicated function
data = np.genfromtxt(input_file, delimiter=',')
t = data[:,0]
Y = data[:,1]
Y_exact = np.tan(t)

#Plot exact vs. approximate solution
plt.figure()
plt.plot(t, Y, "-+", label="Approximate Solution")
plt.plot(t,Y_exact, label="Exact Solution")
plt.legend()
plt.xlabel("t")
plt.ylabel("y")
plt.title("Approximate and Exact Solution")
plt.savefig(output_dir + "/tangent.eps")


#Plot showing evolution of time against number of
plt.figure()
plt.plot(t)
plt.xlabel("idx")
plt.ylabel("t")
plt.title("Time vs. Index of Timestep ")
plt.savefig(output_dir + "/time.eps")

#Plot showing timestep size against time.
plt.figure()
plt.plot(t[0:-1],t[1:] - t[0:-1])
plt.xlabel("t")
plt.ylabel("h")
plt.title("Timestep size vs. Time")
plt.savefig(output_dir + "/timesteps.eps")


