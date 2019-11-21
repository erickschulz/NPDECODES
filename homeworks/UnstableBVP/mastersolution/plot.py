# Python plotting script for NumPDE homework problem "UnstableBVP"
# Call without any arguments
# Generates EPS graphics files for problem
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("h1.txt")
levels = data[:, 0]
h1 = data[:, 1:]

print(data)

# labels = ["Above", "Intersecting", "Below"]
labels = ["curve 1", "curve 2", "curve 3"]
colors = ["b", "r", "k"]

# Plot H1 seminorms 
plt.figure()
for i, (color, label) in enumerate(zip(colors, labels)):
    plt.plot(levels, h1[:, 2*i], color=color, marker="o", label=label)

plt.xlabel("Refinement level $k$")
plt.ylabel("$H_1$ seminorm of the solution $|u_k|_{H_1}$")
plt.title("$H_1$ seminorm for different domains")
plt.legend(loc="best")
plt.savefig("h1_seminorms.eps",bbox_inches="tight")


# Absolute differences of H1 seminorms, logscale in y direction
plt.figure()
for i, (color, label) in enumerate(zip(colors, labels)):
    # exclude last value since obviously it is zero: u_L - u_L
    # add a little something to the y-values to avoid 0 values
    plt.semilogy(levels[:-1], h1[:-1, 2*i + 1] + 1e-16, color=color, marker="o", label=label)

plt.xlabel("Refinement level $k$")
plt.ylabel("$||u_L|_{H_1} - |u_k|_{H_1}|$")
plt.title("Absolute difference of $H_1$ seminorms for different domains")
plt.legend(loc="best")
plt.savefig("h1_differences.eps",bbox_inches="tight")

