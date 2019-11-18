import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("ufinal.csv", delimiter=',')
x = data[0]
ufinal = data[1]

plt.figure()
plt.plot(x, ufinal, '-')
plt.xlabel('x')
plt.ylabel('u(x, T)')
plt.savefig("ufinal.eps")

print("Saved plot as ufinal.eps")