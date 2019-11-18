import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("solution.csv", delimiter=',')
x = data[0]
y = data[1]

plt.figure()
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('u(x,1)')
plt.savefig("solution.eps")

print("Generated solution.eps")
