import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("error.csv", delimiter=',')
h = data[0]
error_short = data[1]
error_long = data[2]

x = [0.01, 0.1]

fig = plt.figure()
plt.loglog(h, error_short, 'o-', label='T = 0.3')
plt.loglog(h, error_long, 'o-', label='T = 3.0')
plt.loglog(x, x, '--', label='slope 1')
plt.xlabel('mesh-width h')
plt.ylabel('error')
plt.legend()
plt.savefig("error.png")

print("Generated error.png")