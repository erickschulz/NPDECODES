import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("convergence.csv", delimiter=',')
M = data[0]

name = ['LaxWendroffRP', 'LaxWendroffSmoothU0', 'GodunovSmoothU0']
color = ['b', 'g', 'c']
n = len(name)

error = []
p = []
for i in np.arange(n):
	error.append(data[i + 1]) 
	p.append(np.polyfit(np.log(M), np.log(error[i]), 1))

fig = plt.figure()
for i in np.arange(n):
	plt.loglog(M, error[i], color[i] + 'o', label = name[i])
	plt.loglog(M, np.exp(p[i][0] * np.log(M) + p[i][1]), color[i] + '--', label='slope = ' + str(np.round(p[i][0], 2)))
plt.xlabel('number of time steps M')
plt.ylabel('error')
plt.legend()
plt.savefig("convergence.eps")

print("Generated convergence.eps")
