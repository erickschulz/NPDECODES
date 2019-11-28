from matplotlib.pyplot import figure, legend, loglog, savefig, xlabel, ylabel
from numpy import exp, genfromtxt, log, polyfit
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
M = data[0]

name = ['LaxWendroffRP', 'LaxWendroffSmoothU0', 'GodunovSmoothU0']
color = ['b', 'g', 'c']
n = len(name)

error = []
p = []
for i in range(n):
	error.append(data[i + 1]) 
	p.append(polyfit(log(M), log(error[i]), 1))

fig = figure()
for i in range(n):
	loglog(M, error[i], color[i] + 'o', label = name[i])
	loglog(M, exp(p[i][0] * log(M) + p[i][1]), color[i] + '--', label='slope = ' + str(round(p[i][0], 2)))
xlabel('number of time steps M')
ylabel('error')
legend()
savefig(output_file)

print('Generated ' + output_file)
