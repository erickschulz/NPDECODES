from matplotlib.pyplot import figure, legend, loglog, savefig, xlabel, ylabel
from numpy import exp, genfromtxt, log, polyfit
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
N = data[:,0]
error = data[:,1]

p = polyfit(log(N),log(error), 1);

figure()
loglog(N, error, 'o', label='error')
loglog(N, exp(p[0] * log(N) + p[1]), '--', label='slope: ' + str(round(p[0], 2)))
xlabel('N = dimension of FE space')
ylabel(r'$|B(u_h)-B(u)|$')
legend()
savefig(output_file)

print('Generated ' + output_file)
