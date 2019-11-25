from numpy import genfromtxt
from matplotlib.pyplot import figure, legend, loglog, savefig, xlabel, ylabel
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
h = data[0]
error_short = data[1]
error_long = data[2]

x = [0.01, 0.1]

fig = figure()
loglog(h, error_short, 'o-', label='T = 0.3')
loglog(h, error_long, 'o-', label='T = 3.0')
loglog(x, x, '--', label='slope 1')
xlabel('mesh-width h')
ylabel('error')
legend()
savefig(output_file)

print('Generated ' + output_file)
