from matplotlib.pyplot import figure, legend, plot, savefig, tight_layout, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
x = data[0]
zeta0 = data[1]
zetaT = data[2]

fig = figure()
plot(x, zeta0, 'k-', label=r'$z(x, 0)$')
plot(x, zetaT, 'c--', label=r'$z(x, T)$')
xlabel(r'$x$')
ylabel(r'$z(x, t)$')
legend(framealpha=1)
tight_layout()
savefig(output_file)

print('Generated ' + output_file)
