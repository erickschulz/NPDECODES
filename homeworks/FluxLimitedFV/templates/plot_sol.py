from matplotlib.pyplot import figure, plot, savefig, title, xlabel, ylabel
from numpy import genfromtxt, linspace
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
x = data[0]
sol = data[1]

fig = figure()
title('Solution')
plot(x, sol, '-')
xlabel('x')
ylabel('sol')
savefig(output_file)

print('Generated ' + output_file)
