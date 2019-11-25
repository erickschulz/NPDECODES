from matplotlib.pyplot import figure, legend, plot, savefig, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
x = data[0]
solution_short = data[1]
solution_long = data[2]

fig = figure()
plot(x, solution_short, '-', label='T = 0.3')
plot(x, solution_long, '--', label='T = 3.0')
xlabel('x')
ylabel('u(x, T)')
legend()
savefig(output_file)

print('Generated ' + output_file)
