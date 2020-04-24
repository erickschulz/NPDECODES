from matplotlib.pyplot import figure, plot, savefig, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
x = data[0]
ufinal = data[1]

figure()
plot(x, ufinal, '-')
xlabel('x')
ylabel('u(x, T)')
savefig(output_file)

print('Generated ' + output_file)
