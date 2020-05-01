from matplotlib.pyplot import figure, legend, plot, savefig, title, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
t = data[0]
norm = data[1]

fig = figure()
title('Norm of the Approximate Solution')
plot(t, norm, '-', label='norm')
xlabel('t')
ylabel(r"$\vert u(t)\vert_{L^2(\mathbb{R}^2;\mathbb{C})})$") #
legend()
savefig(output_file)

print('Generated ' + output_file)
