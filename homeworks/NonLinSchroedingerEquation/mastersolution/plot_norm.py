from matplotlib.pyplot import figure, plot, savefig, title, xlabel, ylabel, yticks
from numpy import genfromtxt, linspace
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
t = data[0]
norm = data[1]

fig = figure()
title('Norm of the Approximate Solution')
plot(t, norm, '-')
xlabel('t')
ylabel(r"$\Vert u(t)\Vert_{L^2(\mathbb{R}^2;\mathbb{C})}$")
yticks(linspace(0.0, 2.0 * max(norm), 11, endpoint=True))
savefig(output_file)

print('Generated ' + output_file)
