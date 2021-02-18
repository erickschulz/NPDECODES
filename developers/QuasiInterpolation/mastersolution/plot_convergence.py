from numpy import genfromtxt
from matplotlib.pyplot import figure, legend, loglog, savefig, title, xlabel, ylabel
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])
title_string = str(argv[3])

data = genfromtxt(input_file, delimiter=',')
h = data[0]
l2_error = data[1]
h1_error = data[2]

y2 = [hi**2 for hi in h]
y1 = [hi for hi in h]

fig = figure()
title(r"{}".format(title_string))
loglog(h, l2_error, 'go-', label=r'$L^2(\Omega)$-error', basex=2.0, basey=2.0)
loglog(h, y2, 'g--', label=r'$\mathcal{O}(h^2)$', basex=2.0, basey=2.0)
loglog(h, h1_error, 'co-', label=r'$H^1(\Omega)$-error', basex=2.0, basey=2.0)
loglog(h, y1, 'c--', label=r'$\mathcal{O}(h^1)$', basex=2.0, basey=2.0)
xlabel(r'meshwidth $h$')
ylabel(r'error $||u-Q_h(u)||$')
legend(framealpha=1.0) # avoid warning because .eps cannot store transparency
savefig(output_file)

print('Generated ' + output_file)
