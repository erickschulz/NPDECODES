from numpy import genfromtxt
from matplotlib.pyplot import figure, legend, loglog, savefig, xlabel, ylabel
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
N = data[:,0]

fig = figure()
for i in range(1, len(data[0])):
  loglog(N, data[:,i], 'o-', label='Assembler ' + str(i))
xlabel('Num dofs')
ylabel('error')
legend()
savefig(output_file)

print('Generated ' + output_file)
