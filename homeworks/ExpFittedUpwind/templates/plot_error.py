from numpy import genfromtxt
from matplotlib.pyplot import figure, legend, loglog, savefig, xlabel, ylabel,show
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, skip_header=1, delimiter=',')
N = data[:,0]
err = data[:,1]

fig=figure()
loglog(N,1/N, '--', label='slope -1')
loglog(N,err, 'o-', label='L2 error')
xlabel('degrees of freedom')
ylabel('error')
legend()
savefig(output_file)
print('Generated ' + output_file)
