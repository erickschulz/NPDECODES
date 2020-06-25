# Small python plotting script for CSV data
from matplotlib.pyplot import figure, plot, savefig, xlabel, ylabel, loglog, title
from numpy import genfromtxt
from sys import argv
# Read command line arguments
if len(argv) < 3:
    print("Usage: " + str(argv[0]) + " <infile> <outfile>")
input_file = str(argv[1])
output_file = str(argv[2])
print("Reading data from ", input_file, ", writing plot to ", output_file)
# Read CSV data using dedicated function
data = genfromtxt(input_file, delimiter=',')
h = data[0] # Here: meshwidths
e = data[1] # Here: errors 

figure()
loglog(h, e)
title("Data from " + input_file)
xlabel(r'mesh_size $h$')
ylabel(r'force error $e$')
savefig(output_file)

print('Generated ' + output_file)
