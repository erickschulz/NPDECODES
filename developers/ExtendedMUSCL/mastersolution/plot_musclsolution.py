# Python script for visualization of FV solution computed with MUSCL scheme
import sys
import matplotlib.pyplot as plt
import numpy as np

# Read filenames
if len(sys.argv) != 3:
    print ("Usage: " + sys.argv[0] + " <infile> <outfile>")
else:
    input_file = str(sys.argv[1])
    output_file = str(sys.argv[2])
    # data 1-periodic on an equidistant mesh 
    vals = np.genfromtxt(input_file, delimiter=',')
    # number of cells
    n = vals.size
    # cell spacing
    h = 1.0/n
    # endpoints and cell centers
    x = np.array([0.0]);
    x = np.append(x,np.linspace(0.5*h,1-0.5*h,n));
    x = np.append(x,np.array([1.0]))
    # Add endpoint values
    y = np.insert(vals,0,0.5*(vals[0]+vals[n-1]))
    y = np.append(y,0.5*(vals[0]+vals[n-1]))
    # Create plot
    plt.figure()
    plt.plot(x,y,'r.-',linewidth=0.5)
    plt.xlabel('x')
    plt.ylabel('u')
    filename = (input_file.split('/'))[-1]
    plt.title(filename + ", n = " + str(n));
    plt.savefig(output_file)
    print('Generated ' + output_file)
    
    
