import sys
import matplotlib.pyplot as plt
import numpy as np

# define a function reading the data from file
# File format: arbitrary number of pairs of lines
# 1st line: <key> = <float value> pairs, comma separated
# 2nd line: comma separated float values
# Here special keys: 'a','b','T'
def readdat(fname):
    params = {} # dictionary for parameters 
    data = []   # list of arrays for data
    # Open the data file
    with open(fname) as f:
        # Read next line of data file into a string variable
        for line in f.readlines():
            # Split string into comma-separated substrings and store them in a list
            tmp = line.split(sep=",")
            # Parameter line: contains '='
            if "=" in line: # parameter line
                tmp2 = [a.split(sep="=") for a in tmp]
                params.update({k.strip() : v.strip() for k, v in tmp2})
            else: # data line
                try:
                    # Check whether data are numbers 
                    float(tmp[0])
                    # Append list of numbers to list of data arrays 
                    data += [[float(v) for v in tmp]]
                except:
                    pass
    return params, data

def plotfvsol(filename):
    params, datas = readdat(filename)
    print("Parameters: ",params)
    # Extract relevant parameters 
    try:
        a = float(params['a'])
        b = float(params['b'])
        T = float(params['T'])
    except:
        print("Missing parameters a, b, T!")
    # Ensure that [a,b] is an interval 
    if a>b:
        a, b = b, a
    # Plot data
    fig, ax = plt.subplots()  # Create a figure containing a single axis
    for i, data in enumerate(datas):
        print("|data[", i, "]| = ",len(data))
        # Number of cell values 
        N = len(data)
        h = (b-a)/N
        x = np.linspace(a+h/2,b-h/2,N)
        plt.plot(x,data,label=str('N={:d}'.format(N)),linewidth=1)
    plt.title(filename + ': solution at T = ' + str(T));    
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.legend()
    plt.show()
    # Save figure
    outfnname = filename.split(sep='.')
    plt.savefig(outfnname[0] + ".eps")
    print("Figure saved in" + outfnname[0] + ".eps")
    plt.close()

if __name__ == "__main__":
    filename = sys.argv[1]
    print ("Reading data from ", filename)
    plotfvsol(filename)
        
