import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sys import argv

def main():
    
    input_file = str(argv[1])
    output_file = str(argv[2])
    
    data = pd.read_csv(input_file, sep=',', header = None)

    num_cells = data[0].tolist()
    l2_errors = data[1].tolist()

    fig = plt.figure()
    plt.plot(num_cells, l2_errors)
    plt.grid()
    plt.xlabel("Number of Cells")
    plt.xscale('log')
    plt.ylabel("L2 Error")
    plt.yscale('log')
    plt.tight_layout()

    plt.savefig(output_file)
    print('Generated ' + output_file)

if __name__ == '__main__':
    main()
