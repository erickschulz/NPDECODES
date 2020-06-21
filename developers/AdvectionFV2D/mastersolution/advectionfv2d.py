import pandas as pd
import matplotlib.pyplot as plt

def main():
    data = pd.read_csv('advectionfv2d.csv', sep=',', header = None)
    
    num_cells = data[0].tolist()
    l2_errors = data[1].tolist()
    clf_thres = data[2].tolist()
    
    #L2 Error Plot
    fig1 = plt.figure()
    plt.title("L2 Errors")
    plt.plot(num_cells, l2_errors)
    plt.xlabel("Number of Cells")
    plt.ylabel("L2 Error")
    plt.yscale('log')
    plt.savefig("advectionfv2d_l2error.png")
    plt.show()
    
    #Threshold Plot
    fig2 = plt.figure()
    plt.title("CFL Threshold")
    plt.plot(num_cells, clf_thres)
    plt.xlabel("Number of Cells")
    plt.ylabel("Threshold")
    plt.savefig("advectionfv2d_cfl_threshold.png")
    
    plt.show()
    
    return 0

if __name__ == '__main__':
    main()
