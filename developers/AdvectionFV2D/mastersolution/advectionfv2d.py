import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def main():
    data = pd.read_csv('advectionfv2d.csv', sep=',', header = None)

    num_cells = data[0].tolist()
    hmin_inv = data[1].tolist()
    l1_errors = data[2].tolist()
    l2_errors = data[3].tolist()
    cfl_thres = data[4].tolist()

    print("\nRates: L2Error - Number Cells")
    for i in range(0, len(num_cells) - 1):
        rate = (np.log(l2_errors[i]) - np.log(l2_errors[i + 1]))/(np.log(num_cells[i + 1]) - np.log(num_cells[i]))
        print(rate)

    print("\nRates: L2Error - Hmin_inv")
    for i in range(0, len(num_cells) - 1):
        rate = (np.log(l2_errors[i]) - np.log(l2_errors[i + 1]))/(np.log(hmin_inv[i + 1]) - np.log(hmin_inv[i]))
        print(rate)

    print("\nRates: L1Error - Number Cells")
    for i in range(0, len(num_cells) - 1):
        rate = (np.log(l1_errors[i]) - np.log(l1_errors[i + 1]))/(np.log(num_cells[i + 1]) - np.log(num_cells[i]))
        print(rate)

    print("\nRates: L1Error - Hmin_inv")
    for i in range(0, len(num_cells) - 1):
        rate = (np.log(l1_errors[i]) - np.log(l1_errors[i + 1]))/(np.log(hmin_inv[i + 1]) - np.log(hmin_inv[i]))
        print(rate)

    fig, axs = plt.subplots(3, 2)

    axs[0, 0].set_title("L2 Errors - Number Cells")
    axs[0, 0].grid()
    axs[0, 0].plot(num_cells, l2_errors)
    axs[0, 0].set_xlabel("Number of Cells")
    axs[0, 0].set_xscale('log')
    axs[0, 0].set_ylabel("L2 Error")
    axs[0, 0].set_yscale('log')

    axs[0, 1].set_title("L2 Errors - 1/hmin")
    axs[0, 1].grid()
    axs[0, 1].plot(hmin_inv, l2_errors)
    axs[0, 1].set_xlabel("1/hmin")
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_ylabel("L2 Error")
    axs[0, 1].set_yscale('log')

    axs[1, 0].set_title("L1 Errors - Number Cells")
    axs[1, 0].grid()
    axs[1, 0].plot(num_cells, l1_errors)
    axs[1, 0].set_xlabel("Number of Cells")
    axs[1, 0].set_xscale('log')
    axs[1, 0].set_ylabel("L1 Error")
    axs[1, 0].set_yscale('log')

    axs[1, 1].set_title("L1 Errors - 1/hmin")
    axs[1, 1].grid()
    axs[1, 1].plot(hmin_inv, l1_errors)
    axs[1, 1].set_xlabel("1/hmin")
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_ylabel("L1 Error")
    axs[1, 1].set_yscale('log')

    axs[2, 0].set_title("Threshold vs. 1/hmin - Number Cells")
    axs[2, 0].grid()
    axs[2, 0].plot(num_cells, cfl_thres)
    axs[2, 0].plot(num_cells, hmin_inv)
    axs[2, 0].set_xlabel("Number Cells")
    axs[2, 0].set_xscale('log')
    axs[2, 0].set_ylabel("L1 Error")
    axs[2, 0].set_yscale('linear')
    axs[2, 0].legend(['Threshold', '1 / hmin'], loc='upper left')
    fig.tight_layout()

    axs[2, 1].set_visible(False)

    plt.savefig("advectionfv2d.png")
    plt.show()
    
    return 0

if __name__ == '__main__':
    main()
