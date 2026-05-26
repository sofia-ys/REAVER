# from data import get_merged_data
import matplotlib.pyplot as plt
import numpy as np

def plot_raan_hist(data, show=False):
    raans = data["RA_OF_ASC_NODE"].values
    plt.hist(raans, bins=18)

    if show:
        plt.show()

def plot_inclination_hist(data, show=False):
    raans = data["INCLINATION"].values
    plt.hist(raans, bins=18)

    if show:
        plt.show()

# if __name__ == '__main__':
#     data = get_merged_data()
#     plot_raan_hist(data)
#     plot_inclination_hist(data)
#     plt.show()
