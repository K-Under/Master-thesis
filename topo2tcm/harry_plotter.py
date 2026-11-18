import matplotlib.pyplot as plt
import numpy as np

"""
This script has some basic plotting functions for quickly checking the data 
"""


def showim(data):
    fig, ax = plt.subplots()
    image = ax.imshow(data[0], cmap="grey", vmin=data.min(), vmax=data.max())

    for z in range(data.shape[0]):
        image.set_data(data[z])
        ax.set_title(f"Z = {z}")
        plt.pause(0.1)

    plt.show()


def plot_data_simple(data):  # simple plot
    x_axis = np.linspace(0, len(data), len(data))
    plt.plot(x_axis, data)
    plt.grid()
    plt.show()
