import matplotlib.pyplot as plt
import numpy as np


def load_data(filename, cols):
    with open(filename) as source:
        lines = source.readlines()
        return [np.array([float(line.split()[i]) for line in lines]) for i in range(cols)]


def reshape_data(arrays, shape):
    return [array.reshape(shape) for array in arrays]


def configure_axis(axis, title, xlabel, ylabel):
    axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)


def plot_contour(ax, X, Y, Z, title, levels, cmap, vmin=None, vmax=None):
    configure_axis(ax, title, "x", "y")
    heatmap = ax.contour(X, Y, Z, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0.5)
    return heatmap


def plot_pcolormesh(ax, X, Y, Z, title, cmap):
    configure_axis(ax, title, "x", "y")
    heatmap = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap)
    return heatmap


def main():
    fig, axs = plt.subplots(5, 2, figsize=(22, 13))

    # Load and reshape data
    x1, y1, psi1, zeta1 = reshape_data(load_data("lab07/psi_zeta_Q-1000.txt", 4), (199, 89))
    x2, y2, psi2, zeta2 = reshape_data(load_data("lab07/psi_zeta_Q-4000.txt", 4), (199, 89))
    x3, y3, psi3 = reshape_data(load_data("lab07/psi_zeta_Q4000.txt", 3), (199, 89))
    x4, y4, u4, v4 = reshape_data(load_data("lab07/u_v_Q-1000.txt", 4), (199, 89))
    x5, y5, u5, v5 = reshape_data(load_data("lab07/u_v_Q-4000.txt", 4), (199, 89))

    # Plot data
    plots = [
        (axs[0][0], x1, y1, psi1, "Funkcja strumienia Psi(x,y) Q=-1000", 40, 'jet', -55, -50),
        (axs[0][1], x1, y1, zeta1, "Funkcja wirowosci Zeta(x,y) Q=-1000", 70, 'jet'),
        (axs[1][0], x4, y4, u4, "Pozioma skladowa predkosci u(x,y) Q=-1000", None, 'jet'),
        (axs[1][1], x4, y4, v4, "Pionowa skladowa predkosci v(x,y) Q=-1000", None, 'jet'),
        (axs[2][0], x2, y2, psi2, "Funkcja strumienia Psi(x,y) Q=-4000", 30, 'jet', -220, -200),
        (axs[2][1], x2, y2, zeta2, "Funkcja wirowosci Zeta(x,y) Q=-4000", 70, 'jet'),
        (axs[3][0], x5, y5, u5, "Pozioma skladowa predkosci u(x,y) Q=-4000", None, 'jet'),
        (axs[3][1], x5, y5, v5, "Pionowa skladowa predkosci v(x,y) Q=-4000", None, 'jet'),
        (axs[4][0], x3, y3, psi3, "Funkcja strumienia Psi(x,y) Q=4000", 30, 'jet', 200, 230),
    ]

    for ax, X, Y, Z, title, levels, cmap, *vlims in plots:
        if levels:
            heatmap = plot_contour(ax, X, Y, Z, title, levels, cmap, *vlims)
        else:
            heatmap = plot_pcolormesh(ax, X, Y, Z, title, cmap)
        plt.colorbar(heatmap, ax=ax)

    # Adjust layout and show
    fig.delaxes(axs[4][1])
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.show()


if __name__ == "__main__":
    main()
