from matplotlib import pyplot as plt
import numpy as np

DELTA = 0.1
NX = 150
DT = 0.05
NT = 1000

alpha = 0.00
beta_values = [0.00, 0.10, 1.00]
xf = 2.5

u_filenames = []
E_filenames = []
beta_counter = 0
for beta in beta_values:
    E_filenames.append(f"e_a_{alpha:.2f}_b_{beta:.2f}.dat")
    u_filenames.append(f"u_a_{alpha:.2f}_b_{beta:.2f}.dat")
# if (beta == 1):
#     E_filenames.append(f"e_a_{1.00:.2f}_b_{beta:.2f}_xf_{xf:.2f}.dat")
#     u_filenames.append(f"u_a_{1.00:.2f}_b_{beta:.2f}_xf_{xf:.2f}.dat")
E_filenames.append("e_a_1.00_b_1.00_xf_2.50.dat")
u_filenames.append("u_a_1.00_b_1.00_xf_2.50.dat")

E_values = []
u_values = []
for idx, (u_filename, E_filename) in enumerate(zip(u_filenames, E_filenames)):
    if (idx != 3):
        u_values.append(np.loadtxt(u_filename).reshape((NT, NX + 1)).T)
        E_values.append(np.loadtxt(E_filename))

u_values_xf = np.loadtxt(u_filenames[-1]).reshape(NT, NX + 1).T
E_values_xf = np.loadtxt(E_filenames[-1])


def plot():
    fig, axs = plt.subplots(3, 2, figsize=(18, 9))
    for idx, E_value in enumerate(E_values):
        axs[0, 0].plot(E_value, label=f"alpha: {alpha}, beta: {beta_values[idx]}")
    axs[0, 0].set_title("E(t)")
    axs[0, 0].legend()

    axs[0, 1].plot(E_values_xf, label=f"alpha: {1:.1f}, beta: {beta_values[idx]}, xf: {xf}")
    axs[0, 1].set_title("E(t)")
    axs[0, 1].legend()

    im1 = axs[1, 0].imshow(u_values[0])
    axs[1, 0].set_title(f"alpha: {alpha}, beta: {beta_values[0]}")
    fig.colorbar(im1, ax=axs[1, 0])

    im2 = axs[1, 1].imshow(u_values[1])
    axs[1, 1].set_title(f"alpha: {alpha}, beta: {beta_values[1]}")
    fig.colorbar(im2, ax=axs[1, 1])

    im3 = axs[2, 0].imshow(u_values[2])
    axs[2, 0].set_title(f"alpha: {alpha}, beta: {beta_values[2]}")
    fig.colorbar(im3, ax=axs[2, 0])

    im4 = axs[2, 1].imshow(u_values_xf)
    axs[2, 1].set_title(f"alpha: {1}, beta: {beta_values[2]}, xf: {xf}")
    fig.colorbar(im4, ax=axs[2, 1])

    plt.show()

if __name__ == "__main__":
    plot()