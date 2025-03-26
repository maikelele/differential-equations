import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

S_global_06 = pd.read_csv("lab04/global_relaxation_06.csv")
S_global_10 = pd.read_csv("lab04/global_relaxation_10.csv")
S_local_10 = pd.read_csv("lab04/local_relaxation_10.csv")
S_local_14 = pd.read_csv("lab04/local_relaxation_14.csv")
S_local_18 = pd.read_csv("lab04/local_relaxation_18.csv")
S_local_19 = pd.read_csv("lab04/local_relaxation_19.csv")

theta_06 = np.loadtxt("lab04/theta_06.csv", delimiter=",")
theta_10 = np.loadtxt("lab04/theta_10.csv", delimiter=",")

V_06 = np.loadtxt("lab04/potential_V_06.csv", delimiter=",")
V_10 = np.loadtxt("lab04/potential_V_10.csv", delimiter=",")

fig, axes = plt.subplots(3, 2, figsize=(14, 18))

axes[0, 0].set_xscale("log")
axes[0, 0].plot(S_global_06["Iteration"], S_global_06["S_global_06"], label="Global ω=0.6, {}".format(S_global_06["Iteration"].iloc[-1]))
axes[0, 0].plot(S_global_10["Iteration"], S_global_10["S_global_10"], label="Global ω=1.0, {}".format(S_global_10["Iteration"].iloc[-1]))
axes[0, 0].set_title("Relaksacja globalna")
axes[0, 0].set_xlabel("nr iteracji")
axes[0, 0].set_ylabel("S")
axes[0, 0].legend()
axes[0, 0].grid()

axes[0, 1].set_xscale("log")
axes[0, 1].plot(S_local_10["Iteration"], S_local_10["S_local_10"], label="Local ω=1.0, {}".format(S_local_10["Iteration"].iloc[-1]))
axes[0, 1].plot(S_local_14["Iteration"], S_local_14["S_local_14"], label="Local ω=1.4, {}".format(S_local_14["Iteration"].iloc[-1]))
axes[0, 1].plot(S_local_18["Iteration"], S_local_18["S_local_18"], label="Local ω=1.8, {}".format(S_local_18["Iteration"].iloc[-1]))
axes[0, 1].plot(S_local_19["Iteration"], S_local_19["S_local_19"], label="Local ω=1.9, {}".format(S_local_19["Iteration"].iloc[-1]))
axes[0, 1].set_title("Relaksacja lokalna")
axes[0, 1].set_xlabel("nr iteracji")
axes[0, 1].set_ylabel("S")
axes[0, 1].legend()
axes[0, 1].grid()

im1 = axes[1, 0].imshow(theta_06.transpose(), origin="lower", cmap="plasma")
axes[1, 0].set_title("Blad: Relaksacja globalna w=0.6")
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
fig.colorbar(im1, ax=axes[1, 0])

im2 = axes[1, 1].imshow(theta_10.transpose(), origin="lower", cmap="plasma")
axes[1, 1].set_title("Blad: Relaksacja globalna w=1.0")
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("y")
fig.colorbar(im2, ax=axes[1, 1])

im3 = axes[2, 0].imshow(V_06.transpose(), origin="lower", cmap="plasma")
axes[2, 0].set_title("Potencjal: Relaksacja globalna w=0.6")
axes[2, 0].set_xlabel("x")
axes[2, 0].set_ylabel("y")
fig.colorbar(im3, ax=axes[2, 0])

im4 = axes[2, 1].imshow(V_10.transpose(), origin="lower", cmap="plasma")
axes[2, 1].set_title("Potencjal: relaksacja globalna w=1.0")
axes[2, 1].set_xlabel("x")
axes[2, 1].set_ylabel("y")
fig.colorbar(im4, ax=axes[2, 1])

fig.subplots_adjust(wspace=1, hspace=90)
plt.tight_layout()
plt.show()
