import numpy as np
import matplotlib.pyplot as plt

def load_potential(filename):
    data = np.loadtxt(filename)
    x, y, V = data[:, 0], data[:, 1], data[:, 2]
    grid_size = int(max(x)) + 1
    potential = np.zeros((grid_size, grid_size))
    for i in range(len(V)):
        potential[int(x[i]), int(y[i])] = V[i]
    return potential

fig, axs = plt.subplots(2, 3, figsize=(12, 8))
for idx, (grid_size, file, scale) in enumerate([(50, "potential_50.dat", 10), 
                                                 (100, "potential_100.dat", 10), 
                                                 (200, "potential_200.dat", 10), 
                                                 (100, "potential_eps1.dat", 0.8), 
                                                 (100, "potential_eps2.dat", 0.8), 
                                                 (100, "potential_eps10.dat", 0.8)]):
    potential = load_potential(file)
    im = axs[idx // 3, idx % 3].imshow(potential.T, cmap="seismic", origin="lower", vmin=-scale, vmax=scale)
    
    fig.colorbar(im, ax=axs[idx // 3, idx % 3], label="Potential (V)")
    
    axs[idx // 3, idx % 3].set_title(f"Potential Map (nx = ny = {grid_size})")
    axs[idx // 3, idx % 3].set_xlabel("x")
    axs[idx // 3, idx % 3].set_ylabel("y")

plt.tight_layout()
plt.show()
