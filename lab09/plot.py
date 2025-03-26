from matplotlib import pyplot as plt
import numpy as np

NX = NY = 40

iter_values = [100, 200, 500, 1000, 2000]
filenames_t = ['t_' + str(iter) + '.dat' for iter in iter_values]
filenames_l = ['l_' + str(iter) + '.dat' for iter in iter_values]

fig, axs = plt.subplots(2, 5, figsize=(12, 8))
for idx, file in enumerate(filenames_t):
    array = np.loadtxt(file).reshape(NX + 1, NY + 1)
    axs[idx // 5, idx % 5].imshow(array.T)
for idx, file in enumerate(filenames_l):
    array = np.loadtxt(file).reshape(NX - 1, NY - 1)
    axs[idx // 5 + 1, idx % 5].imshow(array.T)

plt.show()