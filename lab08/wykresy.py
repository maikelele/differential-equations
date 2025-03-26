import re
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
from PIL import Image

plt.figure()
for file in ('c0.dat', 'c1.dat'):
    data = np.loadtxt(file)
    cx = data
    plt.plot(cx, label=file.split('.')[0])
    plt.suptitle("Calka gestosci")
    plt.xlabel("t")
    plt.ylabel("x")
plt.show()

plt.figure()
for file in ('x0.dat', 'x1.dat'):
    data = np.loadtxt(file)
    xm = data
    plt.plot(xm, label=file.split('.')[0])
    plt.suptitle("Srednie polozenie")
    plt.xlabel("t")
    plt.ylabel("c")
plt.show()

fig, axs = plt.subplots(1, 2, figsize=(30,15))
for i, file in enumerate(('vx.dat', 'vy.dat')):
    V = np.loadtxt(file).reshape(401, 91).T
    if i == 0:
        axs[i].set_title("vx")
    if i == 1:
        axs[i].set_title("vy")
    axs[i].set_xlabel("x")
    axs[i].set_ylabel("y")
    heatmap= axs[i].pcolormesh(V, shading='auto', cmap='jet')
    
    plt.colorbar(heatmap)

plt.show()

output_u0_gif = "output_u0.gif"

image_files_u0 = []

for i,file in enumerate(sorted(glob.glob("u0*.dat"), key=os.path.getctime)):
    U = np.loadtxt(file).reshape((401, 91)).T  
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(U, shading='auto', cmap='jet', vmin=0, vmax=20)
    plt.colorbar(label="u(x,y)")
    plt.title(f"Frame {i}: {file}")
    plt.xlabel("x")
    plt.ylabel("y")
    
    image_file = f"u0_frame_{i}.png"
    plt.savefig(image_file)
    image_files_u0.append(image_file)
    plt.close()

frames = [Image.open(img) for img in image_files_u0]
frames[0].save(output_u0_gif, save_all=True, append_images=frames[1:], duration=200, loop=0)


output_u1_gif = "output_u1.gif"

image_files_u1 = []

for i,file in enumerate(sorted(glob.glob("u1*.dat"), key=os.path.getctime)):
    U = np.loadtxt(file).reshape((401, 91)).T  
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(U, shading='auto', cmap='jet')
    plt.colorbar(label="u(x,y)")
    plt.title(f"Frame {i}: {file}")
    plt.xlabel("x")
    plt.ylabel("y")
    
    image_file = f"u1_frame_{i}.png"
    plt.savefig(image_file)
    image_files_u1.append(image_file)
    plt.close()

frames = [Image.open(img) for img in image_files_u1]
frames[0].save(output_u1_gif, save_all=True, append_images=frames[1:], duration=200, loop=0)

for img_u0, img_u1 in zip(image_files_u0, image_files_u1):
    os.remove(img_u0)
    os.remove(img_u1)

plt.plot()
plt.grid(False)
plt.axis('Off')
information = "Gify znajdują się w plikach output_u0.gif i output_u1.gif"
plt.title(information)
plt.show()