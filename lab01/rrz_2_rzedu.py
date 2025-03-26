#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from matplotlib import pyplot as plt

dt = 1e-4
R = 100
L = 0.1
C = 0.001
omega0 = 1 / (np.sqrt(L*C))
T0 = 2 * np.pi / omega0
t_end = 4 * T0
Q0 = 0
I0 = 0
omegav_values = [0.5 * omega0, 0.8 * omega0, 1.0 * omega0, 1.2 * omega0]

def V(omegav, t):
    return 10 * np.sin(omegav * t)

def metoda_rk4(dt, t_end, omegav):
    t_values = np.arange(0, t_end + dt, dt)
    Q_values = [Q0]
    I_values = [I0]
    
    for i in range(1, len(t_values)):
        Q_k1 = I_values[i-1]
        I_k1 = V(omegav, t_values[i-1])/L - Q_values[i-1] / (L * C) - (R/L) * I_values[i-1]
        
        Q_k2 = I_values[i-1] + (dt/2) * I_k1
        I_k2 = V(omegav, t_values[i-1] + dt/2) / L - (1/ (L*C)) * (Q_values[i-1] + (dt/2) * Q_k1) - R/L * (I_values[i-1] + (dt/2) * I_k1)

        Q_k3 = I_values[i-1] + dt/2 * I_k2
        I_k3 = V(omegav, t_values[i-1] + dt/2) / L - (1/ (L*C)) * (Q_values[i-1] + (dt/2) * Q_k2) - R/L * (I_values[i-1] + (dt/2) * I_k2)

        Q_k4 = I_values[i-1] + dt * I_k3
        I_k4 = V(omegav, t_values[i]) / L - (1/ (L*C)) * (Q_values[i-1] + dt * Q_k3) - R/L * (I_values[i-1] + dt * I_k3)
        
        Q_values.append(Q_values[i-1] + dt/6 * (Q_k1 + 2 * Q_k2 + 2* Q_k3 + Q_k4))
        I_values.append(I_values[i-1] + dt/6 * (I_k1 + 2 * I_k2 + 2 * I_k3 + I_k4))

    return t_values, Q_values, I_values

fig, axs = plt.subplots(2, 1, figsize=(12, 8))

for omegav in omegav_values:
    t_values, Q_values, I_values = metoda_rk4(dt, t_end, omegav)
    axs[0].plot(t_values, I_values, label=f'I(t) (omegav={omegav})')
    axs[1].plot(t_values, Q_values, label=f'Q(t) (omegav={omegav})')

axs[0].legend()
axs[0].grid()
axs[0].set_xlabel("Czas (t)")
axs[0].set_ylabel("I(T)")

axs[1].legend()
axs[1].grid()
axs[1].set_xlabel("Czas (t)")
axs[1].set_ylabel("Q(T)")

plt.show()


