#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt

y0 = 1
lam = -1
t_end = 5
time_steps = [0.01, 0.1, 1.0]

def rozwiazanie_analityczne(t, lam):
    return np.exp(lam * t)

def metoda_eulera(y0, lam, dt, t_end):
    t_values = np.arange(0, t_end + dt, dt)
    y_values = [y0]

    for i in range(1, len(t_values)):
        y_values.append(y_values[i - 1] + dt * lam * y_values[i - 1])
    return t_values, y_values

fig, axs = plt.subplots(1, 2, figsize=(14, 6))
t_analytical = np.linspace(0, t_end, 500)
y_analytical = rozwiazanie_analityczne(t_analytical, lam)
axs[0].plot(t_analytical, y_analytical, 'k--', label='Rozwiazanie analityczne')

for dt in time_steps:
    t_values, y_values = metoda_eulera(y0, lam, dt, t_end)
    axs[0].plot(t_values, y_values, label=f'Rozwiązanie numeryczne (Δt={dt})')

axs[0].set_xlabel('Czas (t)')
axs[0].set_ylabel('y(t)')
axs[0].set_title('Rozwiązania numeryczne i analityczne')
axs[0].legend()
axs[0].grid(True)

for dt in time_steps:
    t_values, y_values = metoda_eulera(y0, lam, dt, t_end)
    y_exact = rozwiazanie_analityczne(t_values, lam)
    error = y_values - y_exact
    axs[1].plot(t_values, error, label=f'Błąd (Δt={dt})')

axs[1].set_xlabel('Czas (t)')
axs[1].set_ylabel('Błąd δ(t)')
axs[1].set_title('Zmiany błędu globalnego')
axs[1].legend()
axs[1].grid(True)

# plt.tight_layout()
# plt.show()


def metoda_rk2(y0, lam, dt, t_end):
    t_values = np.arange(0, t_end + dt, dt)
    y_values = [y0]

    for i in range(1, len(t_values)):
        k1 = lam * y_values[i-1]
        k2 = lam * (y_values[i-1] + dt * k1) 
        y_values.append(y_values[i-1] + 0.5 * dt * (k1 + k2))

    return t_values, y_values

fig, axs = plt.subplots(1, 2, figsize=(14, 6))
t_analytical = np.linspace(0, t_end, 500)
y_analytical = rozwiazanie_analityczne(t_analytical, lam)
axs[0].plot(t_analytical, y_analytical, 'k--', label='Rozwiazanie analityczne')

for dt in time_steps:
    t_values, y_values = metoda_rk2(y0, lam, dt, t_end)
    axs[0].plot(t_values, y_values, label=f'Rozwiązanie numeryczne (Δt={dt})', linewidth=2)

axs[0].set_xlabel('Czas (t)')
axs[0].set_ylabel('y(t)')
axs[0].set_title('Rozwiązania numeryczne i analityczne')
axs[0].legend()
axs[0].grid(True)

for dt in time_steps:
    t_values, y_values = metoda_rk2(y0, lam, dt, t_end)
    y_exact = rozwiazanie_analityczne(t_values, lam)
    error = y_values - y_exact
    axs[1].plot(t_values, error, label=f'Błąd (Δt={dt})')

axs[1].set_xlabel('Czas (t)')
axs[1].set_ylabel('Błąd δ(t)')
axs[1].set_title('Zmiany błędu globalnego')
axs[1].legend()
axs[1].grid(True)

# plt.tight_layout()
# plt.show()

def metoda_rk4(y0, lam, dt, t_end):
    t_values = np.arange(0, t_end + dt, dt)
    y_values = [y0]
    
    for i in range(1, len(t_values)):
        k1 = lam * y_values[i-1]
        k2 = lam * (y_values[i-1] + (dt/2) * k1)
        k3 = lam * (y_values[i-1] + (dt/2) * k2)
        k4 = lam * (y_values[i-1] + dt*k3)
        y_values.append(y_values[i-1] + dt/6 * (k1 + 2 * k2 + 2 * k3 + k4))

    return t_values, y_values

fig, axs = plt.subplots(1, 2, figsize=(14, 6))
t_analytical = np.linspace(0, t_end, 500)
y_analytical = rozwiazanie_analityczne(t_analytical, lam)
axs[0].plot(t_analytical, y_analytical, 'k--', label='Rozwiazanie analityczne')

for dt in time_steps:
    t_values, y_values = metoda_rk4(y0, lam, dt, t_end)
    axs[0].plot(t_values, y_values, label=f'Rozwiązanie numeryczne (Δt={dt})', linewidth=2)

axs[0].set_xlabel('Czas (t)')
axs[0].set_ylabel('y(t)')
axs[0].set_title('Rozwiązania numeryczne i analityczne')
axs[0].legend()
axs[0].grid(True)

for dt in time_steps:
    t_values, y_values = metoda_rk4(y0, lam, dt, t_end)
    y_exact = rozwiazanie_analityczne(t_values, lam)
    error = y_values - y_exact
    axs[1].plot(t_values, error, label=f'Błąd (Δt={dt})')

axs[1].set_xlabel('Czas (t)')
axs[1].set_ylabel('Błąd δ(t)')
axs[1].set_title('Zmiany błędu globalnego')
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.show()

