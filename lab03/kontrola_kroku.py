
import numpy as np
import matplotlib.pyplot as plt

alpha = 5.0
x0 = 0.01
v0 = 0.0
dt0 = 1.0
S = 0.75
p = 2  
t_max = 40
TOL_1 = 1e-2
TOL_2 = 1e-5


def f(x, v):
    return v


def g(x, v, alpha):
    return alpha * (1 - x**2) * v - x


def metoda_trapezow(xn, vn, dt, alpha):
    x_next, v_next = xn, vn  
    tol = 1e-10

    while 1:
        F = x_next - xn - dt / 2 * (f(xn, vn) + f(x_next, v_next))
        G = v_next - vn - dt / 2 * (g(xn, vn, alpha) + g(x_next, v_next, alpha))

        a11 = 1
        a12 = -dt / 2
        a21 = -dt / 2 * (-2 * alpha * x_next * v_next - 1)
        a22 = 1 - dt / 2 * alpha * (1 - x_next**2)

        det = a11 * a22 - a12 * a21

        dx = (-F * a22 + G * a12) / det
        dv = (-a11 * G + a21 * F) / det

        x_next += dx
        v_next += dv

        if abs(dx) < tol and abs(dv) < tol:
            break

    return x_next, v_next


def metoda_rk2(xn, vn, dt, alpha):
    k1x = f(xn, vn)
    k1v = g(xn, vn, alpha)

    k2x = f(xn + dt * k1x, vn + dt * k1v)
    k2v = g(xn + dt * k1x, vn + dt * k1v, alpha)

    x_next = xn + dt / 2 * (k1x + k2x)
    v_next = vn + dt / 2 * (k1v + k2v)

    return x_next, v_next


def kontrola_kroku(method, x0, v0, dt0, t_max, alpha, tol):
    t_array = []
    x_array = []
    v_array = []
    dt_array = []

    t = 0.0
    dt = dt0
    xn = x0
    vn = v0
    
    while t < t_max:
        # dwa kroki dt
        twosteps_x_n1, twosteps_v_n1 = method(xn, vn, dt, alpha)
        twosteps_x_n2, twosteps_v_n2 = method(twosteps_x_n1, twosteps_v_n1, dt, alpha)

        #jeden krok 2 * dt
        onestep_x_n2, onestep_v_n2 = method(xn, vn, 2 * dt, alpha)

        E_x = (twosteps_x_n2 - onestep_x_n2) / (2**p - 1)
        E_v = (twosteps_v_n2 - onestep_v_n2) / (2**p - 1)

        if max(abs(E_x), abs(E_v)) < tol:
            t = t + 2 * dt
            xn = twosteps_x_n2
            vn = twosteps_v_n2
            t_array.append(t)
            dt_array.append(dt)
            x_array.append(xn)
            v_array.append(vn)

        dt = (S * tol / (max(abs(E_x), abs(E_v))))**(1 / (p + 1)) * dt

    return (t_array), (x_array), (v_array), (dt_array) 




t_rk2_tol1, x_rk2_tol1, v_rk2_tol1, dt_rk2_tol1 = kontrola_kroku(
    metoda_rk2, x0, v0, dt0, t_max, alpha, TOL_1
)
t_rk2_tol2, x_rk2_tol2, v_rk2_tol2, dt_rk2_tol2 = kontrola_kroku(
    metoda_rk2, x0, v0, dt0, t_max, alpha, TOL_2
)
t_trap_tol1, x_trap_tol1, v_trap_tol1, dt_trap_tol1 = kontrola_kroku(
    metoda_trapezow, x0, v0, dt0, t_max, alpha, TOL_1
)
t_trap_tol2, x_trap_tol2, v_trap_tol2, dt_trap_tol2 = kontrola_kroku(
    metoda_trapezow, x0, v0, dt0, t_max, alpha, TOL_2
)

plt.figure(figsize=(12, 8))

plt.subplot(4, 2, 1)
plt.plot(t_rk2_tol1, x_rk2_tol1, label="TOL=10e-2")
plt.plot(t_rk2_tol2, x_rk2_tol2, label="TOL=10-5")
plt.title("metoda RK2")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.legend()
plt.grid()

plt.subplot(4, 2, 2)
plt.plot(t_rk2_tol1, v_rk2_tol1, label="TOL=1e-2")
plt.plot(t_rk2_tol2, v_rk2_tol2, label="TOL=1e-5")
plt.title("metoda RK2")
plt.xlabel("t")
plt.ylabel("v(t)")
plt.legend()
plt.grid()

plt.subplot(4, 2, 3)
plt.plot(t_rk2_tol1, dt_rk2_tol1, label="TOL=1e-2")
plt.plot(t_rk2_tol2, dt_rk2_tol2, label="TOL=1e-5")
plt.title("metoda RK2")
plt.xlabel("t")
plt.ylabel("dt(t)")
plt.legend()
plt.grid()

plt.subplot(4, 2, 4)
plt.plot(x_rk2_tol1, v_rk2_tol1, label="TOL=1e-2")
plt.plot(x_rk2_tol2, v_rk2_tol2, label="TOL=1e-5")
plt.title("metoda RK2")
plt.xlabel("x(t)")
plt.ylabel("v(t)")
plt.legend()
plt.grid()

plt.subplot(4, 2, 5)
plt.plot(t_trap_tol1, x_trap_tol1, label="TOL=1e-2")
plt.plot(t_trap_tol2, x_trap_tol2, label="TOL=1e-5")
plt.title("metoda trapezow")
plt.xlabel("x(t)")
plt.ylabel("v(t)")
plt.legend()
plt.grid()

plt.subplot(4, 2, 6)
plt.plot(t_trap_tol1, v_trap_tol1, label="TOL=1e-2")
plt.plot(t_trap_tol2, v_trap_tol2, label="TOL=1e-5")
plt.title("metoda trapezow")
plt.xlabel("t")
plt.ylabel("v(t)")
plt.legend()
plt.grid()

plt.subplot(4, 2, 7)
plt.plot(t_trap_tol1, dt_trap_tol1, label="TOL=1e-2")
plt.plot(t_trap_tol2, dt_trap_tol2, label="TOL=1e-5")
plt.title("metoda trapzow")
plt.xlabel("t")
plt.ylabel("dt(t)")
plt.legend()
plt.grid()

plt.subplot(4, 2, 8)
plt.plot(x_trap_tol1, v_trap_tol1, label="TOL=1e-2")
plt.plot(x_trap_tol2, v_trap_tol2, label="TOL=1e-5")
plt.title("metoda trapezow")
plt.xlabel("x(t)")
plt.ylabel("v(t)")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()


