import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


# Inputs
a = 1  # Start
b = 2  # Stop
h = 0.05  # Step size
w0 = 0.5  # Initial value


def f(t, y):  # IVP function
    return -3 * t * y - y**2


def t():  # Calculte list of time stamps
    n = int((b - a) / h)
    t_list = np.zeros(n + 1)
    for i in range(n + 1):
        t_list[i] = a + h * i
    return t_list


def forward_euler():
    n = int((b - a) / h)
    t_list = t()
    w = np.zeros(n + 1)  # Initialize empty solution list
    w[0] = w0  # Initial value
    for j in range(1, n + 1):
        w[j] = w[j - 1] + h * f(t_list[j - 1], w[j - 1])
    return w


def backward_euler():
    n = int((b - a) / h)
    t_list = t()
    w = np.zeros(n + 1)  # Initialize empty solution list
    w[0] = w0  # Initial value
    for i in range(1, n + 1):
        def func(y): return y - w[i - 1] - h * f(t_list[i], y)
    guess = forward_euler()[i]
    y = fsolve(func, guess)[0]
    w[i] = y
    return w


def trapzoidal():
    n = int((b - a) / h)
    t_list = t()
    w = np.zeros(n + 1)  # Initialize empty solution list
    w[0] = w0  # Initial value
    for i in range(1, n + 1):
        def func(y): return y - w[i - 1] - (h / 2) * \
            (f(t_list[i - 1], w[i - 1]) + f(t_list[i], y))
    guess = forward_euler()[i]
    y = fsolve(func, guess)[0]
    w[i] = y
    return w


def leap_frog():
    n = int((b - a) / h)
    t_list = t()
    w = np.zeros(n + 1)
    w[0] = w0
    w[1] = forward_euler()[1]
    for i in range(2, n + 1):
        w[i] = w[i - 2] + 2 * h * f(t_list[i - 1], w[i - 1])
    return w


def print_table():
    return pd.DataFrame({'Time Step': t(), 'Forward Euler': forward_euler(),
                         'Backward Euler': backward_euler(), 'Trapzoidal': trapzoidal(), 'Leap Frog': leap_frog()})


def rk2():
    w = np.zeros(n + 1)
    t_list = t()
    w[0] = w0
    for i in range(1, n + 1):
        k1 = f(t_list[i - 1], w[i - 1])
        k2 = f(t_list[i - 1] + h, w[i - 1] + h * k1)
        w[i] = w[i - 1] + k / 2 * (k1 + k2)
    return w


def rk4():
    w = np.zeros(n + 1)
    t_list = t()
    w[0] = w0
    for i in range(1, n + 1):
        k1 = h * f(t_list[i - 1], w[i - 1])
        k2 = h * f(t_list[i - 1] + h / 2, w[i - 1] + k1 / 2)
        k3 = h * f(t_list[i - 1] + h / 2, w[i - 1] + k2 / 2)
        k4 = h * f(t_list[i - 1] + h, w[i - 1] + k3)
        w[i] = w[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return w
