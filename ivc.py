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

def AdamsBash2():
  w2 = np.zeros(n+1)
  w2[0] = w[0]
  w2[1] = w[1]
  for i in range(2, n+1):
    w2[i] = w2[i-1]+h/2*(3*f(t[i-1], w2[i-1])-f(t[i-2], w2[i-2]))
  return w2

def AdamsBash3():
  w3 = np.zeros(n+1)
  w3[0] = w[0]
  w3[1] = w[1]
  w3[2] = w[2]
  for i in range(3, n+1):
    w3[i] = w3[i-1] + h/12*(23*f(t[i-1], w3[i-1])- 16*f(t[i-2], w3[i-2])+5*f(t[i-3], w3[i-3]))
  return w3

def AdamsBash4():
  w4 = np.zeros(n+1)
  w4[0] = w[0]
  w4[1] = w[1]
  w4[2] = w[2]
  w4[3] = w[3]
  for i in range(4, n+1):
    w4[i] = w4[i-1] + h/24 *(55*f(t[i-1], w4[i-1])- 59*f(t[i-2], w4[i-2])+37*f(t[i-3], w[i-3])- 9*f(t[i-4], w[i-4]))
  return w4

def AdamsBash5():
  w5 = np.zeros(n+1)
  w5[0] = w[0]
  w5[1] = w[1]
  w5[2] = w[2]
  w5[3] = w[3]
  w5[4] = w[4]
  for i in range(5, n+1):
    w5[i] =  w5[i-1] + h/720*(1901*f(t[i-1], w5[i-1])-2774*f(t[i-2], w5[i-2])+2616*f(t[i-3], w5[i-3])- 1274*f(t[i-4], w5[i-4])+
                              251*f(t[i-5], w5[i-5]))
  return w5
