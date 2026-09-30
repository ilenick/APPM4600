import numpy as np

def f(x, y):
    return 3*x**2 - y**2

def g(x, y):
    return 3*x*y**2 - x**3 - 1

x = 1.0
y = 1.0

for n in range(21):

    print(n, x, y)

    fx = f(x, y)
    gx = g(x, y)

    x_new = x - (1/6)*fx - (1/18)*gx
    y_new = y - (1/6)*gx

    x = x_new
    y = y_new