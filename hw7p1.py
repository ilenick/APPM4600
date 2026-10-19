
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BarycentricInterpolator

def f(x):
    return 1 / (1 + (10*x)**2)

N = 20
h = 2 / (N - 1)

x = -1 + np.arange(N) * h
y = f(x)

p = BarycentricInterpolator(x, y)

xf = np.linspace(-1, 1, 1001)

plt.plot(x, y, 'o', label='Data points')
plt.plot(xf, p(xf), label='Interpolating polynomial')
plt.plot(xf, f(xf), label='f(x)')

plt.xlabel('x')
plt.ylabel('y')
plt.title("Interpolation Plot")
plt.grid(True)
plt.legend()
plt.show()

print("Maximum value of p(x):", np.max(p(xf)))