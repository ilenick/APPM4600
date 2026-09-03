import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(1.920, 2.080, 160)

def p(x):
    return (x-2)**9
    
def f(x):
    return x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304*x - 512

y = p(x)
g = f(x)

plt.plot(x, y)
plt.plot(x, g)
plt.xlabel('x')
plt.ylabel('p(x)')
plt.show()