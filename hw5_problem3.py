import numpy as np

def f(x):
    return x[0]**2 + 4*x[1]**2 + 4*x[2]**2 - 16


def G(x):
    grad = np.array([
        2*x[0],
        8*x[1],
        8*x[2]
    ])

    d = f(x) / np.dot(grad, grad)

    return x - d*grad

x0 = np.array([1.0, 1.0, 1.0])

xstar = x0.copy()

for k in range(20):
    xstar = G(xstar)
print("x* =", xstar)
print()
x = x0.copy()

print(" n       ||e_n||              ||e_(n+1)||/||e_n||^2")
print("-"*65)

for n in range(6):
    e_n = np.linalg.norm(x - xstar)
    x_next = G(x)
    e_next = np.linalg.norm(x_next - xstar)
    ratio = e_next / e_n**2
    print(f"{n:2d}   {e_n:.12e}        {ratio:.12e}")
    x = x_next