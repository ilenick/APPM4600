
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb

def f(x):
    return 1 / (1 + (10*x)**2)

def barycentric_interpolation(x, y, xf):
    N = len(x)

    
    w = (-1.0)**np.arange(N) * comb(N - 1, np.arange(N))

    pf = np.empty_like(xf, dtype=float)

    for k, z in enumerate(xf):
        
        idx = np.where(np.isclose(z, x, rtol=0, atol=1e-14))[0]

        if len(idx) > 0:
            pf[k] = y[idx[0]]
        else:
            terms = w / (z - x)
            pf[k] = np.sum(terms * y) / np.sum(terms)

    return pf


xf = np.linspace(-1, 1, 1001)


for N in [20, 100]:
    x = np.linspace(-1, 1, N)
    y = f(x)

    pf = barycentric_interpolation(x, y, xf)

    plt.figure(figsize=(9, 5))
    plt.plot(xf, f(xf), label='f(x)')
    plt.plot(xf, pf, label=f'Barycentric interpolant, N={N}')
    plt.plot(x, y, 'o', markersize=3, label='Data points')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Barycentric Lagrange Interpolation (N={N})')
    plt.grid(True)
    plt.legend()
    plt.ylim(-2, 2)
    plt.show()

    print(f"N={N}: max |p(x)| on fine grid = {np.max(np.abs(pf)):.6g}")
