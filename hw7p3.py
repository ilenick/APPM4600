
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 1 / (1 + (10*x)**2)

def barycentric_interpolation(x, y, xf):
    N = len(x)

    
    w = np.ones(N)
    for j in range(N):
        w[j] = 1 / np.prod(x[j] - np.delete(x, j))

    
    w = w / np.max(np.abs(w))

    pf = np.empty_like(xf)

    for k, z in enumerate(xf):
        
        idx = np.where(np.isclose(z, x, atol=1e-14, rtol=0))[0]

        if len(idx) > 0:
            pf[k] = y[idx[0]]
        else:
            terms = w / (z - x)
            pf[k] = np.sum(terms * y) / np.sum(terms)

    return pf

xf = np.linspace(-1, 1, 1001)

for N in [20, 100, 200]:

    j = np.arange(1, N + 1)
    x = np.cos((2*j - 1) * np.pi / (2*N))
    y = f(x)


    pf = barycentric_interpolation(x, y, xf)

    plt.figure(figsize=(9, 5))
    plt.plot(xf, f(xf), label='f(x)')
    plt.plot(xf, pf, label=f'Chebyshev interpolant, N={N}')
    plt.plot(x, y, 'o', markersize=3, label='Chebyshev nodes')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Interpolation at Chebyshev Nodes (N={N})')
    plt.grid(True)
    plt.legend()
    plt.show()

    print(f"N={N}, max interpolation error = "
          f"{np.max(np.abs(pf - f(xf))):.6e}")