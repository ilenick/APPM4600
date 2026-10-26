import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import (
    BarycentricInterpolator,
    KroghInterpolator,
    CubicSpline
)

def f(x):
    return 1 / (1 + x**2)

def df(x):
    return -2*x / (1 + x**2)**2

ns = [5, 10, 15, 20]
xx = np.linspace(-5, 5, 4001)
yy = f(xx)

for n in ns:
    x = np.linspace(-5, 5, n)
    y = f(x)

    # Lagrange polynomial
    lagrange = BarycentricInterpolator(x, y)

    # Hermite polynomial: repeat each node and provide
    # the function value followed by its first derivative
    hermite_nodes = np.repeat(x, 2)
    hermite_values = np.column_stack(
        (y, df(x))
    ).ravel()
    hermite = KroghInterpolator(
        hermite_nodes, hermite_values
    )

    # Cubic splines
    natural = CubicSpline(x, y, bc_type="natural")
    clamped = CubicSpline(
        x, y,
        bc_type=((1, df(x[0])), (1, df(x[-1])))
    )

    methods = [
        ("Lagrange", lagrange),
        ("Hermite", hermite),
        ("Natural spline", natural),
        ("Clamped spline", clamped)
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"Interpolation with n = {n} nodes")

    for ax, (name, method) in zip(axes.flat, methods):
        ax.plot(xx, yy, label="f(x)")
        ax.plot(xx, method(xx), label=name)
        ax.plot(x, y, "o", label="Nodes")
        ax.set_xlim(-5, 5)
        ax.set_ylim(-1, 2)
        ax.set_title(name)
        ax.grid(True)
        ax.legend()

    plt.tight_layout()
    plt.show()

    print(f"\nn = {n}")
    for name, method in methods:
        error = np.max(np.abs(method(xx) - yy))
        print(f"{name:16s}: max error = {error:.6g}")