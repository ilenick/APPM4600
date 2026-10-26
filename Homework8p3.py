
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
    # Chebyshev-Lobatto nodes include both endpoints.
    # This lets both cubic splines cover [-5, 5].
    x = np.sort(
        5 * np.cos(np.arange(n) * np.pi / (n - 1))
    )
    y = f(x)

    # Lagrange polynomial
    lagrange = BarycentricInterpolator(x, y)

    # Hermite polynomial: function values and derivatives
    hermite_nodes = np.repeat(x, 2)
    hermite_values = np.column_stack(
        (y, df(x))
    ).ravel()
    hermite = KroghInterpolator(
        hermite_nodes, hermite_values
    )

    # Natural and clamped cubic splines
    natural = CubicSpline(x, y, bc_type="natural")
    clamped = CubicSpline(
        x, y,
        bc_type=((1, df(-5)), (1, df(5)))
    )

    methods = [
        ("Lagrange", lagrange),
        ("Hermite", hermite),
        ("Natural spline", natural),
        ("Clamped spline", clamped)
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"Chebyshev interpolation, n = {n}")

    print(f"\nn = {n}")
    for ax, (name, method) in zip(axes.flat, methods):
        approx = method(xx)
        error = np.max(np.abs(approx - yy))

        ax.plot(xx, yy, label="f(x)")
        ax.plot(xx, approx, label=name)
        ax.plot(x, y, "o", label="Nodes")
        ax.set_title(name)
        ax.grid(True)
        ax.legend()

        print(f"{name:16s}: max error = {error:.6g}")

    plt.tight_layout()
    plt.show()
