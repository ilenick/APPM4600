import numpy as np

tol = 1e-4
Nmax = 100

def F(z):
    x, y = z
    return np.array([
        x**2 + y**2 - 4,
        np.exp(x) + y - 1
    ])

def J(z):
    x, y = z
    return np.array([
        [2*x, 2*y],
        [np.exp(x), 1]
    ])

z0 = np.array([1.0, -1.0])

# Lazy Newton: freeze the Jacobian at the initial guess
def lazy_newton(z0):
    z = z0.copy()
    J0 = J(z0)

    for n in range(Nmax + 1):
        r = np.linalg.norm(F(z))
        if r < tol:
            return z, n, r

        z = z - np.linalg.solve(J0, F(z))

        if not np.all(np.isfinite(z)) or np.linalg.norm(z) > 1e10:
            return z, n + 1, np.inf

    return z, Nmax, np.linalg.norm(F(z))


# Good Broyden: update the Jacobian approximation
def broyden(z0):
    z = z0.copy()
    B = J(z0)

    for n in range(Nmax + 1):
        r = np.linalg.norm(F(z))
        if r < tol:
            return z, n, r

        s = np.linalg.solve(B, -F(z))
        znew = z + s

        if not np.all(np.isfinite(znew)) or np.linalg.norm(znew) > 1e10:
            return znew, n + 1, np.inf

        y = F(znew) - F(z)
        B = B + np.outer(y - B @ s, s) / (s @ s)
        z = znew

    return z, Nmax, np.linalg.norm(F(z))


for name, method in [
    ("Lazy Newton", lazy_newton),
    ("Broyden", broyden)
]:
    z, n, r = method(z0)
    print(f"{name}:")
    print(f"  Iterations: {n}")
    print(f"  Approximation: {z}")
    print(f"  Residual: {r:.3e}")