import numpy as np

A = np.array([[0.16, -0.17],
              [0.52, -0.26]])

def F(z):
    x, y = z
    return np.array([
        3*x**2 + 4*y**2 - 1,
        y**3 - 8*x**3 - 1
    ])

def JF(z):
    x, y = z
    return np.array([
        [6*x, 8*y],
        [-24*x**2, 3*y**2]
    ])

def fixed_point(z):
    return z - A @ F(z)

def newton(z):
    return z - np.linalg.solve(JF(z), F(z))

initial_guesses = [
    (-0.5, 0.25),
    (-1.5, 1.25),
    (-3.0, 3.0)
]

for z0 in initial_guesses:
    print(f"\nInitial guess: {z0}")

    for name, method in [
        ("Fixed point", fixed_point),
        ("Newton", newton)
    ]:
        z = np.array(z0, dtype=float)

        for n in range(100):
            if np.linalg.norm(F(z)) < 1e-4:
                break

            try:
                z_new = method(z)
            except np.linalg.LinAlgError:
                break

            if not np.all(np.isfinite(z_new)) or np.linalg.norm(z_new) > 1e10:
                break

            z = z_new

        print(f"{name}: after {n} iterations, z = {z}")
        print(f"Residual norm: {np.linalg.norm(F(z)):.3e}")