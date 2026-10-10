import numpy as np
import matplotlib.pyplot as plt

A = np.array([[0.16, -0.17],
              [0.52, -0.26]])

def G(z):
    x, y = z
    F = np.array([
        3*x**2 + 4*y**2 - 1,
        y**3 - 8*x**3 - 1
    ])
    return z - A @ F

initial_guesses = [
    (-0.5, 0.25),
    (-1.5, 1.25),
    (-3, 3)
]

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for ax, z0 in zip(axes, initial_guesses):
    z = np.array(z0, dtype=float)
    xs, ys = [z[0]], [z[1]]

    for n in range(30):
        z = G(z)
        xs.append(z[0])
        ys.append(z[1])

        if not np.all(np.isfinite(z)) or np.linalg.norm(z) > 1e6:
            break

    ax.plot(xs, ys, 'o-', markersize=3)
    ax.set_title(f"Initial guess: {z0}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

plt.tight_layout()
plt.show()