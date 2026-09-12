import numpy as np

def g(x):
    return -np.sin(2*x) + (5/4)*x - (3/4)

tol = 1e-10
max_iter = 1000

guesses = np.linspace(-2, 8, 100)

roots = []

for x in guesses:
    for n in range(max_iter):
        x_new = g(x)

        if abs(x_new - x) < tol:
            break

        x = x_new

    # Only keep converged roots
    if n < max_iter - 1:
        already_found = False
        for r in roots:
            if abs(x_new - r) < 1e-8:
                already_found = True
        if not already_found:
            roots.append(x_new)
            print(f"Root: {x_new:.12f}, Iterations: {n+1}")