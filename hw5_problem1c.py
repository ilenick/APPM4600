import numpy as np

def f(x, y):
    return 3*x**2 - y**2

def g(x, y):
    return 3*x*y**2 - x**3 - 1

def jacobian(x, y):
    return np.array([
        [6*x, -2*y],
        [3*y**2 - 3*x**2, 6*x*y]
    ])

x = 1.0
y = 1.0

tol = 1e-10
max_iter = 100

print(f"{'n':>3} {'x_n':>15} {'y_n':>15} {'||F||':>15}")

for n in range(max_iter + 1):

    F = np.array([f(x, y), g(x, y)])

    print(f"{n:3d} {x:15.10f} {y:15.10f} {np.linalg.norm(F):15.6e}")

    # Check convergence
    if np.linalg.norm(F) < tol:
        break

    J = jacobian(x, y)

    # Solve J * delta = F
    delta = np.linalg.solve(J, F)

    # Newton step
    x = x - delta[0]
    y = y - delta[1]