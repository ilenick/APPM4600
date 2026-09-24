import math

def f(x):
    return x**6 - x - 1

def df(x):
    return 6*x**5 - 1

# Actual root (to calculate error)
r = 1.1347241384015194

# -------------------------
# Newton's Method
# -------------------------

x = 2
tol = 1e-12
max_iter = 20

print("NEWTON'S METHOD")
print("n        x_n                  error")
print("-" * 50)

for n in range(max_iter):
    error = abs(x - r)
    print(f"{n:<8}{x:<22.12f}{error:.12e}")

    if error < tol:
        break

    x = x - f(x) / df(x)


# -------------------------
# Secant Method
# -------------------------

x0 = 2
x1 = 1

print("\nSECANT METHOD")
print("n        x_n                  error")
print("-" * 50)

for n in range(max_iter):
    error = abs(x1 - r)
    print(f"{n:<8}{x1:<22.12f}{error:.12e}")

    if error < tol:
        break

    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))

    x0 = x1
    x1 = x2