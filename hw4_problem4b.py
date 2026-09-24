import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**6 - x - 1


def df(x):
    return 6*x**5 - 1

alpha = 1.1347241384015194

# Newton's Method

x = 2
newton_x = []

max_iter = 10

for k in range(max_iter):
    newton_x.append(x)

    x = x - f(x) / df(x)

newton_x = np.array(newton_x)

#  errors
newton_errors = np.abs(newton_x - alpha)

# Secant Method

x0 = 2
x1 = 1

secant_x = [x0, x1]

for k in range(max_iter - 2):

    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))

    secant_x.append(x2)

    x0 = x1
    x1 = x2

secant_x = np.array(secant_x)

# errors
secant_errors = np.abs(secant_x - alpha)

print("NEWTON'S METHOD")
print("-" * 55)
print(" k          x_k                  |x_k - alpha|")
print("-" * 55)

for k in range(len(newton_x)):
    print(f"{k:<5} {newton_x[k]:<22.12f} {newton_errors[k]:.12e}")

print("\nSECANT METHOD")
print("-" * 55)
print(" k          x_k                  |x_k - alpha|")
print("-" * 55)

for k in range(len(secant_x)):
    print(f"{k:<5} {secant_x[k]:<22.12f} {secant_errors[k]:.12e}")

newton_x_error = newton_errors[:-1]
newton_y_error = newton_errors[1:]

secant_x_error = secant_errors[:-1]
secant_y_error = secant_errors[1:]

newton_mask = (
    (newton_x_error > 1e-14) &
    (newton_y_error > 1e-14)
)

secant_mask = (
    (secant_x_error > 1e-14) &
    (secant_y_error > 1e-14)
)


newton_x_error = newton_x_error[newton_mask]
newton_y_error = newton_y_error[newton_mask]

secant_x_error = secant_x_error[secant_mask]
secant_y_error = secant_y_error[secant_mask]

# Calculate slopes

newton_slope = np.polyfit(
    np.log(newton_x_error),
    np.log(newton_y_error),
    1
)[0]

secant_slope = np.polyfit(
    np.log(secant_x_error),
    np.log(secant_y_error),
    1
)[0]


print("\n" + "=" * 55)
print("CONVERGENCE ORDERS")
print("=" * 55)

print(f"Newton slope = {newton_slope:.6f}")
print(f"Expected Newton order = 2")

print()

print(f"Secant slope = {secant_slope:.6f}")
print(f"Expected Secant order = {(1 + np.sqrt(5))/2:.6f}")

plt.figure(figsize=(8, 6))

plt.loglog(
    newton_x_error,
    newton_y_error,
    'o-',
    label=f"Newton (slope = {newton_slope:.3f})"
)

plt.loglog(
    secant_x_error,
    secant_y_error,
    's-',
    label=f"Secant (slope = {secant_slope:.3f})"
)

plt.xlabel(r'$|x_k-\alpha|$')
plt.ylabel(r'$|x_{k+1}-\alpha|$')

plt.title('Convergence Order')

plt.legend()

plt.grid(True, which="both")

plt.show()