import math

#second method

def f(x):
    return math.exp(3*x) - 27*x**6 + 27*x**4*math.exp(x) - 9*x**2*math.exp(2*x)

def df(x):
    return (
        3*math.exp(3*x)
        - 162*x**5
        + 108*x**3*math.exp(x)
        + 27*x**4*math.exp(x)
        - 18*x*math.exp(2*x)
        - 18*x**2*math.exp(2*x)
    )


m = 3
x = 4.0

tol = 1e-10
max_iter = 100

for k in range(max_iter):
    x_new = x - m * f(x) / df(x)

    print(f"Iteration {k+1}: x = {x_new:.12f}, |f(x)| = {abs(f(x_new)):.3e}")

    if abs(x_new - x) < tol:
        x = x_new
        iterations = k + 1
        break

    x = x_new

else:
    iterations = max_iter
    print("Modified Newton's method did not converge.")

print("\nRoot:", x)
print("f(root):", f(x))
print("Number of iterations:", iterations)

#third method

def h(x):
    return f(x) / df(x)

# Modified Newton's method
x = 4.0
m = 3
tol = 1e-10
max_iter = 100

for k in range(max_iter):
    x_new = x - m * h(x)

    print(f"Iteration {k+1}: x = {x_new:.12f}, |f(x)| = {abs(f(x_new)):.3e}")

    if abs(x_new - x) < tol:
        x = x_new
        iterations = k + 1
        break

    x = x_new
else:
    iterations = max_iter
    print("Modified Newton's method did not converge.")

print("\nRoot:", x)
print("f(root):", f(x))
print("Number of iterations:", iterations)

