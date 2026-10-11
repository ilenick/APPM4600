
import numpy as np

def F(v):
    x, y, z = v
    return np.array([
        x + np.cos(x*y*z) - 1,
        (1-x)**0.25 + y + 0.05*z**2 - 0.15*z - 1,
        -x**2 - 0.1*y**2 + 0.01*y + z - 1
    ])

def J(v):
    x, y, z = v
    return np.array([
        [1-y*z*np.sin(x*y*z),
         -x*z*np.sin(x*y*z),
         -x*y*np.sin(x*y*z)],
        [-0.25*(1-x)**(-0.75),
         1, 0.1*z-0.15],
        [-2*x, -0.2*y+0.01, 1]
    ])

def newton(v, tol=1e-6, Nmax=100):
    v = np.array(v, dtype=float)

    for k in range(Nmax + 1):
        r = np.linalg.norm(F(v))
        if r < tol:
            return v, k, r

        step = np.linalg.solve(J(v), -F(v))
        alpha = 1.0

        # Backtrack if necessary to stay in the domain
        # and reduce the residual.
        while alpha > 1e-12:
            vn = v + alpha*step
            if (vn[0] < 1 and np.all(np.isfinite(vn))
                    and np.linalg.norm(F(vn)) < r):
                break
            alpha *= 0.5
        else:
            return v, k, r

        v = vn

    return v, Nmax, np.linalg.norm(F(v))

def steepest_descent(v, tol=1e-6, Nmax=100000):
    v = np.array(v, dtype=float)

    for k in range(Nmax + 1):
        fv = F(v)
        r = np.linalg.norm(fv)

        if r < tol:
            return v, k, r

        # Gradient of phi = 1/2 ||F||^2
        grad = J(v).T @ fv
        phi = 0.5 * r**2
        alpha = 1.0

        # Armijo backtracking line search
        while alpha > 1e-14:
            vn = v - alpha*grad

            if vn[0] < 1 and np.all(np.isfinite(vn)):
                fn = F(vn)
                if (np.all(np.isfinite(fn)) and
                    0.5*np.dot(fn, fn)
                    <= phi - 1e-4*alpha*np.dot(grad, grad)):
                    break

            alpha *= 0.5
        else:
            return v, k, r

        v = vn

    return v, Nmax, np.linalg.norm(F(v))


# Same initial guess for all methods
v0 = np.array([0.5, 0.5, 1.5])

# Method 1: Newton, tolerance 1e-6
v1, n1, r1 = newton(v0)

# Method 2: Steepest descent, tolerance 1e-6
v2, n2, r2 = steepest_descent(v0, tol=1e-6)

# Method 3: Steepest descent to 5e-2,
# then Newton to 1e-6
v_coarse, n_sd, r_sd = steepest_descent(
    v0, tol=5e-2
)
v3, n_newton, r3 = newton(v_coarse, tol=1e-6)

print("Method 1: Newton")
print("Solution:", v1)
print("Iterations:", n1, "Residual:", r1)

print("\nMethod 2: Steepest descent")
print("Solution:", v2)
print("Iterations:", n2, "Residual:", r2)

print("\nMethod 3: Steepest descent + Newton")
print("Steepest descent result:", v_coarse)
print("Steepest descent iterations:", n_sd)
print("Newton iterations:", n_newton)
print("Total iterations:", n_sd + n_newton)
print("Final solution:", v3)
print("Final residual:", r3)
