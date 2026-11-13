
import numpy as np
from scipy.special import gamma


def CompTrap(a, b, n, f):

    h = (b-a)/n
    xnode = a + np.arange(n+1)*h

    I_trap = h*(0.5*f(xnode[0]) +
                np.sum(f(xnode[1:-1])) +
                0.5*f(xnode[n]))

    return I_trap


def driver():

    # Values of x to test
    xvalues = [2, 4, 6, 8, 10]

    # Integration parameters
    a = 0
    b = 40
    h = 0.001
    n = int((b-a)/h)

    print("Step size:", h)
    print("Number of subintervals:", n)
    print("Integration interval:", [a, b])
    print()

    print("{:<8} {:<18} {:<18} {:<15}".format(
        "x", "Trap Approx", "Scipy Gamma", "Relative Error"))

    for x in xvalues:

        # Gamma integrand
        f = lambda t: t**(x-1)*np.exp(-t)

        # Numerical approximation
        I_trap = CompTrap(a, b, n, f)

        # Exact value from scipy
        I_exact = gamma(x)

        # Relative error
        rel_error = abs(I_trap-I_exact)/abs(I_exact)

        print("{:<8} {:<18.10f} {:<18.10f} {:.6e}".format(
            x, I_trap, I_exact, rel_error))


driver()