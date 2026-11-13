
import numpy as np
from scipy.special import gamma
from numpy.polynomial.laguerre import laggauss

# Number of Gauss-Laguerre nodes
n = 5

# Obtain nodes and weights
nodes, weights = laggauss(n)

# Values of x to test
xvalues = [2, 4, 6, 8, 10]

print("x       Gauss-Laguerre    SciPy Gamma       Relative Error")

for x in xvalues:

    # Gauss-Laguerre approximation
    I_GL = np.sum(weights * nodes**(x-1))

    # Exact Gamma function
    I_exact = gamma(x)

    # Relative error
    rel_error = abs(I_GL-I_exact)/abs(I_exact)

    print(f"{x:<7} {I_GL:<17.10f} {I_exact:<17.10f} {rel_error:.6e}")
