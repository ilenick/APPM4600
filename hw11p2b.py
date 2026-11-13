
import matplotlib.pyplot as plt
import numpy as np


def driver():

    f = lambda x: np.sin(x)
    a = 0
    b = np.pi

    # Exact integral
    I_ex = 2

    nrefine = 5

    # Romberg trapezoidal rule
    RomT = rombergTrap(a, b, nrefine, f)

    # Romberg Simpson's rule
    RomS = rombergSimp(a, b, nrefine, f)

    # Number of panels
    ntestT = 2**np.arange(0, nrefine+1)
    ntestS = 2*2**np.arange(0, nrefine+1)

    # Absolute errors along the diagonal
    errT = abs(I_ex - np.diag(RomT))
    errS = abs(I_ex - np.diag(RomS))

    print("Romberg Trapezoidal Errors:")
    print(errT)

    print("\nRomberg Simpson Errors:")
    print(errS)

    # Plot comparison
    plt.figure()
    plt.semilogy(ntestT, errT, 'o-', label='Romberg Trapezoidal')
    plt.semilogy(ntestS, errS, 's-', label='Romberg Simpson')

    plt.xlabel('Number of subintervals')
    plt.ylabel('Absolute error')
    plt.title('Romberg Integration Comparison')
    plt.legend()
    plt.grid(True)
    plt.show()


def rombergTrap(a, b, nrefine, f):

    ntest = 2**np.arange(0, nrefine+1)
    nn = len(ntest)

    Rom = np.zeros((nn, nn))

    # First column: composite trapezoidal rule
    for j in range(0, nn):
        Rom[j, 0] = CompTrap(a, b, ntest[j], f)

    # Richardson extrapolation
    for j in range(1, nn):
        for i in range(j, nn):
            Rom[i, j] = Rom[i, j-1] + (
                Rom[i, j-1] - Rom[i-1, j-1]
            )/(4**j - 1)

    return Rom


def rombergSimp(a, b, nrefine, f):

    # Simpson requires an even number of panels
    ntest = 2*2**np.arange(0, nrefine+1)

    nn = len(ntest)

    Rom = np.zeros((nn, nn))

    # First column: composite Simpson's rule
    for j in range(0, nn):
        Rom[j, 0] = CompSimp(a, b, ntest[j], f)

    # Richardson extrapolation
    for j in range(1, nn):
        denom = 2**(2*j+2) - 1

        for i in range(j, nn):
            Rom[i, j] = Rom[i, j-1] + (
                Rom[i, j-1] - Rom[i-1, j-1]
            )/denom

    return Rom


def CompTrap(a, b, n, f):

    h = (b-a)/n
    xnode = a + np.arange(0, n+1)*h

    I_trap = h*(0.5*f(xnode[0]) +
                np.sum(f(xnode[1:-1])) +
                0.5*f(xnode[n]))

    return I_trap


def CompSimp(a, b, n, f):

    h = (b-a)/n
    xnode = a + np.arange(0, n+1)*h

    I_simp = (h/3)*(f(xnode[0]) +
                    f(xnode[n]) +
                    4*np.sum(f(xnode[1:-1:2])) +
                    2*np.sum(f(xnode[2:-1:2])))

    return I_simp


driver()
