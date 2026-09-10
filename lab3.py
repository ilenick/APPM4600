import numpy as np
import matplotlib.pyplot as plt

Nmax = 100
tol = 10e-6
f1 = lambda x: 1+0.5*np.sin(x)

def fixedpt(f, x0, tol, Nmax):
    x = np.zeros(Nmax)
    count = 0
    while (count < Nmax):
        count = count +1
        x1 = f(x0)
        x[count] = x1
        if (abs(x1-x0) < tol):
            xstar = x1
            ier = 0
            return [xstar,ier]
        x0 = x1

    xstar = x1
    ier = 1
    return [x]

x = fixedpt(f1, 1, tol, Nmax)

print(x)

################################################ 
# PSUEDOCODE FOR STEPHENSON METHOD

#def p(n)
#def b = g(p(n))
#def c = g(b)
# for (n in Nmax) {
#    p(n+1) = a - (b-a)**2 / (c - 2*b + a)
# }
# return p(n+1)

## We need a sequence p(n), and a function g(x) for us to plug values into. We should output the value that p(n) has converged to.
################################################

n = 0
def g(x):
    g = (10 / (x+4))**(1/2)
    return

def p(n):
    1 + n
    return

pn = p(n) # just for testing 

def b(pn):
    g(pn)
    return

def c(g(x))
    return g(b)

def stephenson(g(x), p(n), a, b, c):
    p(n+1) = a - (b-a)^2 / (c - 2*b + a)
    return