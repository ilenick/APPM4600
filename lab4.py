import numpy as np
import matplotlib.pyplot as plt

f = lambda x: np.exp(x**2 + 7*x - 30) - 1
fp = lambda x: (2*x+7)*np.exp(x**2+7*x-30)

def newton(f, fp, bmp, tol, Nmax):
    p = np.zeros(Nmax+1)
    p[0]=bmp
    for it in range(Nmax):
        p1 = bmp-f(bmp)/fp(bmp)
        p[it+1] = p1
        if (abs(bmp)<tol):
            pstar = p1
            info = 0
            return [p, pstar, info, it]
        pstar = p1
        info = 1
        return [p, pstar, info, it]
    return

def bisection(f, a, b, fp):
    fa = f(a)
    fb = f(b)
    if (fa*fb>0):
        ier = 1
        astar = a
        return [astar, ier]

    if (fa == 0):
        astar = a
        ier = 0
        return [astar, ier]

    if (fb == 0):
        astar = b
        ier = 0
        return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a) > abs(fp(d))):
        fd = f(d)
        if (fd == 0):
            astar = d
            ier = 0
            return [astar, ier]
        if (fa*fd < 0):
            b = d
        else:
            a = d
            fa = fd
        d = 0.5*(a+b)
        count = count + 1
        print('abs(d-a) = ', abs(d-a))

    astar = d
    ier = 0
    print('count = ', count)
    return astar

#bi_root = bisection(f, 2, 4.5, 2.89)
#print(bi_root)

#new_root = newton(f, fp, 4.5, 0.0002, 100000)
#print(new_root)

biroot = bisection(f, 2, 4.5, fp)
newroot = newton(f, fp, biroot, 0.001, 10000)
print(newroot)