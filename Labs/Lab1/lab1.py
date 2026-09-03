import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
import math


#######################
#3.2 excersises
#######################



#3.2.1
#x = np.linspace(0, 2*np.pi, 10)
#y = np.arange(0, 2*np.pi, 10)





# 3.2.2
#x[0:3] #each entry for linspace is x[a] where a is the index
#y[0:3] #same method for accessing each index of y





#3.2.3
#print('The first 3 entries of x are: ', x[0:3])





#3.2.4/3.2.5
#w = 10**(-np.linspace(1, 10, 10))
        #entries are w[0], w[1], ..., w[9]
#x = np.linspace(1, 10, 10)
#s = 3*(10**(-np.linspace(1, 10, 10)))


#plt.semilogy(x, w)
#plt.semilogy(x, s)
#plt.xlabel('x')
#plt.ylabel('w')
#plt.title('Plot of two vectors in semilogy')
#plt.show()





####################
#4.2 excersises
####################

#4.2.1

def driver():
    n = 100
    x = np.linspace(0,np.pi,n)
# this is a function handle. You can use it to define
# functions instead of using a subroutine like you
# have to in a true low level language.
    f = lambda x: x**2 + 4*x + 2*np.exp(x)
    g = lambda x: 6*x**3 + 2*np.sin(x)
    y = np.sin(x)
    w = np.cos(x)
# evaluate the dot product of y and w
    dp = dotProduct(y,w,n)
# print the output
    print('the dot product is : ', dp)
    return


def dotProduct(x,y,n):
# Computes the dot product of the n x 1 vectors x and y
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    return dp
driver()

#4.2.1
# I changed the vectors to y = sin(x), w=cos(X)

A = np.array([
    [1, 2],
    [3, 4]
])

def matrixVectorProduct(X, Y):
    mvp = []

    return mvp

#This is all I finished in lab, I will complete it later though!