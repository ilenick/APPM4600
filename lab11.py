import numpy as np
import matplotlib.pyplot as plt

def lgwt(N,a,b):
    N = N-1
    N1 = N+1
    N2 = N+2
    eps = np.finfo(float).eps
    xu = np.linspace(-1,1,N1)

    # Initial guess
    y = np.cos((2*np.arange(0,N1)+1)*np.pi/(2*N+2))+(0.27/N1)*np.sin(np.pi*xu*N/N2)

    # Legendre-Gauss Vandermonde Matrix
    L = np.zeros((N1,N2))

    # Compute the zeros of the N+1 Legendre Polynomial
    # using the recursion relation and the Newton-Raphson method
    y0 = 2.
    one = np.ones((N1,))

    # Iterate until new points are uniformly within epsilon of old points
    while np.max(np.abs(y-y0)) > eps:

        L[:,0] = one
        L[:,1] = y
        for k in range(2,N1+1):
            L[:,k] = ((2*k-1)*y*L[:,k-1]-(k-1)*L[:,k-2])/k

        lp = N2*(L[:,N1-1]-y*L[:,N2-1])/(1-y**2)

        y0 = y
        y = y0-L[:,N2-1]/lp

    # Linear map from [-1,1] to [a,b]
    x = (a*(1-y)+b*(1+y))/2

    # Compute the weights
    w = (b-a)/((1-y**2)*lp**2)*(N2/N1)**2

    return x,w

def driver():
    
    f = lambda x: np.sin(1/x)
    a = 0.1
    b = 2
    
    # exact integral
    I_ex = 2
    
    N = 5
    # for simpson's n must be even.
    # n+1 = number of pts.
    ntest = np.arange(2,N+1,step=2)

    errorT = np.zeros(len(ntest))
    errorS = np.zeros(len(ntest))

    for j in range(0,len(ntest)):
        n = ntest[j]

        I_trap = CompTrap(a,b,n,f)
        errorT[j] = abs(I_ex-I_trap)

        I_simp = CompSimp(a,b,n,f)
        errorS[j] = abs(I_ex-I_simp)

    print('n = ', ntest[-1])
    print('trapezoidal absolute error = ', errorT[-1])
    print('Simpson absolute error = ', errorS[-1])

    # observed order of convergence from the last two runs
    rateT = np.log(errorT[-1]/errorT[-2])/np.log(ntest[-2]/ntest[-1])
    rateS = np.log(errorS[-1]/errorS[-2])/np.log(ntest[-2]/ntest[-1])
    print('observed order (trap) = ', rateT)
    print('observed order (Simpson) = ', rateS)

    # reference slopes O(h^2) and O(h^4), with h = (b-a)/n
    h = (b-a)/ntest
    plt.loglog(ntest,errorT,'o-',label='Composite trapezoidal')
    plt.loglog(ntest,errorS,'s-',label='Composite Simpson')
    plt.loglog(ntest,errorT[0]*(h/h[0])**2,'k--',label='$O(h^2)$')
    plt.loglog(ntest,errorS[0]*(h/h[0])**4,'k:',label='$O(h^4)$')
    plt.xlabel('n (number of subintervals)')
    plt.ylabel('absolute error')
    plt.legend()
    plt.show()

# get lgwts routine and numpy

def composite_simpson(a, b, f, N):
    if N < 3 or N % 2 == 0:
        raise ValueError("N must be odd and at least 3")


    h = (b - a) / (N - 1)
    x = np.linspace(a, b, N)
    y = f(x)


    I = h/3 * (
        y[0] + y[-1]
        + 4*np.sum(y[1:-1:2])
        + 2*np.sum(y[2:-1:2])
    )


    return I

def eval_gauss_quad(M,a,b,f):
  """
  Non-adaptive numerical integrator for \int_a^b f(x)w(x)dx
  Input:
    M - number of quadrature nodes
    a,b - interval [a,b]
    f - function to integrate
  
  Output:
    I_hat - approx integral
    x - quadrature nodes
    w - quadrature weights

  Currently uses Gauss-Legendre rule
  """
  x,w = lgwt(M,a,b)
  I_hat = np.sum(f(x)*w)
  return I_hat,x,w

def adaptive_quad(a,b,f,tol,M,method):
  """
  Adaptive numerical integrator for \int_a^b f(x)dx
  
  Input:
  a,b - interval [a,b]
  f - function to integrate
  tol - absolute accuracy goal
  M - number of quadrature nodes per bisected interval
  method - function handle for integrating on subinterval
         - eg) eval_gauss_quad, eval_composite_simpsons etc.
  
  Output: I - the approximate integral
          X - final adapted grid nodes
          nsplit - number of interval splits
  """
  # 1/2^50 ~ 1e-15
  maxit = 5
  left_p = np.zeros((maxit,))
  right_p = np.zeros((maxit,))
  s = np.zeros((maxit,1))
  left_p[0] = a; right_p[0] = b;
  # initial approx and grid
  s[0],x,_ = method(M,a,b,f);
  # save grid
  X = []
  X.append(x)
  j = 1;
  I = 0;
  nsplit = 1;
  while j < maxit:
    # get midpoint to split interval into left and right
    c = 0.5*(left_p[j-1]+right_p[j-1]);
    # compute integral on left and right spilt intervals
    s1,x,_ = method(M,left_p[j-1],c,f); X.append(x)
    s2,x,_ = method(M,c,right_p[j-1],f); X.append(x)
    if np.max(np.abs(s1+s2-s[j-1])) > tol:
      left_p[j] = left_p[j-1]
      right_p[j] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j] = s1
      left_p[j-1] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j-1] = s2
      j = j+1
      nsplit = nsplit+1
    else:
      I = I+s1+s2
      j = j-1
      if j == 0:
        j = maxit
  return I,np.unique(X),nsplit



# get adaptive_quad routine and numpy from adaptive_quad.py
# get plot routines
import matplotlib.pyplot as plt

# specify the quadrature method 
# (eval_gauss_quad, eval_composite_trap, eval_composite_simpsons)
method = eval_gauss_quad

# interval of integration [a,b]
a = 0.; b = 1.
# function to integrate and true values
# TRYME: uncomment and comment to try different funcs
#        make sure to adjust I_true values if using different interval!
f = lambda x: np.log(x)**2; I_true = 2; labl = "log^2(x)"
#f = lambda x: 1./(np.power(x,(1./5.))); I_true = 5./4.; labl = '1/(x^{1/5})'
#f = lambda x: np.exp(np.cos(x)); I_true = 2.3415748417130531; labl = 'exp(cos(x))'
#f = lambda x: x**20; I_true = 1./21.; labl = 'x^{20}'
# below is for a=0.1, b = 2
#f = lambda x: np.sin(1./x); I_true = 1.1455808341; labl = 'sin(1/x)'

# absolute tolerance for adaptive quad 
tol = 1e-14
# machine eps in numpy
eps = np.finfo(float).eps

# number of nodes and weights, per subinterval 
Ms = np.arange(2,15); nM = len(Ms)
# storage for error
err_old = np.zeros((nM,))
err_new = np.zeros((nM,))


# loop over quadrature orders
# compute integral with non adaptive and adaptive
# compute errors for both 
for iM in range(nM):
  M = Ms[iM]; 
  # non adaptive routine 
  # Note: the _,_ are dummy vars/Python convention 
  # to store uneeded returns from the routine
  I_old,_,_ = method(M,a,b,f)
  # adaptive routine
  I_new,X,nsplit = adaptive_quad(a,b,f,tol,M,method)
  err_old[iM] = np.abs(I_old-I_true)/I_true
  err_new[iM] = np.abs(I_new-I_true)/I_true 
  # clean the error for nice plots
  if err_old[iM] < eps:
    err_old[iM] = eps 
  if err_new[iM] < eps:
    err_new[iM] = eps
  # save grids for M = 2
  if M == 2:
    mesh = X
  
# plot the old and new error for each f and M
fig,ax = plt.subplots(1,2)
ax[0].semilogy(Ms,err_old,'ro--')
ax[0].set_ylim([1e-16,2]);
ax[0].set_xlabel('$M$')
ax[0].set_title('Non-adaptive')
ax[0].set_ylabel('Relative Error');
ax[1].semilogy(Ms,err_new,'ro--',label=labl)
ax[1].set_ylim([1e-16,2]);
ax[1].set_xlabel('$M$')
ax[1].set_title('Adaptive')
ax[1].legend()
plt.show()

# plot the adaptive mesh for M=2
fig,ax = plt.subplots(1,1)
ax.semilogy(mesh,f(mesh),'ro',label=labl)
ax.legend()
plt.show()


        
def CompTrap(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    
    I_trap = h*f(xnode[0])*1/2
    
    for j in range(1,n):
         I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])
    
    return I_trap     

def CompSimp(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    I_simp = f(xnode[0])

    nhalf = n/2
    for j in range(1,int(nhalf)+1):
         # odd part
         I_simp = I_simp +4*f(xnode[2*j-1])
         # even part (interior only; x_n is added below)
         if j < nhalf:
             I_simp = I_simp+2*f(xnode[2*j])
    I_simp= I_simp + f(xnode[n])
    
    I_simp = h/3*I_simp
    
    return I_simp     

driver()    