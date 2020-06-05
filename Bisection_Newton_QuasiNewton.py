#
# name:
# email: 
#

import numpy as np
import math  

#----------------------------------------------------
def bisection(f, x1, x2, tol=1e-14, itmax=200, SHOW=False):
    """ find a solution of f(x)=0 on the interval [x1,x2] by the bisection method. 
        return the results in a tuple (x, k, f(x)), 
                where x is the root, k is the number of iterations.
        (your code should not exceed itmax iterations. if SHOW is passed as True, then
         your code should print out intermediate results)
    """  
    if SHOW:  print("\n-------- bisection ------------")

    if f(x1) * f(x2) > 0:
      print("Out of bound values")
      return x1, 0, f(x1) 
    xm = x1 
    for i in range(itmax):
      xm = (x1+x2)/2
      if(x2-x1 <= tol): 
        break
      if (f(xm) == 0.0): 
        break
      if (f(xm) * f(x1) < 0): 
          x2 = xm 
      else: 
          x1 = xm  
      print("k = " + str(i) + " x = " + str(xm) + " f(x) = " + str(f(xm)))
    return xm, i, f(xm)
 
def newton(f, fderivative, x, tol=1e-14, itmax=100, SHOW=False):   
    """ Newton's method, starting from an initial estimate x of the root,  
        (your code should not exceed itmax iterations.if SHOW is passed as True, then
         your code should print out intermediate results)
    """
    if SHOW:  print("\n-------- Newton's method ------------")   
    for i in range(itmax): 
        x -= f(x) / fderivative(x)
        if((f(x) / fderivative(x)) <= tol): 
            break 
        print("k = " + str(i) + " x = " + str(x) + " f(x) = " + str(f(x)))
    return x, i, f(x)

def quasi_newton(f, x, h=1e-4, tol=1.0e-14, FD='CFD', itmax=120, SHOW=False):   
    """ Quasi Newton's method, starting from an initial value x.
        The derivative is automatically computed via numerical differentiation,
        it defaults to using central finite difference scheme.
        The default for the grid length used for FD is h=1e-4.
        (your code should not exceed itmax iterations. if SHOW is passed as True, then
         your code should print out intermediate results)
    """
    if SHOW:  print("\n-------- Quasi Newton's method,  FD scheme ={} ------------".format(FD))        
    fder = 1
    if(FD == 'FFD'):
        fder = f(x+h) - f(x)
    elif(FD == 'BFD'):
        fder = f(x) - f(x-h)
    else:
        fder = f(x+(h/2)) - f(x-(h/2))
    h1 = f(x) / fder 
    for i in range(itmax):
        if(abs(h1) <= tol): 
            break
        if(FD == 'FFD'):
            fder = f(x+h) - f(x)
        elif(FD == 'BFD'):
            fder = f(x) - f(x-h)
        else:
            fder = f(x+(h/2)) - f(x-(h/2))
        h1 = f(x) / fder
        x -= h1 
        print("k = " + str(i) + " x = " + str(x) + " f(x) = " + str(f(x)))
    return x, i, f(x) 

def solver_tests():
    '''test your code on the equations listed in project pdf, remember to pass SHOW=True to your code
       to print out results at each iteration.
       for this problem, do not introduce any new variables below.
       you only need to fill the lines that have '=' but with empty right-hand-side (using proper
       function names above and the variables that have been assigned)
    '''

    #f1(x) =  x^9 - 8x^6 + 5x + 10 
    #define this function using variable name f1 (try to use lambda function for all problems here)
    f1 = lambda x: math.pow(x,9) - 8*math.pow(x,6) + 5*x + 10
    #add code to call bisection to solve f1(x)=0 on [-1, 1]
    (x, k, fx) = bisection(f1, -1, 1)
    print('\n By bisection: x={},  k={},  fx={}\n'.format(x, k, fx))
    
    #add code to call newton to solve f1(x)=0, set initial root as -1
    #you need to manually get the derivative of f1, pass it as a function named f1der 
    f1der = lambda x: 9*math.pow(x,8)- 48*math.pow(x,5) + 5 
    (x, k, fx) = newton(f1, f1der, -1)  
    print('\n By Newton: x={},  k={},  fx={}\n'.format(x, k, fx))    

    #add code to call quasi_newton to solve f1(x)=0, set initial root as -1
    (x, k, fx) = quasi_newton(f1, -1, .5)
    print('\n By quasi_Newton: x={},  k={},  fx={}\n'.format(x, k, fx))    

    #f2(x) =  (x+1)^4 - x e^(sin(x)) + 5x -8
    #define this function using variable name f2
    f2 = lambda x: math.pow(x+1,4) - x * math.exp(math.sin(x)) + 5*x - 8
    #add code to call bisection to solve f2(x)=0 on [-3, 3]
    (x, k, fx) = bisection(f2, -3, 3)
    print('\n By bisection: x={},  k={},  fx={}\n'.format(x, k, fx))
    
    #add code to call newton to solve f2(x)=0, set initial root as 3
    #you need to manually get the derivative of f1, pass it as a function named f2der 
    f2der = lambda x: - x * math.exp(math.sin(x)) * math.cos(x) - math.exp(math.sin(x)) + 4 * math.pow(x+1,3) + 5
    (x, k, fx) = newton(f2, f2der, 3) 
    print('\n By Newton: x={},  k={},  fx={}\n'.format(x, k, fx))    

    #add code to call quasi_newton to solve f2(x)=0, set initial root as 3
    (x, k, fx) = quasi_newton(f2, 3, .5)
    print('\n By quasi_Newton: x={},  k={},  fx={}\n'.format(x, k, fx))    

    #f3(x) = e^(-x^2) + x^3 - 100  
    #interestingly, for this function, the default tol lead to trouble for newton's method!
    #so you set tol=1e-13 instead of the default 1e-14
    #define this function using variable name f3
    f3 = lambda x: math.exp(math.pow(-x,2)) + math.pow(x,3) - 100
    #add code to call bisection to solve f3(x)=0 on [-10, 10], call using tol=1e-13
    (x, k, fx) = bisection(f3, -10, 10, tol=1.0e-13)
    print('\n By bisection: x={},  k={},  fx={}\n'.format(x, k, fx))
    
    #add code to call newton to solve f3(x)=0, set initial root as 10, call using tol=1e-13
    #you need to manually get the derivative of f3, pass it as a function named f3der 
    f3der = lambda x: 3 * math.pow(x,2) - 2*x*math.exp(math.pow(-x,2))
    (x, k, fx) = newton(f3, f3der, 3, tol=1.0e-13)
    print('\n By Newton: x={},  k={},  fx={}\n'.format(x, k, fx))    

    #add code to call quasi_newton to solve f3(x)=0, set initial root as 10, 
    #call using tol=1e-13,  h=0.3, and FD='FFD'
    (x, k, fx) = quasi_newton(f2, 3, h=.3, FD='FFD', tol=1.0e-13)
    print('\n By quasi_Newton: x={},  k={},  fx={}\n'.format(x, k, fx))    

    #add code to call quasi_newton to solve f3(x)=0, set initial root as 10, 
    #call using tol=1e-13,  h=0.3, and FD='BFD'
    (x, k, fx) = quasi_newton(f2, 3, h=.3, FD='BFD', tol=1.0e-13)
    print('\n By quasi_Newton: x={},  k={},  fx={}\n'.format(x, k, fx))
    
if __name__=='__main__':
    
    solver_tests()

