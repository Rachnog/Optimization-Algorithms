'''
    Little pack for 1-dimensional optimization algorithms
    NOTE: Svenn algorithm is a first step of DSC-algoritm
'''

__author__ = 'Alex Gonchar'

import numpy as np
import matplotlib.pylab as plt
from sympy import *


def fun(x):
    '''
    Input your symbolic function here
    '''
    #return x*(2*x - 3)
    return x**2 + 2/x
    #return (100-x)**2
    #return x**2 - 6*x


def deriv(point):
    '''
        Computes derivative in a point
    '''    
    x = Symbol('x')
    y = fun(x)
    #print str(y)
    d = y.diff()
    #print str(d)
    f = lambdify(x, d, 'numpy')
    return f(point)
 
    
def deriv2(point):
    x = Symbol('x')
    y = deriv(x)
    d = y.diff()
    #print str(d)
    f = lambdify(x, d, 'numpy')
    return f(point)
    
    
def plotFun():
    '''
    Helpful function for plotting your function
    '''  
    x = np.arange(-5,5)
    plt.plot(x, fun(x))


def svenn(x0, delta):
    print "Starting Svenn algorithm..."
    x1 = x0
    x2 = x1+delta
    if fun(x1) < fun(x2):
        delta = -delta
        x1 = x0
        x2 = x1+delta
    while fun(x2) < fun(x1):
        delta*=2
        x1 = x2
        x2 = x1+delta
        print "x1: " + str(x1) + " , f(x1): " + str(fun(x1))
        print "x2: " + str(x2) + " , f(x2): " + str(fun(x2)) 
    a = x2 - 3*delta/2
    b = x2 - delta/2
    if a<b:
        print "a: " + str(a)
        print "b: " + str(b)
    else:
        temp = b
        b = a
        a = temp
        print "a: " + str(a)
        print "b: " + str(b)
    return [a , b]    
    

def dsc(x0, delta):    
    svenn_res = svenn(x0, delta)
    print "Svenn stage done, quadratic approximation once now..."
    a = svenn_res[0]
    b = svenn_res[1]
    
    x1 = a
    x3 = b
    x2 = (abs(b) + abs(a))/2
    
    print "x1 = %8.5f, x2 = %8.5f, x3 = %8.5f" % (x1, x2, x3)
    
    xs = x2 + ((x3-x2)*(fun(x1) - fun(x3)))/(2*(fun(x1) - 2*fun(x2) + fun(x3)))
    print "x* =  %8.5f" % xs
    return [a, b, (a+b)/2, xs]
    
def dsc_pauell(x0, delta, eps):
    dsc_res = dsc(x0, delta)
    print "Pauell stage..."
    a = dsc_res[0]
    b = dsc_res[1]
    xmin = dsc_res[2]
    xs = dsc_res[3]
    
    while abs(xmin - xs) >= eps and abs(fun(xmin) - fun(xs)) >= eps:
        if xs < xmin:
            b = xmin
            xmin = xs
        elif xs > xmin:
            a = xmin
            xmin = xs
        else:
            break
        
        f = [fun(a), fun(xmin), fun(b)]
        print "a = %8.5f, min = %8.5f b = %8.5f" % (a, xmin, b)
        a1 = (f[1] - f[0]) / (xmin - a)
        a2 = ((f[2] - f[0]) / (b - a) - a1) / (b - xmin)
        
        xs = (a + xmin)/2 - a1/(2*a2)
        print "x* = %8.5f" % xs
        
        f_min = min(f)
        if f[0] == f_min:
            xmin = a
        elif f[1] == f_min:
            xmin = xmin
        else:
            xmin = b
        print "a = %8.5f, min = %8.5f b = %8.5f" % (a, xmin, b)    
        
    
def dix(a, b, delta):
    print "Starting dichotomy method..."
    while round(abs(b-a), 3) > abs(delta):
        x = (a + b - delta)/2
        y = (a + b + delta)/2 
        print "i x: %8.5f, y: %8.5f,  L: %8.5f" % (x, y, abs(b-a))
        if fun(x) < fun(y):
            b = y
        else:
            a = x 
    print "[%8.5f; %8.5f] = %8.5f" % (a, b, (a+b)/2)        


def gold(a, b, delta):
    print "Starting gold search..."
    l = b - a
    x = a + 0.382*l
    y = a + 0.618*l
    while abs(b-a) >= delta:
        if (fun(x) < fun(y)):
            print "fun(x1) in %8.5f: %8.5f < fun(x2) in %8.5f: %8.5f" % (x, fun(x), y, fun(y))
            b = y
            y = x
            x = a + b - y
            print "a: %8.5f; b: %8.5f; L: %8.5f" % (a, b, (b-a))
        else:
            print "fun(x1) in %8.5f: %8.5f > fun(x2) in %8.5f: %8.5f" % (x, fun(x), y, fun(y))
            a = x
            x = y
            y = a + b - x
            print "a: %8.5f; b: %8.5f; L: %8.5f" % (a, b, (b-a))
    print "[%8.5f; %8.5f] => %8.5f" % (a, b, (a+b)/2)  
               
            
def bolzano(a, b, delta):
    print "Starting Bolzano method..."
    c = (a+b)/2
    print "a = %8.5f, b = %8.5f, c = %8.5f, f'(%8.5f) = %8.5f" % (a, b, c, c, deriv(c))
    while (abs(deriv(c)) > delta):
        if deriv(c) < 0:
            a = c
            c = (a+b)/2    
            print "a = %8.5f, b = %8.5f, c = %8.5f, f'(%8.5f) = %8.5f" % (a, b, c, c, deriv(c))
        else:
            b = c
            c = (a+b)/2    
            print "a = %8.5f, b = %8.5f, c = %8.5f, f'(%8.5f) = %8.5f" % (a, b, c, c, deriv(c))    
    print "Answer: a = %8.5f, b = %8.5f" % (a, b)


def newton(a, delta):
    xn = a
    xn1 = xn - deriv(xn)/deriv2(xn)
    while abs(deriv(xn)) >= delta:
        xn = xn1
        xn1 = xn - deriv(xn)/deriv2(xn)
        print "xk = %8.5f, f'(xk) = %8.5f" % (xn, deriv(xn))          


def chords(a, b, eps):
    '''
        thanks to @FlyingPirate
    '''
    xApprox = b - ((deriv(b) * (b - a)) / (deriv(b) - deriv(a)))
    while deriv(xApprox) > eps:
        if signum(deriv(xApprox)) * signum(deriv(a)) > 0:
            a = xApprox            
        else:
            b = xApprox
        print "[%8.5f, %8.5f]"  % (a, b)   
        xApprox = b - ((deriv(b) * (b - a)) / (deriv(b) - deriv(a)))        
    return xApprox
    

def signum(x): 
    if x < 0: 
        return -1 
    elif x == 0: 
        return 0 
    else: 
        return 1 

        
    
      
    
            
        
    
    
    
    
            
            
        
    
    
    
           

        



