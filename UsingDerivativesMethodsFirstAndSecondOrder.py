# -*- coding: utf-8 -*-
"""
Created on Sun Apr 06 14:01:22 2014

@author: Alex
"""
import math
import numpy as np
import matplotlib.pylab as pl

count = 0

def calcX(x0, grad, lmb):
    '''
        x0 - gradient*lambda
    '''
    return sub(x0, mults(grad, lmb))

def svenn(x0, grad, lmb, delta):
    """
        One-dimensional Svenn search
    """
    print "Svenn stage..."
    f0 = fun(calcX(x0, grad, lmb))
    if f0 < fun(calcX(x0, grad, lmb+delta)):
        delta = -delta
    x1 = lmb + delta
    f1 = fun(calcX(x0, grad, x1))
    while f1 < f0:
        delta *= 2
        lmb = x1
        x1 = lmb + delta
        f0 = f1
        f1 = fun(calcX(x0, grad, x1))
    a = lmb + delta/2
    b = lmb - delta/2        
    if a > b:
        temp = b
        b = a
        a = temp     
    print "svenn a: " + str(a)
    print "svenn b: " + str(b)    
    return [a , b]
        

def gold(a, b, eps, x0, grad):
    """
        One-dimensional gold search
    """
    l = b - a
    x1 = a + 0.382*l
    x2 = a + 0.618*l
    while l > eps:
        if fun(calcX(x0, grad, x1)) < fun(calcX(x0, grad, x2)):
            b = x2
            x2 = x1
            l = b - a
            x1 = a + 0.382*l
        else:
            a = x1
            x1 = x2
            l = b - a
            x2 = a + 0.618*l
    print "gold a: " + str(a)
    print "gold b: " + str(b)    
    return [a, b]            
 
           
def calcLambda(x0, grad, eps, lmb):
    line = svenn(x0, grad, lmb, 0.1)
    line = gold(line[0], line[1], eps, x0, grad)
    lmb = (line[0] + line[1])/2
    return lmb    
    

def plot(points, col):
    n = 256
    x = np.linspace(-12, 12, n)
    y = np.linspace(-12, 12, n)
    X, Y = np.meshgrid(x, y)
    
    xs = []
    ys = []
    
    pl.contourf(X, Y, fun([X, Y]), 8, alpha=.75, cmap='jet')
    C = pl.contour(X, Y, fun([X, Y]), 8, colors='black', linewidth=.5) 
    
    for i in range(len(points)):
        xs.append(points[i][0])
        ys.append(points[i][1])
    
    pl.plot(xs, ys, marker='o', linestyle='--', color=str(col), label='Square')

def add(x, y):
    res = []
    for i in xrange(len(x)):
        res.append(x[i] + y[i])
    return res

def sub(x, y):
    res = []
    for i in xrange(len(x)):
        res.append(x[i] - y[i])
    return res        

def mults(x, n):
    res = []
    for i in xrange(len(x)):
        res.append(x[i]*n)
    return res
    
def derivative(x, n):
    h = []
    for i in xrange(len(x)):
        if i == n:
            h.append(0.000000000001)
        else:
            h.append(0)
    return (fun([x[0] + h[0], x[1] + h[1]]) - fun([x[0] - h[0], x[1] - h[1]]))/(2*h[n])

def derivative2(x, a, b):
    ai = []
    aj = []
    for i in xrange(len(x)):
        if i == a:
            ai.append(0.001)
        else:
            ai.append(0)
    for j in xrange(len(x)):
        if j == b:
            aj.append(0.001)
        else:
            aj.append(0)
    return (fun(add(x, add(ai, aj))) - fun(add(x, ai)) - fun(add(x, aj)) + fun(x))/ (ai[a]**2)     

def gradient(x):
    grad = []
    for i in xrange(len(x)):
        grad.append(derivative(x, i))
    return grad 
    
def norm(s1):
    normas = 0
    for i in xrange(len(s1)):
        normas += s1[i]**2
    return math.sqrt(normas)    

def hesse(x):
    h = []
    for i in xrange(len(x)):
        for j in xrange(len(x)):
            h.append(derivative2(x, i, j))
    return h        
    
def fun(x):
    return (x[0] - 6)**4 - x[0]*x[1] + 3*x[1]**4
    #return 4*(x[0]-5)**2 + (x[1] - 6)**2  
    #return (1-x[0])**2 + 100*(x[1] - x[0]**2)**2

def fixedGradient(x0, eps1, eps2):
    """
        Just gradient descent with fixed lamba
        Lambda will be reduced, if next iteration gives worse (bigger) value
    """
    print "fixed gradient goes!"
    xs = []
    xs.append(x0)
    x0 = mults(x0, -1)
    x1 = [-x0[0], -x0[1]]
    lmb = 1.0
    while norm(sub(x1, x0))/norm(x0) >= eps2:
        x0 = x1
        grad = gradient(x0)
        if norm(grad) < eps1:
            break;
        x1 = sub(x0, mults(mults(grad, (1 / norm(grad))), lmb))        
        if fun(x1) > fun(x0):
            lmb /= 2.0
            print "lmb reduced = " + str(lmb)
        print x1
        xs.append(x1)
        plot(xs, 'r') 
    return x1     


def fastestDescent(x0, eps1, eps2):
    """
        Gradient descent with lambda. calculated with one-deimensional search
    """
    print "fastest descent goes!"
    xs = []
    xs.append(x0)
    x0 = mults(x0, -1)
    x1 = [-x0[0], -x0[1]]
    lmb = 0.1
    while norm(sub(x1, x0))/norm(x0) >= eps1:
        x0 = x1
        grad = gradient(x0)
        if norm(grad) < eps1:
            break
        lmb = calcLambda(x0, grad, eps2, lmb)
        x1 = calcX(x0, grad, lmb)
        print x1
        xs.append(x1)
    print "end fastest descent"
    plot(xs, 'g')
    return x1    


def partan(x0, eps1, eps2):
    """
        Parallel tangents method (advanced)
    """
    print "SPartan goes!"
    grad = gradient(x0)
    lmb = 0.1
    xs = []
    xs.append(x0)
    while (norm(grad) >= eps1):
        print "x0 = " + str(x0)
        lmb = calcLambda(x0, grad, eps2, lmb)
        x1 = calcX(x0, grad, lmb)
        print "x1 = " + str(x1)
        grad = gradient(x1)
        lmb = calcLambda(x1, grad, eps2, lmb)
        x2 = calcX(x1, grad, lmb)
        print "x2 = " + str(x2)
        grad = sub(x2, x0)
        lmb = calcLambda(x2, grad, eps2, lmb)
        x2 = calcX(x2, grad, lmb)
        print "x3 = " + str(x2)
        grad = gradient(x2)
        lmb = calcLambda(x2, grad, eps2, lmb)
        x0 = calcX(x2, grad, lmb)
        print "x4 = " + str(x0)
        grad = sub(x0, x1)
        lmb = calcLambda(x0, grad, eps2, lmb)
        x0 = calcX(x0, grad, lmb)
        print "x5 = " + str(x0)
        grad = gradient(x0)
        print "-------------------------"
        xs.append(x0)    
    plot(xs, 'orange')    
    return x0    

def newton(x0, eps):
    """
        Newton second-order optimization method
    """
    xs = []
    xs.append(x0)
    lmb = 0.1
    while norm(gradient(x0)) > eps:
        hesseMatrix = hesse(x0)
        hesseMatrix = np.reshape(hesseMatrix, (-1, len(x0)))
        print "hesse: " 
        print hesseMatrix
        grad = gradient(x0)
        print "grad" + str(grad)
        invHesse = np.linalg.pinv(np.array(hesseMatrix))
        print "hesse inv " 
        print invHesse
        transGrad = np.transpose(np.array([grad])) 
        print "transpose grad"
        print transGrad
        res = np.dot(invHesse, transGrad)
        print "res "
        print res
        normNow = norm(np.transpose(res)[0])
        print "norm" 
        print normNow
        fin = mults(np.transpose(res)[0], 1/normNow)
        print "fin "
        print fin
        lmb = calcLambda(x0, fin, eps, lmb)
        print "lambda " + str(lmb)
        x0 = calcX(x0, fin, lmb)
        print "!!!!" + str(x0)
        xs.append(x0)
    plot(xs, 'yellow')    

def main():
    point = [-3,3]
    fixedGradient(point, 0.01, 0.01)
    fastestDescent(point, 0.001, 0.1)
    partan(point, 0.01, 0.01)
    newton(point, 0.01)


if __name__ == '__main__':
   main()         
        