# -*- coding: utf-8 -*-
"""
Created on Sat May 17 16:35:36 2014

@author: Alex
"""

import math
import numpy as np
import matplotlib.pylab as pl


count = 0
r = 1

def incCount():
    global count
    count = count+1

def zeroCount():
    global count
    count = 0


def changeR(value):
    global r
    r = value
    

def calcX(x0, grad, lmb):
    '''
        x0 - gradient*lambda
    '''
    return sub(x0, mults(grad, lmb))

def svenn(x0, grad, lmb, delta):
    """
        One-dimensional Svenn search
    """
    #print "Svenn stage..."
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
    #print "svenn a: " + str(a)
    #print "svenn b: " + str(b)    
    return [a , b]


def dsc(x0, grad, lmb, delta):
    svenn_res = svenn(x0, grad, lmb, delta)
    x1 = svenn_res[0]
    x3 = svenn_res[1]
    x2 = (x1 + x3)/2
    f1 = fun(calcX(x0, grad, x1))
    f2 = fun(calcX(x0, grad, x2))
    f3 = fun(calcX(x0, grad, x3))
    xApprox = x2 + ((x3 - x2) * (f1 - f3)) / (2 * (f1 - 2 * f2 + f3))
    return [x1, x2, x3, xApprox]


def dscPowell(x0, grad, eps, lmb, delta):
    dsc_res = dsc(x0, grad, lmb, delta)
    a = dsc_res[0]
    xmin = dsc_res[1]
    b = dsc_res[2]
    xApprox = dsc_res[3]

    while abs(xmin-xApprox) >= eps or abs(fun(calcX(x0, grad, xmin)) - fun(calcX(x0, grad, xApprox))) >= eps:
        if xApprox < xmin:
            b = xmin
        else:
            a = xmin
        xmin = xApprox
        funcRes =  [fun(calcX(x0, grad, a)), fun(calcX(x0, grad, xmin)), fun(calcX(x0, grad, b))]
        a1 = (funcRes[1] - funcRes[0]) / (xmin - a)
        a2 = ((funcRes[2] - funcRes[0]) / (b - a) - a1) / (b - xmin)
        xApprox = (a + xmin) / 2 - a1 / (2 * a2)
    return xmin          

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

def plot3D(points, col):
    n = 256
    x = np.linspace(-100, 100, n)
    y = np.linspace(-100, 100, n)
    z = np.linspace(-100, 100, n)
    X, Y, Z = np.meshgrid(x, y, z)
    
    xs = []
    ys = []
    zs = []
    
    #pl.contourf(X, Y, Z, fun([X, Y, Z]), 8, alpha=.75, cmap='jet')
    #C = pl.contour(X, Y, Z, fun([X, Y, Z]), 8, colors='black', linewidth=.5) 
    
    for i in range(len(points)):
        xs.append(points[i][0])
        ys.append(points[i][1])
        zs.append(points[i][2])
    
    pl.plot(xs, ys, marker='o', linestyle='--', color=str(col), label='Square')    

def plot(points, col):
    n = 256
    x = np.linspace(-100, 100, n)
    y = np.linspace(-100, 100, n)
    X, Y = np.meshgrid(x, y)
    
    xs = []
    ys = []
    
    pl.contourf(X, Y, fun([X, Y]), 8, alpha=.75, cmap='jet')
    pl.contour(X, Y, fun([X, Y]), 8, colors='black', linewidth=.5) 
    
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
    incCount() 
    return ourFun(x) + constraintFun(x, r)

def ourFun(x):
    """
        Целевая функция для минимизации
    """
    #return (x[0] - 4)**2 + (x[1] - 4)**2
    #return 4*(x[0]-5)**2 + (x[1] - 6)**2 
    #return (1-x[0])**2 + 100*(x[1] - x[0]**2)**2
    return (10*(x[0] - x[1])**2 + (x[0] - 1)**2)**4
    #return 4*(x[0]-4)**2 + x[0]*x[1] + 3*x[1]**2
    

def constraintFun(x, r):
    """
        Часть функция для минимизации, где мы работаем с ограничениями
    """
    line1 = 3*(x[1] - 1)**2 + (x[0] - 1)**2 - 0.5
    line2 = (x[1] - 1)**2
    return r*sqrtCut(line1)**2 + r*sqrtCut(line2)**2 
    #return 0.5*(line + abs(line))**2
    #return r*(1/line)
 
def sqrtCut(x):
    """
        Квадрат срезки
    """
    if x <= 0:
        return x
    else:
        return 0

def conjugatedDirectionPowell(x0, eps, eps2):
    zeroCount()
    iteration = 0
    si = np.eye(len(x0))
    xn = mults(x0, -1)
    xs = []
    xs.append(x0)
    while norm(sub(xn, x0)) > 0.1 * eps:
        iteration+=1
        x0 = xn
        xi = x0
        fi = []

        for i in xrange(len(x0)):
            lmb = dscPowell(xi, si[i], eps2, 0, 0.1)
            xn = calcX(xi, si[i], lmb)
            fi.append(fun(xi) - fun(xn))
            xi = xn
        for i in xrange(len(si)-1):
            si[i] = si[i+1]

        si[len(si)-1] = sub(xn, x0)
        lmb =  dscPowell(xn, si[len(si) - 1], eps2, 0, 0.1)
        xn = calcX(xn, si[len(si) - 1], lmb)
        #print xn
        xs.append(xn)
    #plot(xs, 'red')       
    #print "calculations: " + str(count)
    return xn  

def main():
    point = [0.8,0.7]
    minimum = conjugatedDirectionPowell(point, 0.01, 0.01)   
    for i in range(10):   
        minimum = conjugatedDirectionPowell(point, 0.001, 0.001)           
        print "penalty: " + str(constraintFun(minimum, r) ) + " R = " + str(r)
        changeR(r*10)
        point = minimum
    print point    
    print ourFun(minimum)
    #print "Count " + str(count)
    
    #line = []
    #for i in xrange(-10, 10):
    #   line.append([i, cf([i, i])])
    #plot(line, 'green')    

if __name__ == "__main__":
    main()