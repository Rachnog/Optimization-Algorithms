# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 21:12:42 2014

@author: Alex
"""

import matplotlib.pylab as pl
import numpy as np
import math


def norm(s1):
    return math.sqrt(s1[0]**2 + s1[1]**2)


def fun(dot):
    return (dot[0] - 2)**2 + 3*dot[1]**2

    
def hj():
    x0 = [4, 9]
    start = x0[:]
    deltax = [0.6, 0.8]
    xses = []
    xses.append(x0)
    epsilon1 = epsilon2 = 0.01
    n = 0

    while True:
        print "WE START WITH" + str(x0)
        print "DELTA X: (%8.5f, %8.5f)" %(deltax[0], deltax[1])
        
        x0 = search(x0, deltax) # Считаем икс 1
        xses.append(x0)
        print "X1 ==> fun(%8.5f, %8.5f): %8.5f" % (x0[0], x0[1], fun(x0))
        
        error_point = norm([x0[0] - xses[len(xses)-2][0], x0[1] - xses[len(xses)-2][1]])/norm(x0)
        error_func = abs(fun(x0) - fun(xses[len(xses)-2]))/abs(fun(xses[len(xses)-2]))
        print "||x^2 - x^1|| / ||x^2|| = %8.5f" % error_point
        print "|f(x^2) - f(x^1)| / |x^2| = %8.5f" % error_func
        
        if error_point < epsilon1 and error_func < epsilon2:
            print "DONE"
            print "N=%d" % (n+1)
            break 

        if fun(x0) > fun(xses[len(xses)-2]):
            print "fun(x%d) > fun(x%d) (%8.5f > %8.5f)" %(n+1, n, fun(x0), fun(xses[len(xses)-2]))
            print "REDUCING STEP"
            print "SET x0 TO (%8.5f, %8.5f)" %(xses[len(xses)-3][0], xses[len(xses)-3][1])
            x0 = xses[len(xses)-3]
            start = x0[:]
            deltax[0] /= 2
            deltax[1] /=2
            xses = []
            xses.append(x0)
            n = 0
            print "---------------------------------------------------"
            continue        
        
        xr = [2*x0[0] - start[0], 2*x0[1] - start[1]]
        print "XR fun(%8.5f, %8.5f): %8.5f" % (xr[0], xr[1], fun(xr))
        
        if fun(search(xr[:], deltax)) < fun(x0):
            start = x0[:]
            x0 = search(xr[:], deltax)
            print "X2 ==> fun(%8.5f, %8.5f): %8.5f" % (x0[0], x0[1], fun(x0))
        n+=1
        print "---------------------------------------------------"
                
    plot(xses) 
    for i in range(len(xses)):
        print str(i) + " " + str(xses[i][0]) + ", " + str(xses[i][1]) + " --> " + str(fun(xses[i]))


def search(x0, deltax):
    x1 = [x0[0] + deltax[0], x0[1]]
    if fun(x0) > fun(x1):
        x0 = x1[:]
    x2 = [x0[0] - deltax[0], x0[1]]
    if fun(x0) > fun(x2):
        x0 = x2[:]
    x3 = [x0[0], x0[1] + deltax[1]]
    if fun(x0) > fun(x3):
        x0 = x3[:]
    x4 = [x0[0], x0[1] - deltax[1]]
    if fun(x0) > fun(x4):
        x0 = x4[:]
    return x0    
    
        
def plot(points):
    '''
        Plotting 2D function and way search
    '''
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
    
    pl.plot(xs, ys, marker='o', linestyle='--', color='r', label='Square')            