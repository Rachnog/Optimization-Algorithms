# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 19:25:20 2014

Rosenbrock's optimization algorithm with visualization
Some helpful outputs for better intuition are commented ('#')
 
@author: Alex
"""
 
import numpy as np
import math
import matplotlib.pylab as pl

count = 0
    
def norm(s1):
    normas = 0
    for i in xrange(len(s1)):
        normas += s1[i]**2
    return math.sqrt(normas)      

def incCount():
    global count
    count = count+1

def perpendicular(a):
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b


def fun(x):
    incCount()
    #return 4*(x[0]-5)**2 + (x[1]-6)**2
    #return (x[0] - variant)**2 - x[0]*x[1] + 3*x[1]**2
    #return (x[0] - 3)**2 + x[0] * x[1] + 3 *x[1]**2
    #return (1-x[0])**2 + 100*(x[1] - x[0]**2)**2  
    return (10*(x[0] - x[1])**2 + (x[0] - 1)**2)**4

def constraintFunCircle(x):
    '''
        Все ограничения подаем в виде >=0 и меняем знаки
    '''
    return -(x[0]+1.2)**2 - (x[1])**2 + 1    


def constraintFunLine(x):
    """
        Не меняем знак???
    """
    return 1 - x[0] - 2*x[1]


def constraintFunFignya(x):
    return -(x[1]**3 + x[0]**2-3*x[1])**4 - (x[0]**3-x[1]**2-3*x[0]) - x[1]**2 + 3

def outCircle(x):
    return -(x[0]-5)**2 - (x[1]-5)**2 + 10

def innerCircle(x):
    return 2*(x[0]+0.5)**2 + (x[1]+0.7)**2 - 1

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
     
 
def rozenbrock(x0, lmb, epsilon1, epsilon2, alpha, beta):
    """
        x0 - start point, define as [x, y]
        lmb - delta x, define as [0.1, 0.1]
        epsilon1, epsilon2 - estimates errors
        alpha, beta - coefficients for increasing\decreasing step
        success = array for saving YES\NO answers
    """

    start = x0
    s1 = [1, 0]  # s1 vector
    s2 = [0, 1]  #s2 vector
    N = 0  # iterations count
    success1 = []
    success2 = []
    points = []
    points.append(x0)
    iterations = 0
    
    temp = []

    while True:
        temp = x0        
        x1 = [x0[0] + lmb[0]*s1[0], x0[1] + lmb[0]*s1[1]]

        if fun(x1) < fun(x0):
            success1.append("Y")
            lmb[0] *= alpha
            x0 = x1[:]
 
        elif fun(x1) > fun(x0):
            success1.append("N")
            lmb[0] *= beta
 
        x2 = [x0[0] + lmb[1]*s2[0], x0[1] + lmb[1]*s2[1]]
        
       
        if fun(x2) < fun(x0):
            success2.append("Y")
            lmb[1] *= alpha
            x0 = x2[:]
        elif fun(x2) > fun(x0):
            success2.append("N")
            lmb[1] *= beta
        
        if constraintFunCircle(x0) >= 0:
            lmb[0] /= 2
            lmb[1] /= 2
            print "lmb reduced cause of constraint"
            print lmb
            x0 = temp[:]
            continue
 
        if success1[len(success1)-1] == "N" and success2[len(success2)-1] == "N" and "Y" in success1 and "Y" in success2:
            print "Our stop point: [%8.5f, %8.5f]" % (x0[0], x0[1])
            s1 = [x0[0] - start[0], x0[1] - start[1]]
            print "S1" + str(s1)
            norma = norm(s1)
            s1[0] /= norma
            s1[1] /= norma
            s2 = perpendicular(s1)
            print "F(X1) = %8.5f" % fun(x0)
            
            error_point = norm([x0[0] - start[0], x0[1] - start[1]])/norm(x0)
            error_func = abs(fun(x0) - fun(start))/abs(fun(start))

            print "||x^2 - x^1|| / ||x^2|| = %8.5f" % error_point
            print "|f(x^2) - f(x^1)| / |x^2| = %8.5f" % error_func

            if error_point < epsilon1 or error_func < epsilon2 and N>1:
                print "DONE"
                print "N = %3d" % N
                break
            
            points.append(x0)
            start = x0[:]
            success1 = []
            success2 = []
            N+=1
                    
            print "___________________________________________________________________________"
        iterations+=1 
    
      
    #xOne = dsc_pauell(x0, 0.1, epsilon2)
    #print xOne    
    plot(points)   
    print "FUNCTIONS COUNT = " + str(count)
    print "RESULT " + str(x0)


def main():
    x0 = [4, 4]
    rozenbrock(x0, [0.3, 0.1], 0.001, 0.001, 2.9, -0.5)    
    
    
if __name__ == '__main__':
   main()  

   