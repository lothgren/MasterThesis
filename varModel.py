import scipy.integrate as sp
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import statistics as st

a = 0.3
s = 0.6
k = 3
c = 10  #-(a+s-k-1) 

mu = 50
g = 10
t1 = 5.9
t2 = 4.2
eta = 0.3


def tau(e):
    return t2 * e/(1+eta*e) + t1


def fixedPoints():
    try:
        e1 = (g-t1)/(t2+eta*(t1-g))
        x1 = 1
    except:
        e1 = None
        x1 = None
    e2 = 1/(2*mu*eta) * (c*eta-mu-2/g*(eta*t1+t2) + ((c*eta-mu-2/g*(eta*t1+t2))**2 - 4*mu*eta*(2*t1/g-c))**(1/2) )
    x2 = (c-mu*e2)/2

    return [x1, x2], [e1, e2]


def _f(t,x):
    dx1 = x[0]*(1-x[0])*(2*x[0] + mu*x[1] + a+s-k-1)
    dx2 = (g*(1-x[0])-0.1)*x[1]

    return [dx1, dx2]

def _g(t,x):
    dx1 = x[0]*(1-x[0])*(2*x[0] + mu*x[1]-c)
    dx2 = (tau(x[1]) - g*x[0])*x[1]
    return [dx1, dx2]


def solver(x0,N,h,func=_f):
    tSpan = [0,N*h]
    tVec = np.arange(0,N*h,h)
    sol = sp.solve_ivp(func, tSpan, x0, t_eval=tVec ,method="RK45", dense_output=True)
    return [sol.t, sol.y[0], sol.y[1]]


def barrier():
    x = np.arange(0,1.01,0.01)
    e1 = []
    e2 = []
    for i in x:
        e1.append( (g*i-t1)/(t2+eta*(t1-g*i)) )
        e2.append( (c-2*i)/mu )
    plt.plot(x,e1,"k--",linewidth=0.3)
    plt.plot(x,e2,"k--", linewidth=0.3)

def simCurves():
    N = 1000
    h = 0.01
    xFix,efix = fixedPoints()
    plt.figure()
    plt.ylim(-1,efix[0]+1)
    for i in range(1):
        t,x,e = solver([xFix[1]+0.01,efix[1]+0.01],N,h,func=_g)
        plt.plot(x,e,linewidth=0.7)
    barrier()
    #plt.plot(1-tau/g,1/mu*(2*tau/g + (-c) - 1),"rx")
    
    plt.plot(xFix[0],efix[0], "rx")
    plt.plot(xFix[1],efix[1], "bx")
    
    

    plt.xlabel("Procentage behaving environmentall friendly, x")
    plt.ylabel("Environmental impact, $\epsilon$")
    plt.show()


simCurves()