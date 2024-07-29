#Plotting of model

import numpy as np
import matplotlib.pyplot as plt

a = 0.3
sig = 0.6
k = 3
ask = a+sig-k #must be less than 1

mu = 0.6
tau = 0.1
g = 10

def xDot(x,e):
    """(11) of Frieswijk"""
    if ask >= 1:
        raise ValueError
    dx = x*(1-x)* (2*x+mu*e+ask-1)
    return dx

def eDot(x,e):
    """(11) of Frieswijk"""
    de = e*(g*(1-x)-tau)
    return de

def eulForward(x0,e0,N,h):
    """Does standar euler forward on intial values (x0,e0), N steps with step length h"""
    x = [x0]
    e = [e0]
    for i in range(N):
        x.append(x[i]+h*xDot(x[i],e[i]))
        e.append(e[i]+h*eDot(x[i],e[i]))
    
    return (x,e)


def rk4(x0,e0,N,h):
    x = [x0]
    e = [e0]
    for i in range(N):
        k1 = (xDot(x[i],e[i]),eDot(x[i],e[i]))
        k2 = (xDot(x[i]+h/2*k1[0],e[i]+h/2*k1[1]), eDot(x[i]+h/2*k1[0],e[i]+h/2*k1[1]) )
        k3 = (xDot(x[i]+h/2*k2[0],e[i]+h/2*k2[1]), eDot(x[i]+h/2*k2[0],e[i]+h/2*k2[1]) )
        k4 = (xDot(x[i]+h*k3[0],e[i]+h*k3[1]), eDot(x[i]+h*k3[0],e[i]+h*k3[1]) )
        x.append( x[i] + h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]) )
        e.append( e[i] + h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]) )

    return (x,e)


def simCurves(N,h):
    for i in range(1):
        x,e = rk4((5+i*2)/10,0.01,N,h)
        plt.plot(x,e)
        print((x[-1],e[-1]))

    plt.plot(1-tau/g,1/mu*(2*tau/g + (-ask) - 1),"rx")
    
    plt.xlabel("Procentage behaving environmentall friendly, x")
    plt.ylabel("environmental impact, $\epsilon$")


def testCons(N,h):
    """Testing the sensitivity on a + s - k constant"""
    global ask
    t = np.arange(0,N*h+h,h)
    cons = np.arange(-10,0,1)
    for i in cons:
        ask = i
        x,e=rk4(0.5,0.5,N,h)
        plt.plot(t,e, linewidth=.5)

    plt.xlabel("time, t")
    plt.ylabel("environmental impact, $\epsilon$")


def testMu(N,h):
    """Testing the sensitivity on a + s - k constant"""
    global mu
    t = np.arange(0,N*h+h,h)
    cons = np.arange(0.1,10.1,1)
    print(cons)
    for i in cons:
        mu = i
        x,e=rk4(0.5,0.5,N,h)
        plt.plot(t,e, linewidth=.5)

    plt.xlabel("time, t")
    plt.ylabel("environmental impact, $\epsilon$")


def main():
    N = 1000000
    h = 0.0001
    
    plt.figure()
    simCurves(N,h)
    #testCons(N,h)
    #testMu(N,h)
    plt.show()

main()

