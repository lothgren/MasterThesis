#Plotting of model

import numpy as numpy
import matplotlib.pyplot as plt

mu = 0.6
a = 0.3
sig = 0.6
k = 3
tau = 0.1
g = 10

def xDot(x,e):
    """(11) of Frieswijk"""
    dx = x*(1-x)* (2*x+mu*e+a+sig-k-1)
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


def main():
    for i in range(5):
        x,e = eulForward((i+3)/10,0.1,10000,0.001)
        plt.plot(x,e)

    plt.show()

main()
