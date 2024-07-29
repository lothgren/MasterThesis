#Non-complete graphs

import scipy.integrate as sp
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import statistics as st


mu = [ 0.6, 0.7, 0.7]
c = [10, 3, 2.1]
g = [10, 5, 3]
tau = 2

def _f(t,x):
    dx = []
    de = 0
    for i in range(len(x)-1):
        dx.append(x[i]*(1-x[i])*(2*x[i] + mu[i]*x[-1] -c[i]))
        de += g[i]*(1-x[i])
    dx.append((de -tau)*x[-1])

    return dx

def solver(x0,N,h,func=_f):
    tSpan = [0,N*h]
    tVec = np.arange(0,N*h,h)
    sol = sp.solve_ivp(func, tSpan, x0, t_eval=tVec ,method="RK45", dense_output=True)
    return sol.t, sol.y



def simCurves():
    N = 1000
    h = 0.01    
    t, y = solver([0.4, 0.4, 0.5],N,h,func=_f)
    for i in range(len(y)-1):
        #plt.figure()
        plt.plot(y[i],y[-1],linewidth=0.5)
        plt.plot(y[i][-1],y[-1][-1],"rx")
        print(y[i][-1])
        
        plt.xlabel("Procentage behaving environmentall friendly, x")
        plt.ylabel("Environmental impact, $\epsilon$")
        plt.title("$\mu =$ {}, $\gamma=$ {}".format(mu[i],g[i]))
    
    plt.plot([0,1],[c[1]/mu[1],c[1]/mu[1]],"k--",linewidth=0.5)
    plt.plot([0,1],[c[0]/mu[0],c[0]/mu[0]],"k--",linewidth=0.5)
    plt.plot([1-tau/g[0],1-tau/g[0]],[0,80],"k--",linewidth=0.5)
    plt.plot([1-tau/g[1],1-tau/g[1]],[0,80],"k--",linewidth=0.5)
    plt.show()



simCurves()