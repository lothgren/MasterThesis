import scipy.integrate as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random as rnd
import statistics as st

a = 0.3
s = 0.6
k = 3
c = -(a+s-k-1) #must be larger than 2??

mu = 0.6
tau = 3
g = 10


def _f(t,x):
    dx1 = x[0]*(1-x[0])*(2*x[0] + mu*x[1] - c)
    dx2 = (g*(1-x[0])-tau)*x[1]
    return [dx1, dx2]

def xLeaves(t,x):
    if x[0] > 0:
        return x[0]-1
    else:
        x[0]
xLeaves.terminal = True

def eLeaves(t,x):
    return x[1]
eLeaves.terminal = True


def solver(x0,T,h=np.inf,func=_f):
    tSpan = [0,T]
    tVec = np.arange(0,T,0.01)
    bd = [eLeaves]
    sol = sp.solve_ivp(func, tSpan, x0, t_eval=tVec ,method="RK45", dense_output=True, max_step = h, events=bd)
    for i in range(len(bd)):
        if len(sol.t_events[i]) != 0:
            print("e has left!")
    return [sol.t, sol.y[0], sol.y[1]]


def genInit(eRange = [0,1]):
    """Generates random initial values, x0 in (0,1), e0 in eRange"""
    x0 = rnd.random()
    e0 = eRange[0] + (eRange[1]-eRange[0])* rnd.random()
    if x0 == 0 or x0 == 1 or e0 == 0:
        raise Exception("Initial values on boundary")
    return [x0, e0]


def simCurves(e0, m):
    """Plots k curves """
    T = 30
    h = 0.01

    for i in range(m):
        t,x,e = solver([(6+i)/10,e0],T,h)
        df = pd.DataFrame(data = {"x": x, "e": e}, index=t)
        sns.lineplot(data = df, x = "x", y="e", sort=False)
    
    plt.plot(1-tau/g,1/mu*(2*tau/g + k -s -a - 1),"r*")
    
    plt.xlabel("Procentage behaving environmentall friendly, x")
    plt.ylabel("Environmental impact, $\epsilon$")
    plt.savefig("Plots/someCurves2")


def distOfMaxE(N):
    """Gives the distribution of max value of epsilon for random initial values (x_0,e_0) in [0,1]^2"""
    eLst = []
    for i in range(N):
        t,x,e = solver(genInit(),500,0.1)
        eLst.append(max(e))
    df = {"Max $\epsilon$": eLst}
    sns.displot(df, x = "Max $\epsilon$", kde = True)
    plt.savefig("Plots/distMaxE")


def testMu():
    global mu
    N, h = 1000, 0.01
    lst = []
    plt.figure()
    for i in range(20):
        mu = (i+1)/10
        l = []
        for j in range(10000):
            t,x,e = solver(genInit(),N,h)
            lst.append([mu, max(e)])
    df = pd.DataFrame(data = np.array(lst), columns=["mu", "e"])
    sns.lineplot(df, x = "mu", y = "e", err_style="band",errorbar="sd")
    plt.xlabel("$\mu$")
    plt.ylabel("Max $\epsilon$")
    plt.savefig("Plots/muSens")


def testTau():
    global tau    
    T, h = 50, 0.01
    tauList = np.arange(0.001,g-0.001,(g-tau)/20)
    lst = []

    for i in tauList:
        tau = i
        for j in range(1000):
            t,x,e = solver(genInit(),T,h)
            lst.append([tau,max(e)])
        print(tau)
    df = pd.DataFrame(data = np.array(lst), columns=["tau","e"])
    sns.lineplot(df, x = "tau", y = "e", err_style="band", errorbar="sd")
    plt.xlabel(r"$\tau$")
    plt.ylabel("$\epsilon$")
    #plt.savefig("Plots/tauSens")


def testG():
    global g    
    N, h = 1000, 0.01
    gList = np.arange(tau+0.001,20,(g-tau)/20)
    lst = []

    for i in gList:
        g = i
        for j in range(10000):
            t,x,e = solver(genInit(),N,h)
            lst.append([g,max(e)])
    df = pd.DataFrame(np.array(lst), columns=["g","e"])
    sns.lineplot(df, x="g", y="e", err_style="band", errorbar="sd")

    plt.xlabel("$\gamma$")
    plt.ylabel("Max $\epsilon$")
    plt.savefig("Plots/gSense")


def testC(N):
    global c
    cList = np.arange(2+0.01,10,0.2)
    lst = []

    for i in cList:
        c = i
        for j in range(N):
            t,x,e = solver(genInit(),1000,0.01)
            lst.append([c,max(e)])
    df = pd.DataFrame(np.array(lst), columns=["c","e"])
    sns.lineplot(df, x="c", y="e", err_style="band", errorbar="sd")

    plt.xlabel("Constant " + r"$c$")
    plt.ylabel("Max $\epsilon$")
    plt.savefig("Plots/cSense")


def testE0():
    N, h = 600, 0.1
    e0List = np.arange(0.001,20,0.05)
    eList = []
    for e0 in e0List:
        l = []
        for i in range(100):
            t,x,e = solver([rnd.random(), e0],N,h)
            l.append(max(e))
        eList.append(st.mean(l))
    
    plt.figure()
    plt.plot(e0List,eList)
    plt.xlabel("$\epsilon_0$")
    plt.ylabel("Max $\epsilon$")



def main():
    sns.set_theme()

    #testE0()
    #simCurves(0.6,1) 
    #distOfMaxE(1000)
    #testMu()
    testTau()
    #testG()
    #testC(1000)

    plt.show()
    #help(sp.solve_ivp)

main()


