from sinP_approx import TrigP
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

P=3
LAMBDA = 1
y0 = 0
x0 = 0
xn = 13
h=0.005
def pr(msg):
    TrigP.pr(msg)



def Qfunc(t):
    t = t%1.3
    return t**2

def ODE(xn, dny0, Lambda):
    y = dny0
    return 1 + ((P/(P-1))-1)*(Lambda - Qfunc(xn) - 1)*np.abs(TrigP.getcosp(y,P))**P




def RK_Method(ODE, x0, dny0, h, xn, Lambda):
    numIteration = round((xn-x0) / h)
    yn = [dny0]
    for i in range(numIteration-1):
        k1 = h * ODE(x0, dny0, Lambda)
        k2 = h * ODE(x0 + 1/2 * h, dny0 + 1/2 * k1, Lambda)
        k3 = h * ODE(x0 + 1/2 * h, dny0 + 1/2 * k2, Lambda)
        k4 = h * ODE(x0 + h, dny0 + k3, Lambda)
        dny0 = dny0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
        yn.append(dny0)
        x0 = x0+h
    return yn

def plotRKM(x0, y0, h, xn,col):
    numIteration = round((xn-x0) / h)
    xAxis = []
    yAxis = []
    for i in range(numIteration) :
        xAxis.append(x0 + i * h)
        yAxis.append(np.cos(x0 + i*h))
     #RK_Method(ODE, x0, y0, h, xn)
    plt.plot(xAxis, yAxis,color=col)
    return yAxis

def rotationNum():
    xaxis = []
    yaxis = []
    for i in range(500):
        xaxis.append(0.01*i)
        theta = RK_Method(ODE, 0, 0, 0.005, 13, xaxis[i])
        yaxis.append((theta[len(theta)-1] - 0)/13)
    plt.plot(xaxis, yaxis)



rotationNum()
plt.show()



