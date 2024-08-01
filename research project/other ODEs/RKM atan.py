import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

x0 = 0.0
dny0 = np.array([0, 1])
h = 0.0005
xn = 20
def ODE(x, dny0):
    y, dy = dny0
    return np.array([dy, -2*x*dy**2])

def func(x):
    return np.arctan(x)


def RK_Method(ODE, x0, dny0, h, xn):
    numIteration = round(xn / h)
    yn = [dny0[0]]
    for i in range(numIteration-1):
        k1 = h * ODE(x0, np.array(dny0))
        k2 = h * ODE(x0 + 0.5 * h, dny0 + 0.5 * k1)
        k3 = h * ODE(x0 + 0.5 * h, dny0 + 0.5 * k2)
        k4 = h * ODE(x0 + h, dny0 + k3)
        dny0 = dny0 + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4)
        yn.append((dny0[0] - func(x0)) / func(x0))
        x0 = x0 + h
    return yn

# def plotFunc(x0,h,xn):
#     numIteration = round(xn / h)
#     xAxis = []
#     yAxis = []
#     for i in range(numIteration):
#         xAxis.append(x0 + i*h)
#         yAxis.append(func(x0 + i*h))
#     plt.plot(xAxis, yAxis,color="red")


def plotRKM(x0, dny0, h, xn):
    numIteration = round(xn / h)
    xAxis = []
    for i in range(numIteration) :
        xAxis.append(x0 + i * h)
    yAxis = RK_Method(ODE, x0, dny0, h, xn)
    plt.plot(xAxis, yAxis,color='purple',)




plotRKM(x0, dny0, h, xn)
#plotFunc(x0,h,xn)
plt.xlim(left=0)
plt.show()

