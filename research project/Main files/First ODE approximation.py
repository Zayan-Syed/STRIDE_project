import matplotlib.axes
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import matplotlib

x0 = 0
#Initial value of the original function
y0 = 1
#Initial value of the derivative of the function
dy0 = 0
#Step factor
h = 0.00001
#input for the final output that we are looking for
xn = 2 * np.pi

#Colors used for plotting alternate graphs
EULER_COL= 'red'
RK_COL = 'blue'
COS_COL = 'green'
SCI_RK_COL = 'yellow'
RKB_COL = 'purple'

#Differential equation that maps to first derivative
def dydx(x, y, dy):
    return dy
#Differential equation that maps to second derivative
def d2ydx2(x, y, dy):
    return -1 * y

#Takes in an input vector [y,dy] and returns a vector of differential equations 
def ODE(x, dny0):
    #y represents the expression y', dy represents dy' or y''
    y, dy = dny0
    #Returns the differential equations: y' = dy and y'' = -y
    return np.array([dy, -1*y])

#implementation of Euler method
def Euler_method(x0, y0, dy0, h, xn):
    numIteration = round(xn / h)
    yn = [y0]
    for i in range(numIteration-1):
        #yn+1 term
        y1 = y0 + h * dydx(x0,y0, dy0)
        #y'n+1 term
        dy1 = dy0 + h * d2ydx2(x0, y0, dy0)
        y0 = y1
        dy0 = dy1
        #Takes error by doing y-cos/cos
        yn.append((y1 - np.cos(x0 + i*h)) / np.cos(x0 + i*h))
    return yn

def plotEulerM(x0, y0, dy0, h, xn):
    numIteration = round(xn / h)
    xAxis = []
    for i in range(numIteration) :
        xAxis.append(x0 + i * h)
    yAxis = Euler_method(x0, y0, dy0, h, xn)
    plt.plot(xAxis, yAxis,color=EULER_COL)


def RK_Method(ODE, x0, dny0, h, xn):
    numIteration = round(xn / h)
    yn = [dny0[0]]
    for i in range(numIteration-1):
        k1 = h * ODE(x0, np.array(dny0))
        #dny0 is a numpy vector property so the expression below applies additiion 
        # and scalar multiplication to all elements
        k2 = h * ODE(x0 + 1/2 * h, dny0 + 1/2 * k1)
        k3 = h * ODE(x0 + 1/2 * h, dny0 + 1/2 * k2)
        k4 = h * ODE(x0 + h, dny0 + k3)
        dny0 = dny0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
        #The 0th element contains the value of y at the given point so we take that value
        yn.append(dny0[0])
        x0 = x0+h
    return yn



def plotRKM(x0, y0, dy0, h, xn):
    numIteration = round(xn / h)
    xAxis = []
    for i in range(numIteration) :
        xAxis.append(x0 + i * h)
    yAxis = RK_Method(ODE, x0, np.array([y0, dy0]), h, xn)
    for i in range(len(yAxis)):
        if np.cos(x0+i*h) == 0:
            yAxis[i] = None
            #100*np.abs(yAxis[i] - np.cos(x0+i*h)) / 0.00000001
        else:
            yAxis[i] = 100*np.abs(yAxis[i] - np.cos(x0+i*h)) / np.abs(np.cos(x0+i*h))
            yAxis[i] = np.log(yAxis[i])
    plt.plot(xAxis, yAxis,color=RKB_COL)
    plt.xlabel("x")
    plt.ylabel("Log of Error taken compared to cos x")


def plotCos(x0,h,xn):
    numIteration = round(xn / h)
    xAxis = []
    yAxis = []
    for i in range(numIteration):
        xAxis.append(x0 + i*h)
        yAxis.append(np.cos(x0 + i*h))
    plt.plot(xAxis, yAxis,color=COS_COL)


def plotSciRK(x0, y0, dy0, h, xn) :
    t = np.linspace(x0, xn, round(xn / h))
    sol = integrate.solve_ivp(ODE, (0,max(t)), [y0, dy0], 'RK45', t, False, None, True, None)
    plt.plot(t, sol.y[0].tolist(), color=SCI_RK_COL)



plotRKM(x0, y0, dy0, h, xn)
plt.xlim(left=0)
plt.show()




    



