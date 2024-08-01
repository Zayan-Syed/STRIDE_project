import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

class TrigP:
    RANGE = 1000
    P=2

    cosp = None
    def pr(msg):
        print('\n')
        print(msg,'\n')

    #Represents the signed p-1 function
    def phi(t, p):
        if t ==0: return 0
        return t*(1.0*np.abs(t))**(p-2)

    #func represents the nth term in the sequence that will be integrated and 
    # each element represnts the value of the function at each splice
    #xaxis represents the range [0,pi_p/2] 
    def innerIntegral(p, func,xaxis):
        for i in range(len(func)):
            #Applies the signed p-1 function array
            func[i] = TrigP.phi(func[i],p)
        #finds the index that of the xaxis array that contains the value of PI_p/2
        index = round((TrigP.PIp(p)/2) / xaxis[1])
        #integrates the function with the constant limits 0 and PI_p/2
        integPIp2 = integrate.simpson(func[:index+1],x=xaxis[:index+1])
        tot=0
        #Stores the integral at each splice
        integ = []
        integ.append(integPIp2)
        for i in range(len(func)-1):
            #Splits the x and y axis splices into a further 100 elements 
            # so it can be integrated using the simpson method
            xrange= np.linspace(xaxis[i],xaxis[i+1],100)
            yrange=np.linspace(func[i],func[i+1],100)
            tot = tot + (integrate.simpson(yrange,x=xrange))
            integ.append(integPIp2 - tot)
        return integ



        
    def outerIntegral(p, func, xaxis):
        q = p/(p-1)
        for i in range(len(func)):
            func[i] = TrigP.phi(func[i], q)
        tot=0
        integ = []
        #Stores the value at x=0 manually as the value at 0 is known
        integ.append(0)
        for i in range(len(func)-1):
            xrange=np.linspace(xaxis[i],xaxis[i+1],100)
            yrange=np.linspace(func[i],func[i+1],100)
            tot = tot + (integrate.simpson(yrange,x=xrange))
            integ.append(tot)
        return integ
        


    def gensinp(p):
        xaxis = []
        #Splites the region 0,PI_p/2 into 100 pieces
        interval = (TrigP.PIp(p)/2)/TrigP.RANGE
        for i in range(TrigP.RANGE+1):
            #Stores the values x1 to xn
            xaxis.append(i*interval)
        #Creates the initial 0th element of the function seqeunce
        yaxis = [1 for i in range(TrigP.RANGE+1)]
        for i in range(30):
            yaxis = TrigP.innerIntegral(p, yaxis,xaxis)
            yaxis = TrigP.outerIntegral(p,yaxis,xaxis)
        #Takes the value of the sequence function at PI_p/2
        MAX = yaxis[TrigP.RANGE]
        P_ROOT_CONSTANT = (p-1)**(1.0/p)
        for i in range(len(yaxis)):
            yaxis[i] = yaxis[i] / MAX
            yaxis[i] = yaxis[i] * P_ROOT_CONSTANT
        #Reflects the function in the line x=PI_p/2 to get the function in the region [0,PI_p]
        for i in range(1,TrigP.RANGE+1):
            yaxis.append(yaxis[TrigP.RANGE-i])
            xaxis.append(TrigP.PIp(p)/2 + i*interval)
        for i in range(1,len(yaxis)):
            yaxis.append(-1*yaxis[i])
            xaxis.append(TrigP.PIp(p)+i*interval)
        return (xaxis,yaxis)

    def gencosp(p):
        sinp = TrigP.gensinp(p)
        sinpx = sinp[0]
        sinpy = sinp[1]
        pi_indices = []
        pi_indices.append(0) #[0]
        pi_indices.append(round((TrigP.PIp(p)/2) / sinpx[1])) #[1]
        pi_indices.append(round(TrigP.PIp(p) / sinpx[1])) #[2]
        pi_indices.append(round((TrigP.PIp(p) * 3/2) / sinpx[1])) #[3]
        pi_indices.append(round((2*TrigP.PIp(p)) / sinpx[1])) #[4]
        cosp = []
        cosp.extend(sinpy[pi_indices[1]:pi_indices[2]])
        cosp.extend(sinpy[pi_indices[2]:pi_indices[3]])
        cosp.extend(sinpy[pi_indices[3]:pi_indices[4]])
        cosp.extend(sinpy[pi_indices[0]:pi_indices[1]+1])
        return (sinpx, cosp)

    def plot_trip_func(axes):
        plt.plot(axes[0], axes[1], color="purple")


    def PIp(p):
        return (2 * (p-1)**(1/p) * np.pi / p)  / np.sin(np.pi / p)
    
    def getcosp(x,p):
        if p != TrigP.P or TrigP.cosp == None:
            TrigP.cosp = TrigP.gencosp(p)
            TrigP.P = p
        cospx = TrigP.cosp[0]
        cospy = TrigP.cosp[1]
        #Represents 2* PI_p
        PIP_2 = cospx[len(cospx)-1]
        if x > PIP_2:
            x = x%PIP_2
        lowbound = round(x / cospx[1])
        if (lowbound > x):
            lowbound = lowbound-1
        upbound = lowbound + 1
        if upbound >= len(cospx):
            return cospy[len(cospy)-1]
        gradient = (cospy[upbound] - cospy[lowbound]) / (cospx[upbound] - cospx[lowbound])
        return cospy[lowbound] + (x - cospx[lowbound])*gradient 
    
