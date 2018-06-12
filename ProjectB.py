#Investigating asymmetry in particle collisions

import scipy as sp
import numpy as np
from scipy import *
from numpy import *
from matplotlib import pyplot as plt
import scipy.optimize as spo

#Task 1
print "Task 1"
# PART 1
# plotting the functions of sigmaA and sigmaB as a function of COM energy

# step 1; writing out some predefined parameters to minimize mistypes
k=1000000.0    #constant of the collider     
Mz=91.0        #mass of the Z boson in GeV      
Gz=2.5        #width of Z boson defined as per equation 
rs=0.14      #the parameter r describes the measurements at peak of cross-section for Z production 
js=-0.033    #j parameters describe energy dependance of total and forward and backward cross-sections   
ra=0.0027         
ja=0.81 


# step 2; defining the functions sigmaA and sigmaS with the squared COM 
#energy as its argument
def sigS(S):
    result= ((4/3)*pi)*((1/S)+(((S*rs)+(S-Mz**2)*-js)/((S-Mz**2)**2+((Mz**2)*(Gz**2)))))
    return result
    
def sigA(S):
    result2= (pi)*(((S*ra)+(S-Mz**2)*ja)/((S-Mz**2)**2+((Mz**2)*(Gz**2))))
    return result2
    
#step 3; a linear space of COM energies within the accelerators range (20-140GeV),was created-
e = np.linspace(20, 140, 100)
# squaring all the COM energies to get S values,
s = e*e

#step 4; the function SigmaS was performed on the numerical array S, 
# and the results were labelled B
B = sigS(s)
print "COM energies squared- s values", s
print "SigS(s)", B

#step 5; the function SigmaA was performed on the numerical array S, 
# and the results were labelled C
C= sigA(s)
print "COM energies sqaured- s values", s
print "SigA(s)", C

# step 6; plot the functions SigmaA(s) and SigmaS(s) 
# against the COM energies, e
#plt.plot(e, B, 'r-')
#plt.plot (e, C, 'b-')
#plt.legend(('sigma-S', 'sigma-A'),loc='upper right')
#plt.xlabel("COM energy/ GeV")
#plt.ylabel("Sigma-A and Sigma-S values/ GeV^-2")
#plt.show()
#plt.savefig("task1plot1.png")
#plt.clf()

# PART 2
# Plotting of Total Muon Pairs Produced per Day as 
# a Function of COM Energy

# step 1; the expression for DNu/Dcos(theta)
# was integrated between -0.95 and 0.95 (the range of angles), 
#in order to obtain an expression for the number of muon pairs as a function of COM energy


#step 2; the function for Number of muon pairs 
# with the argument, s (COM energy squared) is defined in python
def Nu(x):
    result=2*1000000*sigS(x)*((0.95)+((0.95**3)/3))
    return result

#step 3; the function was carried out on the predefined s values
D= Nu(s)
print s
print D

# step 4; the aqcuired expression for of Nu as a function of 
# COM energy was plotted.

plt.plot(e, D, 'r-')
plt.xlabel("COM energy / GeV")
plt.ylabel("Number of muon pairs produced")
plt.show()
plt.savefig("task1plot2.png")
plt.clf()


#Task 2
print "Task 2"
#simulating the number of muons detected per cos theta bin for an accelerator 
#operating at a fixed COM energy.
#step 1; fix COM energy as 90.0 GeV, then s=8100 GeV^2
sqrtS =90
sfixed=sqrtS**2

#step 2; calculate sigma A and sigma S parameters for s=8100 GeV^2
SS= sigS(sfixed)
print "Sigma-S values", SS
SA= sigA(sfixed)
print "Sigma-A values", SA

#step 3; set dcos(theta) as 0.01
#separate the data points into cos theta bins of width 0.01
dt= 0.01
cost = np.arange(-0.95, 0.95, dt)

#step 4; define the expression for the Expected value of muon pairs per 
#cos theta bin
#this is DNu/Dcos(theta) 

def dNu(theta):
    result3 =dt*k*(SS*(1+theta**2)+SA*theta)
    return result3
    
#step 5; calculate expected number of muons per cos theta bin
D= dNu(cost)
print "cos theta values", cost
print "expected number of muon pairs", D

#step 6; plot expected number of muons as a function of cos theta 
plt.plot(cost, D, 'r-')
plt.xlabel("cos theta")
plt.ylabel("Expected number of muon pairs")
plt.show()
plt.savefig("task2plot1.png")

#step 7; find the number of data points in D
N= len(D)

#step 8; smear using a Poisson random number generator, 
#with mean D, and the number of trials as the number of D values
E = np.random.poisson(D,N)
print "smeared data with poisson noise", E

#step 9; plot the smeared data as a histogram. 
#The red line is the expected value, D

plt.bar(cost,E,width=1e-2)
plt.plot(cost, D, 'r-')
plt.ylabel("Number of muon pairs")
plt.legend(('expected value', 'smeared data'),loc='upper right')
plt.show()
plt.savefig("task2plot2.png")
plt.clf()

#step 10; calculate the error in the experimental values
#this is a counting experiment following Poisson distribution
#thus the error= square root of counts for each interval

err=(E)**0.5
print "Error", err

#step 11; Add these as error bars to the graph
plt.plot(cost, E, 'b-')
plt.xlabel("cos theta")
plt.ylabel("(Experimental) Number of muon pairs")
plt.errorbar(cost, E, xerr=0, yerr= err)
plt.show()
plt.savefig("task2plot3.png")
plt.clf()

#Task 3
print "Task 3"

#Part 1- plotting the actual ratios of SigmaA/SigmaS for comparison
e=np.linspace(20, 140, 190)
s=e*e
SS= sigS(s)
SA= sigA(s)
ratio1= (SA/SS)

plt.plot(e, ratio1, 'r-')
plt.xlabel("Center of Mass Energies/ GeV")
plt.ylabel("ratio of Sigma-A/Sigma-S")
plt.show()
plt.savefig("task3plot1.png")
plt.clf()

#Part 2- fitting the ratio and plotting as a function of COM energies and bin widths
#step 1; define function for data fitting-
def fitfunc(x,a,s):
    func= k*(s*(1+x**2)+a*x)
    return func
#step 2; create array of ratios
Aratios= np.array([], dtype=float)
Aerr=np.array([], dtype=float)

#step 3;create loop that goes through values of e and fits for the ratio 
for q in e:
    s=q*q
    SS=sigS(s)
    SA= sigA(s)
    D=dNu(cost)
    E=np.random.poisson(D,len(D))
    err=(E)**0.5
    x=np.arange(-0.95, 0.95, dt)
    y= np.array([E], dtype=np.float)
    y_err= err
    initial_guess=[SA,SS]
    po,po_cov= spo.curve_fit(fitfunc,cost,E,initial_guess,y_err)
    ratio2 =po[0]/po[1]
    rerr= sp.sqrt(((po_cov[0,0])/po[0])**2+((po_cov[1,1])/po[1])**2) #adding errors in quadrature to get error in ratio
    Aratios= np.append(Aratios, ratio2, axis=None)
    Aerr= np.append(Aerr, rerr, axis=None)

#step 4;plot
plt.plot(e,Aratios, 'b-')
plt.errorbar(e,Aratios, xerr=0, yerr= Aerr)
plt.xlabel("Center of Mass Energies/ GeV")
plt.ylabel("ratio of Sigma-A/Sigma-S")
plt.show()
plt.savefig("task3plot2.png")
plt.clf()

#step5; effect of bin width on ratio
#fix COM energy
sqrtS= 90.0
s= sqrtS**2 

  
#step 6; calculating simulated data with poisson noise
def dNu(dt,theta):
    result3 =dt*k*(SS*(1+theta**2)+SA*theta)
    return result3

dt=np.array([0.01,0.05,0.1])
for x in dt:
    print "bin width=", x
    cost = np.arange(-0.95, 0.95, x)
    D= dNu(x,cost)
    E=np.random.poisson(D,len(D))
    err=(E)**0.5

#step 7; define function for data fitting-
    def fitfunc(x,a,s):
        func= k*(s*(1+x**2)+a*x)
        return func
    x=np.arange(-0.95, 0.95, x)
    y= np.array([E], dtype=np.float)
    y_err= err

#step 8; make an initial guess and use curve fitting function 
# to guess parameters and their uncertainties
    initial_guess=[-0.004637,0.042374]
    po,po_cov= spo.curve_fit(fitfunc,cost,E,initial_guess,y_err)

    #step 9; print parameters and uncertainties
    print "Values with errors;"
    print "sigma a=", po[0], "+/-", sp.sqrt(po_cov[0,0])
    print "sigma s=", po[1], "+/-", sp.sqrt(po_cov[1,1])

    #step 10; calulate experimental ratio of SigmaA/Sigma S
    ratio2 =po[0]/po[1]
    rerr= sp.sqrt(((po_cov[0,0])/po[0])**2+((po_cov[1,1])/po[1])**2)
    print "ratio from fitted results", ratio2 
    print "uncertainty in ratio", rerr

