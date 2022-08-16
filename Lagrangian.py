import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.optimize
from sympy import *

"""
NOTE : This code works for both Equal and Unequal Stages. But for 2 or more Stages and Unequal Stages requires really large amount of computation and hence time consuming (Few Hours).

Test Values : Take 2 Equal Stages : Structural Efficiency of 0.015, Isp of 240 s and Burnout Velocity of 4000 m/s or Payload Fraction of 0.176
"""

g0 = 9.81

y = Symbol('y')

print("See the Code first to enter Test/Trial values.")
print("WARNING : DO NOT RUN THIS CODE FOR MULTIPLE UNEQUAL STAGES, IT WILL TAKE TOO MUCH COMPUTATIONAL TIME")

N = int(input("\nEnter No. of Stages (2,3,4 or 5 stages): "))
mL = float(input("\nEnter Payload Mass in Kgs: "))

e = np.linspace(0.1,0.4,N)         #Structural Efficiency of Each Stage
Isp = np.linspace(0.1,0.4,N)       #Isp of Each Stage
Vstage = np.linspace(0.1,0.4,N)    #Optimized Velocity Increment of Each Stage
Pistage = np.linspace(0.1,0.4,N)   #Optimized Stage payload Fraction

for i in np.arange(0,N):
    print("\n******STAGE",i+1,"******")
    e[i] = float(input("Enter Structural Efficiency: "))
    Isp[i] = float(input("Enter Isp: "))


flag = input("\nEnter 0 to Maximize Final Burnout Velocity for a given Mission Payload Fraction as Constraint OR Enter anything else to Maximize Mission Payload Fraction for a given Final Burnout Velocity as Constraint: ")
if flag == '0' :
    Pistar = float(input("\nEnter Mission Payload Fraction: "))
    z=1
    exp = 0
    for i in np.arange(0,N):
        z = z*(-1)*y*e[i]/((1-e[i])*(y+g0*Isp[i]))

    exp = z - Pistar
    rootss = solve(exp, y)
    desired_root=0
    for lamb in rootss:
        for i in np.arange(0,N):
            if (-1)*lamb*e[i]/((1-e[i])*(lamb+g0*Isp[i])) < 0 :
                continue
            else:
                desired_root = lamb
                break


    Vstar = 0
    for i in np.arange(0,N):
        Pistage[i]=(-1)*desired_root*e[i]/((1-e[i])*(desired_root+g0*Isp[i]))
        Vstage[i] = (-1)*g0*Isp[i]*math.log(e[i] + (1-e[i])*Pistage[i])
        Vstar = Vstar + Vstage[i]

    print("Maximum Burnout Velocity: ",Vstar)


else :
    Vstar = float(input("\nEnter Final Burnout Velocity: "))
    z=0
    exp = 0
    for i in np.arange(0,N):
        z = z + (-1)*(g0)*Isp[i]*log((e[i]*y*g0*Isp[i])/(1+y*g0*Isp[i]))
    exp = z - Vstar
    rootss = solve(exp, y)
    desired_root=0
    for lamb in rootss:
        for i in np.arange(0,N):
            if (-1)*e[i]/((1-e[i])*(1+lamb*g0*Isp[i])) < 0 :
                continue
            else:
                desired_root = lamb
                break


    Pistar = 1
    for i in np.arange(0,N):
        Pistage[i]=(-1)*e[i]/((1-e[i])*(1+desired_root*g0*Isp[i]))
        Vstage[i] = (-1)*g0*Isp[i]*math.log(e[i] + (1-e[i])*Pistage[i])
        Pistar = Pistar*Pistage[i]

    print("\nMaximum Payload Fraction:",Pistar)



mp = [None]*N
ms = [None]*N
mo = [None]*N

for i in reversed(range(0,N)) :
    if i == (N -1) :
        ms[i] = e[i]*((1 - Pistage[i])/Pistage[i])*mL
        mp[i] = (1 - e[i])*((1 - Pistage[i])/Pistage[i])*mL
        mo[i] = mp[i] + ms[i] + mL
    else :
        ms[i] = e[i]*((1 - Pistage[i])/Pistage[i])*mo[i+1]
        mp[i] = (1 - e[i])*((1 - Pistage[i])/Pistage[i])*mo[i+1]
        mo[i] = mp[i] + ms[i] + mo[i+1]

print("\nThus, Total Lift Off Mass : %.5f Kgs" %mo[0])
