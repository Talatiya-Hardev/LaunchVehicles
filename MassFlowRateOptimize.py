# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 21:01:07 2022

@author: Hardev
"""

"""
This Code is Applicable for Loss-Optimization of only Large-Scale Launch Vehicles
"""


#Programm to Plot Gravity Loss vs Burn Rate and Give Optimum Burnout Parameters

import matplotlib.pyplot as plt
import numpy as np
import math

from ISA import ISA

N = 1000
g0=9.81

def Get_GL_DL(m0,mp,Isp,BR):
    mb = m0 - mp
    Y = mp/m0
    N = 1000
    Vb = g0*Isp*np.log(m0/mb)         #IBV, No Gravity or Drag Effects

    FT = mp/BR
    Vbg = Vb - g0*FT
    t = np.linspace(0,FT,N)
    Vbt = g0*Isp*np.log(m0/(m0 - BR*t)) - g0*t
    Hbg = (m0*g0*Isp/BR)*((1 - Y)*np.log(1 - Y) + Y) - 0.5*g0*FT*FT
    Hbt = (m0*g0*Isp/BR)*((1 - (BR*t/m0))*np.log(1 - (BR*t/m0)) + (BR*t/m0)) - 0.5*g0*t*t
    rhos = np.zeros_like(t)
    #print(Vb,FT,Vbg,Hbg)
    DP = []

    for i in range(len(Hbt)) :
        Temp,rhos[i],P = ISA(Hbt[i],0)
        DP.append(0.5*rhos[i]*Vbt[i]*Vbt[i])


    THDP = np.zeros(shape=(3,N))
    THDP[0] = t
    THDP[1] = Hbt
    THDP[2] = DP

    THDP = THDP.transpose()

    abc = np.zeros(shape=(N,3))
    abc[:,0] = t
    abc[:,1] = Hbt
    abc[:,2] = DP

    #print(THDP[0])

    sorted_array = THDP[np.argsort(THDP[:,2])]
    t_maxDP = sorted_array[-1,0]
    h_maxDP = sorted_array[-1,1]
    DP_maxDP = sorted_array[-1,2]

    #print(t_maxDP,h_maxDP)

    T2,rho_maxDP,P = ISA(h_maxDP,0)
    #T2,rhoat13200m,P = ISA(13200,0)
    V_maxDP = g0*Isp*np.log(m0/(m0 - BR*t_maxDP)) - g0*t_maxDP
    Drag_maxDP = 0.5*rho_maxDP*V_maxDP*V_maxDP*Sref*Cd0
    Mass_maxDP = m0 - BR*t_maxDP
    #print(Mass_maxDP)

    avg_Drag_accn = Drag_maxDP/(2*Mass_maxDP)
    #print(avg_Drag_accn)

    """IMP : If Alt vs DP does not form like the one in the TB,
    then you can't apply Avg_AD = MaxDrag/2M
    """

    Vbgd = Vbg - (avg_Drag_accn*FT)
    Hbgd = Hbg - (0.5*avg_Drag_accn*FT*FT)
    # print(Hbgd)


    gldl = 0.5*Vb*Vb - 0.5*Vbgd*Vbgd - g0*Hbgd

    return gldl,avg_Drag_accn

#############################################################################################################

while True:
    D = input("Use default Launch Vehicle Parameters? : 1 for Yes. 0 for Manual : ")
    if D =="1":
        mp,ms,mL,Isp,g0 = 60000,18000,2000,240,9.81
        Sref,Cd0 = 3.14159265,1
    else :
         mL = float(input("\nEnter Payload Mass :"))
         mp = float(input("Enter Propellant Mass :"))
         ms = float(input("Enter Structural Mass :"))
         Isp = float(input("Specify Average Specific Impulse :"))
         Sref = float(input("Enter Reference Surface Area :"))
         Cd0 = float(input("Enter Drag Coefficient :"))

    m0 = mL+mp+ms
    mb = m0 - mp
    Y = mp/m0
    print("\nPayload Mass : ",mL,"\nPropellant Mass : ",mp,"\nStructural Mass : ",ms,"\nSpecific Impulse : ",Isp,"\nReference Surface Area : ",Sref,"\nDrag Coefficient : ",Cd0)

    BR = np.linspace(300,4000,1000,endpoint=True,dtype = float)
    gldl = np.zeros_like(BR)
    avg_aD = np.zeros_like(BR)

    Vb = g0*Isp*np.log(m0/mb)         #IBV, No Gravity or Drag Effects

    gldl00,avg00 = Get_GL_DL(m0, mp, Isp, 600)

    for i in range(len(BR)) :
        gldl[i],avg_aD[i] = Get_GL_DL(m0,mp,Isp,BR[i])
        if gldl[i] < 0 :                    #NEW LINES ADDED HERE 3
            gldl[i] = float("nan")
            avg_aD[i] = float("nan")


    Loss_BR = np.zeros(shape=(N,2))
    Loss_BR[:,0] = BR
    Loss_BR[:,1] = gldl

    array_sorted = Loss_BR[np.argsort(Loss_BR[:,1])]
    Optimized_BR = array_sorted[0,0]
    Opt_T = mp/Optimized_BR

    Opt_gldl,Opt_aD = Get_GL_DL(m0, mp, Isp, Optimized_BR)

    Vbgopt = Vb - g0*Opt_T
    Hbgopt = (m0*g0*Isp/Optimized_BR)*((1 - Y)*np.log(1 - Y) + Y) - 0.5*g0*Opt_T*Opt_T
    Vbgd = Vb - g0*Opt_T - (Opt_aD*Opt_T)
    Hbgd = (m0*g0*Isp/Optimized_BR)*((1 - Y)*np.log(1 - Y) + Y) - 0.5*g0*Opt_T*Opt_T - (0.5*Opt_aD*Opt_T*Opt_T)

    #print("\n",Optimized_BR, Opt_aD)
    print("\n Optimized Mass Flow Rate: ",Optimized_BR,"Kg/s")
    #print(Vbgopt,Vbgd,Hbgopt,Hbgd)
    print("\nUsing Optimized Mass Flow Rate")
    print("Burnout Velocity: ",Vbgd,"m/s")
    print("Burnout Height: ",Hbgd,"m")



    plt.plot(BR,gldl,color = "red")
    plt.xlabel("Mass Flow Rate (Kg/s)")
    plt.ylabel("Combined Losses (Gravity and Drag) (Joules)")
    plt.title("Combined Losses vs Burn Rate")
    plt.grid()
    plt.show()




    f = input("\nPress any key to continue or 0 to exit the program. ")
    if f == "":
        continue
    elif f == '0' :
        break
