# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 02:10:55 2022

@author: Hardev
"""
import matplotlib.pyplot as plt
import numpy as np
import math
from ISA import ISA

g0 = 9.81
Sref=3.14159265
Cd0=1

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

    T2,rho_maxDP,P = ISA(h_maxDP,0)
    #T2,rhoat13200m,P = ISA(13200,0)
    V_maxDP = g0*Isp*np.log(m0/(m0 - BR*t_maxDP)) - g0*t_maxDP
    Drag_maxDP = 0.5*rho_maxDP*V_maxDP*V_maxDP*Sref*Cd0
    Mass_maxDP = m0 - BR*t_maxDP

    avg_Drag_accn = Drag_maxDP/(2*Mass_maxDP)
    #print(avg_Drag_accn)

    Vbgd = Vbg - (avg_Drag_accn*FT)
    Hbgd = Hbg - (0.5*avg_Drag_accn*FT*FT)
    #print(Hbgd)


    gldl = 0.5*Vb*Vb - 0.5*Vbgd*Vbgd - g0*Hbgd
    
    return gldl,avg_Drag_accn,Vbgd,Hbgd