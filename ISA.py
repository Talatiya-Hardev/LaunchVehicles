#Author : Hardev Talatiya

import matplotlib.pyplot as plt
import numpy as np
import math

def ISA(h,t) :

    if t ==0 :
        h0,g0,T0,rho0,P0,R = 0,9.81,288.15,1.225,101325,287
    else :
        T0 = 288.15 + t
        h0,g0,rho0,P0,R = 0,9.81,1.225,101325,287

    if h <= 11000 :  #Gradient Layer
        L = -0.0065
        T = T0 + L*(h - h0)
        k = g0/(L*R)
        rho = rho0*pow((T/T0),(-(k+1)))
        P = P0*pow((T/T0),(-k))
        return T,rho,P

    elif 11000 < h <= 20100 :     #IsoThermal Layer
        T11,rho11,P11 = ISA(11000,t)
        T = T11
        rho = rho11*math.exp((-g0/(R*T))*(h-11000))
        P = P11*math.exp((-g0/(R*T))*(h-11000))
        return T,rho,P

    elif 20100 < h <= 32000 :   #Gradient Layer
        T21,rho21,P21 = ISA(20100,t)
        L = 0.001
        k = g0/(L*R)
        T = T21 + L*(h - 20100)
        rho = rho21*pow((T/T21),(-(k+1)))
        P = P21*pow((T/T21),(-k))
        return T,rho,P

    elif 32000 < h <= 47000 :  #Gradient Layer
        T32,rho32,P32 = ISA(32000,t)
        L = 0.0028
        k = g0/(L*R)
        T = T32 + L*(h - 32000)
        rho = rho32*pow((T/T32),(-(k+1)))
        P = P32*pow((T/T32),(-k))
        return T,rho,P

    elif 47000 < h <= 51000 :     #IsoThermal Layer
        T47,rho47,P47 = ISA(47000,t)
        T = T47
        rho = rho47*math.exp((-g0/(R*T))*(h-47000))
        P = P47*math.exp((-g0/(R*T))*(h-47000))
        return T,rho,P

    elif 51000 < h <= 71000 :  #Gradient Layer
        T51,rho51,P51 = ISA(51000,t)
        L = -0.0028
        k = g0/(L*R)
        T = T51+ L*(h - 51000)
        rho = rho51*pow((T/T51),(-(k+1)))
        P = P51*pow((T/T51),(-k))
        return T,rho,P

    elif 71000 < h <= 86000 :  #Gradient Layer
        T71,rho71,P71 = ISA(71000,t)
        L = -0.002
        k = g0/(L*R)
        T = T71 + L*(h - 71000)
        rho = rho71*pow((T/T71),(-(k+1)))
        P = P71*pow((T/T71),(-k))
        return T,rho,P

# Right Now This file is being imported somehwere hence below lines are commented. To use this by itself remove the comments from below
# t = float(input("Default ISA (Press 0) or ISA +/- Temperature? : "))
# h = float(input("Enter height in metres : "))
# T,rho,P = ISA(h,t)
# print(T,rho,P)
