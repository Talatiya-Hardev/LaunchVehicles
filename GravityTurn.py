# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 12:49:37 2022

@author: Hardev T
"""

import matplotlib.pyplot as plt
import numpy as np
import math

g0 = 9.81
"""
Masses in tonnes
M0 = Total Mass, Mp = Propellant Mass, V0 = Initial Velocity during Pitch Kick, h0 and x0 = initial altitude and distance covered
theta0 = Pitch Kick/Initial pitch angle
thetab = Burnout pitch angle

"""
#Enter the Values Below
#For Longer Trajectories you will have to zoom it closely to get look of all 3 trajectories (Especially Constant Velocity)

M0,Mp,Isp,V0,h0,x0,theta0 = 80,60,240,90,400,400,2   #In Tonnes and degrees
theta_b = 20         #Enter Pitch angle at Burnout



def ConstPitch(M0,Mp,Isp,V0,h0,theta0,theta_b) :

    N=10000
    T0=0
    theta0 = math.radians(theta0)
    theta_b = math.radians(theta_b)
    theta_t = np.linspace(theta0,theta_b,N)
    q0 = g0*math.sin(theta0)/V0
    t = T0 + (theta_t - theta0)/(q0)
    #theta_t = theta0 + q0*t

    Vt = g0*np.sin(theta_t)/q0

    location = np.zeros([N,2])
    location[:,0] = (g0/(2*q0*q0)) * ( (theta_t - theta0) -  (( np.sin(2*theta_t) - np.sin(2*theta0) )/2) ) + x0
    location[:,1] = (g0/(4*q0*q0)) * (np.cos(2*theta0) - np.cos(2*theta_t)) + h0
    #print(location[:,0])
    plt.plot(location[:,0], location[:,1],label="Constant Pitch")
    #plt.plot(t,location[:,0])
    plt.title("Altitude vs Horizontal Distance ")
   # plt.grid()
    print("Constant Pitch Trajectory: ",t[-1], "s")


def ConstVelocity(M0,Mp,Isp,V0,h0,theta0,T0,theta_b) :
    N=10000
    T0 = 0
    theta0 = math.radians(theta0)
    theta_b = math.radians(theta_b)
    theta_t = np.linspace(theta0,theta_b,N)
    t = T0 + (V0/g0)*np.log(np.tan(theta_t/2)/np.tan(theta0/2))
    #print(t)
    q0 = g0*math.sin(theta0)/V0
    T0 = 0
    dt = t - T0
    #theta_t = 2*np.arctan( math.tan(theta0/2)*np.exp(g0*dt/V0) )
    dtheta_t = theta_t - theta0
    Vt = g0*np.sin(theta_t)/q0
    location = np.zeros([N,2])
    location[:,0] = (V0**2/g0)*dtheta_t + x0
    location[:,1] = (V0**2/g0)*np.log( np.sin(theta_t)/np.sin(theta0) ) + h0
    plt.plot(location[:,0], location[:,1],label="Constant Velocity")
    #plt.plot(t,location[:,0])
    plt.title("Altitude vs Horizontal Distance ")
    #plt.grid()
    print("Constant Velocity Trajectory: ",t[-1], "s")




def ConstTbyM(M0,Mp,Isp,V0,h0,theta0,T0,theta_b,n0) :    # Set your T/M value in the variable n0
           #T/m value
    N=1000
    T0 = 0
    theta0 = math.radians(theta0)
    theta_b = math.radians(theta_b)
    theta_t = np.linspace(theta0,theta_b,N)
    Dr = np.power(np.tan(theta0/2),n0-1) + np.power(np.tan(theta0/2),n0+1)
    k = V0/Dr
    Tterm1 = (np.power(np.tan(theta_t/2),n0-1)/(n0-1)) + (np.power(np.tan(theta_t/2),n0+1)/(n0+1))
    Tterm2 = (np.power(np.tan(theta0/2),n0-1)/(n0-1)) + (np.power(np.tan(theta0/2),n0+1)/(n0+1))
    t = T0 + (k/g0)*(Tterm1 - Tterm2)


    q0 = g0*math.sin(theta0)/V0
    T0 = 0
    #t = np.linspace(T0,T,N)
    dt = t - T0
    #theta_t = 2*np.arctan( math.tan(theta0/2)*np.exp(g0*dt/V0) )
    dtheta_t = theta_t - theta0
    Vt = g0*np.sin(theta_t)/q0
    location = np.zeros([N,2])

    Hterm1 = np.power(np.tan(theta_t/2),2*(n0-1))/(n0-1) - np.power(np.tan(theta_t/2),2*(n0+2))/(n0+2)
    Hterm2 = np.power(np.tan(theta0/2),2*(n0-1))/(n0-1) - np.power(np.tan(theta0/2),2*(n0+2))/(n0+2)
    Xterm1 = np.power(np.tan(theta_t/2),2*n0-1)/(2*n0-1) + np.power(np.tan(theta_t/2),2*n0+1)/(2*n0+1)
    Xterm2 = np.power(np.tan(theta0/2),2*n0-1)/(2*n0-1) + np.power(np.tan(theta0/2),2*n0+1)/(2*n0+1)

    location[:,0] = (2*(k**2)/g0)*(Xterm1 - Xterm2) + x0
    location[:,1] = k**2/(2*g0)*(Hterm1 - Hterm2) + h0
    plt.plot(location[:,0], location[:,1],label="Constant T/M")
    #plt.plot(t,location[:,0])
    plt.title("Altitude vs Horizontal Distance ")
    #plt.grid()
    print("Constant T/M Trajectory: ",t[-1], "s")
T0=0

print("To pitch from",theta0, "degrees to",theta_b,"degrees,the time taken is:")
ConstPitch(M0,Mp,Isp,V0,h0,theta0,theta_b)
ConstVelocity(M0,Mp,Isp,V0,h0,theta0,T0,theta_b)
ConstTbyM(M0,Mp,Isp,V0,h0,theta0,T0,theta_b,1.5)
plt.legend()
plt.grid()
plt.xlabel("Horizontal Distance (m)")
plt.ylabel("Altitude (m)")
plt.show()
