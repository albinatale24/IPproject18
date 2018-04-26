# Compressed Air Cannon
# This code has been developed by Alberto Natale
# University of Southampton 2018

import numpy as n
import matplotlib.pyplot as plt
import xlsxwriter

P0 = 200000                   #in chamber pressure
Pb0 = 100000                  #atmospheric pressure
rmax = 0.8                    #valve coefficient
B = 3.11 * (10**19)           #engineering energy constant
Cv = 1.93                     #valve coefficient flow coefficient
Gg = 1                        #gravitational acceleration
g = 9.81                      #gravity  
T = 293                       #temperature
Z = 1                         #constant
kb = 1.38064852 * (10**-23)   #constant
d = 0.05                      #ball resting distance from valve opening
rad = 0.01                    #radius of ball
V0 = 0.004                    #chamber volume
L = 0.32                      #length of barrel
Patm = 100000                 #atmospheric pressure
m = 0.004                     #mass of the ball
y = 1.4                       #air coefficient
cd = 0                        #drag coefficient  
ad = 1.225                    #air density  
angle = 0                     #angle of fire 
yo = 0                        #heigh of fire  

a = n.radians(angle)            
f = 10*m*g*n.cos(a)*0
A = 3.1416*((rad)**2.0)

def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def vadd(P0):
    """Function given an initial chamber pressure return cannon exit velocity assuming adiabatic expansion"""

    vadd = (((2/m)*(((((P0*V0)/(y-1))*(1-((V0/(A*L+V0))**(y-1)))))-A*L*Patm-L*f-L*m*g*n.sin(a))))**0.5
    return vadd

def viss(P0):
    """Function given an initial chamber pressure return cannon exit velocity assuming isothermal expansion"""
    
    viss = ((2.0/m)*(((P0*V0)*n.log(1.0+((A*L)/V0)))-A*L*Patm-L*f-L*m*g*n.sin(a)))**0.5
    return viss

def graph_comp():
    """Function creates a graph comparing Adiabatic, Isothermal and Valve model"""
    
    p = n.linspace(100000, 600000, 100)
    vad = vadd(P0)
    vis = viss(P0)
    vlv = euler_loop(P0)
    
    plt.figure(figsize=(9,6))
    plt.plot(p,vad, label='Adiabatic')
    plt.plot(p,vis, label='Isothermal')
    plt.plot(p,vlv, label='Valve Model')
    plt.legend(fontsize=14)
    plt.xlabel('Pressure Pa', fontsize=15)
    plt.ylabel('Velocity m/s', fontsize=15)
    plt.title('Exit Velocity over Pressures',fontsize=18, fontweight='bold')
    plt.savefig('Comparatives Velocities Adibatic vs Isothermal vs Valve.png')


def euler_loop(P0):
    """Function given an initial chamber pressure return cannon exit velocity considering valve effect"""
    
    r0 = (P0 - Pb0) / P0
    Q0 = B * P0 * Cv * (1 - ( r0 / (3*rmax))) * ((r0 / (Gg*T*Z))**0.5)
    """metter q a 0"""
    N0 = (P0*V0) / (kb * T)
    Nb0 = (Pb0*A*d) / (kb * T)
    
    """1"""
    dt = 0.000001             # this value sets up the incremental dt 
  
    Ndt = N0 - Q0 * dt
    Nbdt = Nb0 + Q0 * dt
    Pdt = (Ndt*kb*T) / (V0)
    Pbdt = (Nbdt*kb*T) / (A*d)
    rdt = (Pdt - Pbdt) / Pdt
    Qdt = B * Pdt * Cv * (1 - (rdt / (3 * rmax))) * ((rdt / (Gg*T*Z))**0.5)
    xdt = 0
    vdt = 0
    adt = ((((Pbdt-Patm)*A))-(0.5*ad*(vdt**2)*A*cd)-(m*g*n.sin(a))-f)/m
    
    """2"""
    N2dt = Ndt - Qdt * dt
    Nb2dt = Nbdt + Qdt * dt
    P2dt = (N2dt * kb * T) / V0
    x2dt = xdt + vdt * dt + (1/2) * adt * (dt**2)
    v2dt = vdt + adt*dt
    Pb2dt = (Nb2dt * kb * T)/(A * (d + x2dt))
    a2dt = (((Pb2dt-Patm)*A) - (0.5*ad*(vdt**2)*A*cd) - (m*g*n.sin(a))-f) / m
    r2dt = (P2dt - Pb2dt) / P2dt
    Q2dt = B * P2dt * Cv * (1 - (r2dt / (3 * rmax))) * ((r2dt / (Gg*T*Z))**0.5)

    data = [N2dt, Nb2dt, P2dt, x2dt, Pb2dt, r2dt, Q2dt, v2dt, a2dt]  
       #     0       1     2     3     4      5     6    7      8 
            
    while data[3] <= L:
        
        N3dt = data[0] - data[6] * dt
        Nb3dt = data[1] + data[6] * dt
        P3dt = (data[0] * kb * T) / V0
        x3dt = data[3] + data[7] * dt + (1/2) * data[8] * (dt**2)
        v3dt = data[7] + data[8]*dt
        a3dt = (((data[4]-Patm)*A)-(0.5*ad*(data[7]**2)*A*cd)-(m*g*n.sin(a))-f)/m
        Pb3dt = (data[1] * kb * T)/(A * (d + data[3]))
        r3dt = (data[2] - data[4]) / data[2]

        acls.append(a3dt)
        diss.append(x3dt)

        if r3dt <= rmax:
            Q3dt = B * data[2] * Cv * (1 - (data[5] / (3 * rmax))) * ((data[5] / (Gg*T*Z))**0.5)
        else:
            Q3dt = (2/3) * B * data[2] * Cv * ((rmax/(Gg*T*Z))**0.5) 
    
        del data
        data = [N3dt, Nb3dt, P3dt, x3dt, Pb3dt, r3dt, Q3dt, v3dt, a3dt]
        
    return v3dt

def distance(v):
    """Given velocity of firing, function returns distance flew by object """
    
    return ((v*n.cos(a))/g)*(v*n.sin(a)+(((v*n.sin(a))**2)+2*g*yo)**0.5)  
