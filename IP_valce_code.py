import numpy as n
import matplotlib.pyplot as plt
import xlsxwriter

P0 = 200000                  #in chamber pressure
Pb0 = 100000                  #atmospheric pressure

rmax = 0.8                    #valve coefficient
B = 3.11 * (10**19)           #engineering energy constant
Cv = 1.93                     #valve coefficient flow coefficient
Gg = 1                        #gravitational acceleration

T = 293
Z = 1
kb = 1.38064852 * (10**-23)
d = 0.05
rad = 0.01                    #radius of ball
A = 3.1416*((rad)**2.0)
V0 = 0.004
L = 0.32                       #length of barrel
Patm = 100000
m = 0.004
y = 1.4
cd = 0.4

angle = 0
a = n.radians(angle)
g = 9.81
yo = 1
f = 10*m*g*n.cos(a)*0

ad = 1.225

acls = []
diss = []

def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def vadd(P0, A):
    vadd = (((2/m)*(((((P0*V0)/(y-1))*(1-((V0/(A*L+V0))**(y-1)))))-A*L*Patm-L*f-L*m*g*n.sin(a))))**0.5
    return vadd

def vadd2(m, A, Po, Vo, y, L, Pa, f, a, g, yo):
    vadd = (((2/m)*(((((Po*Vo)/(y-1))*(1-((Vo/(A*L+Vo))**(y-1)))))-A*L*Pa-L*f-L*m*g*n.sin(a))))**0.5
    return vadd

def viss(m, A, Po, Vo, y, L, Pa, f, a, g, yo):
    viss = ((2.0/m)*(((Po*Vo)*n.log(1.0+((A*L)/Vo)))-A*L*Pa-L*f-L*m*g*n.sin(a)))**0.5
    return viss

def graph_comp():
    p = n.linspace(100000, 600000, 100)
    show_vv = []
    vad = vadd2(m, A, p, V0, y, L, Patm, f, a, g, yo)
    vis = viss(m, A, p, V0, y, L, Patm, f, a, g, yo)
    
    for i in p:
        show_vv.append(loop(i, A))
    
    plt.figure(figsize=(9,6))
    plt.plot(p,vad, label='Adiabatic')
    plt.plot(p,vis, label='Isothermal')
    plt.plot(p, show_vv, label='Valve Model')
    plt.legend(fontsize=14)
    plt.xlabel('Pressure Pa', fontsize=15)
    plt.ylabel('Velocity m/s', fontsize=15)
    plt.title('Exit Velocity over Pressures',fontsize=18, fontweight='bold')
    plt.savefig('Comparatives Velocities Adibatic vs Isothermal vs Valve.png')


def loop(P0, A):
    r0 = (P0 - Pb0) / P0
    Q0 = B * P0 * Cv * (1 - ( r0 / (3*rmax))) * ((r0 / (Gg*T*Z))**0.5)
    """metter q a 0"""
    N0 = (P0*V0) / (kb * T)
    Nb0 = (Pb0*A*d) / (kb * T)
    
    """1"""
    
    dt = 0.000001
    
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
    return ((v*n.cos(a))/g)*(v*n.sin(a)+(((v*n.sin(a))**2)+2*g*yo)**0.5)  
    
    
vv = loop(P0, A)
dis = ((vv*n.cos(a))/g)*(vv*n.sin(a)+(((vv*n.sin(a))**2)+2*g*yo)**0.5)  

vvad = vadd(P0, A)
disad = ((vvad*n.cos(a))/g)*(vvad*n.sin(a)+(((vvad*n.sin(a))**2)+2*g*yo)**0.5) 

"""
ek = 0.25                  #error coefficient
adjusted_speed = loop(P0, A)*(ek)
dissadj = ((adjusted_speed*n.cos(a))/g)*(adjusted_speed*n.sin(a)+(((adjusted_speed*n.sin(a))**2)+2*g*yo)**0.5) 
"""


print('\n', '\n')   
print(' VALVE', '\n', "Velocity m/s  ", loop(P0, A),'\n',"Distance m    ",  dis,)
#print('\n', "Adj. Velocity m/s  ", adjusted_speed, '\n', "Adj. Distance m    ",  dissadj,)
print('\n', '\n')
print(' ADIABATIC', '\n', "Velocity m/s  ", vadd(P0, A),'\n',"Distance m    ",  disad, '\n')
    
"""

#Pressure Velocity Plot

ps = range(150000, 800000, 5000)
v_val = []
v_ad = []

for i in ps:
    v_val.append(loop(i, A))
    
for i in ps:
    v_ad.append(vadd(i, A))

plt.figure(figsize=(9,6))
plt.plot(ps, v_val, label='Valve Calculated')
plt.plot(ps, v_ad, label='Adiabatic Expansion')
plt.xlabel('Pressure Pa',fontsize=14)
plt.ylabel('Exit Velocity m/s',fontsize=14)
plt.title('Velocity over Pressure', fontsize=16, fontweight='bold')
plt.legend(fontsize=12)    
plt.savefig('Velocity_over_Pressure.png')

   
#Radius Space

rs = list(frange(0.02, 0.05, 0.0005))
As = []
vrad = []
vrval = []

for i in rs:
    As.append(3.1416*(i)**2.0)

for i in As:
    vrad.append(vadd(P0, i))

for i in As:
    vrval.append(loop(P0, i))


plt.figure(figsize=(9,6))
plt.plot(rs,vrval, label='Valve Calculated')
plt.plot(rs,vrad, label='Adiabatic Expansion')
plt.legend(fontsize=12,loc=1)
plt.xlabel('m', fontsize=14)
plt.ylabel('m/s', fontsize=14)
plt.title('Exit Velocity over Ball Radius',fontsize=16, fontweight='bold')
plt.savefig('Radius_Effect_over_Velocity.png')

"""


#Acceleration Space



plt.figure(figsize=(9,6))
plt.plot(diss, acls,'g.', markersize=1)
plt.legend()
plt.xlabel('Distance in Barrel m', fontsize=15)
plt.ylabel('Acceleration m/s2', fontsize=15)
plt.title('Change of Acceleration in the Barrel', fontsize=18)
plt.savefig('Change of Acceleration in the Barrel.png')



workbook = xlsxwriter.Workbook('CannonNew1sp.xlsx')
worksheet = workbook.add_worksheet()


"""pressure_range = range(130000, 400000, 5000)"""

#pressure_range = [150000, 200000, 250000, 300000, 350000]
pressure_range = range(130000, 300000, 2500)


row = 0
for i in pressure_range:
        vsv = loop(i, A)
        worksheet.write(row, 0, vsv)
        worksheet.write(row, 1, distance(vsv))
        worksheet.write(row, 2, i)
        row = row + 1

theoretical = []
practical = [33.6, 45.6, 55.2, 62.4, 69.6]

for i in pressure_range:
    theoretical.append(loop(i, A))
    
plt.figure(figsize=(9,6))
plt.plot(pressure_range,theoretical,'o--',label='Theoretical')
plt.plot(pressure_range,practical,'ro--', label='Experimental')
"""plt.errorbar(pressure_range, practical, yerr=2.4, xerr=10000)"""
plt.legend(fontsize=14,loc=0)
plt.xlabel('Pressure Pa', fontsize=15)
plt.ylabel('Velocity m/s', fontsize=15)
plt.title('Exit Velocity over Pressures',fontsize=18, fontweight='bold')
plt.savefig('Comparatives Velocities B.png')

print(theoretical)