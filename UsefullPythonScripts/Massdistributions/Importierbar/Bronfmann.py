import numpy as np

MP = [[2.25,0.11,0.42,39,77],            #MESSWERTE R, North, South, z0, z1/2 (ZH)
[2.75,0.07,0.49,36,80],
[3.25,0.27,0.38,0,61],
[3.75,0.34,0.63,-8,65],
[4.25,0.91,0.65,1,71],
[4.75,1.23,0.49,-10,72],
[5.25,1.29,0.72,-1,82],
[5.75,1.15,0.63,-4,83],
[6.25,0.95,0.71,-19,73],
[6.75,0.56,0.87,-22,63],
[7.25,0.73,0.93,-14,58],
[7.75,0.64,0.83,-9,72],
[8.25,0.41,0.5,-4,80],
[8.75,0.28,0.37,13,66],
[9.25,0.63,0.29,-4,23],
[9.75,0.22,0.18,-20,147]]

x=np.linspace(1,16,16)
y=x.copy()
north = x.copy()
south = x.copy()
z0=x.copy()
zH=x.copy()
for i in range(len(y)):
    x[i]=MP[i][0]
    north[i]=MP[i][1]
    south[i]=MP[i][2]
    y[i]=(north[i]+south[i])/2
    z0[i]=MP[i][3]
    zH[i]=MP[i][4]

"""Linear Interpoliert und Rescaliert auf R_Sonne = 8.5 kpc"""
def Bronf_R(R):
    R=R/0.85
    if R<=2.25:
        return 0
    else:
        if R>=9.75:
            return 0
        else:
            i = int(R*2 -4.5)+1
            return 1/0.85*(y[i-1]+(y[i]-y[i-1])/(x[i]-x[i-1])*(R-x[i-1]))

def Bronf_z0(R):
    R=R/0.85
    if R<=2.25:
        return 0
    else:
        if R>=9.75:
            return 0
        else:
            i = int(R*2 -4.5)+1
            return 1/0.85*(z0[i-1]+(z0[i]-z0[i-1])/(x[i]-x[i-1])*(R-x[i-1]))

def Bronf_zH(R):
    R=R/0.85
    if R<=2.25:
        return 1000
    else:
        if R>=9.75:
            return 1000
        else:
            i = int(R*2 -4.5)+1
            return 1/0.85*(zH[i-1]+(zH[i]-zH[i-1])/(x[i]-x[i-1])*(R-x[i-1]))

"""Zusammenführen der Komponenten"""
def n_H2	(R,z):
    return Bronf_R(R)*np.exp(-(z-Bronf_z0(R))**2*np.log(2)/(Bronf_zH(R))**2) 
