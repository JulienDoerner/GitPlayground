import numpy as np
import matplotlib.pyplot as plt
from pylab import *


"""Außen"""
RS=8500.
def a(R):
    if R<=RS:
        return 1
    else:
        return R/RS

def n_A_H2(R,z):	# Dichte molekularer Wasserstoff
	return np.exp(-((R-4500)**2-(RS-4500)**2)/(2900**2))*(R/RS)**-0.58*np.exp(-(z/81*a(R)**0.58)**2)

def n_A_HI_cold(R,z):	# Kalter Atomarer Wasserstoff
	T1 = 0.859*np.exp(-(z/(127*a(R)))**2)
	T2 = 0.047*np.exp(-(z/(318*a(R)))**2)
	T3 = 0.094*np.exp(-(np.abs(z)/(403*a(R))))
	return 0.304/(a(R)**2)*(T1+T2+T3)
	
def n_A_HI_warm(R,z):	# Warmer Atomarer Wasserstoff
	T1 = (1.745 - 0.289/a(R))*np.exp(-(z/(127*a(R)))**2)
	T2 = (0.473 - 0.070/a(R))*np.exp(-(z/(318*a(R)))**2)
	T3 = (0.283 - 0.142/a(R))*np.exp(-abs(z)/(403*a(R)))
	return 0.226/a(R)*(T1+T2+T3)

def n_A_HI(R,z):		# Atomarer Wasserstoff insg.
	return n_A_HI_cold(R,z) + n_A_HI_warm(R,z)

def n_A_HII_warm(R,z): 	# Warmer Ionisierter Wasserstoff
	T1 = 0.0237*np.exp(-(R**2-RS**2)/37000**2)*np.exp(-abs(z)/1000)
	T2 = 0.0013*np.exp(-((R-4000)**2 - (RS-4000)**2)/2000**2)*np.exp(-abs(z)/150)
	return T1 + T2
	
def n_A_HII_hot(R,z):	# Heisser Ionisierter Wasserstoff
	T1 = 0.12*np.exp(-(R-RS)/4900) 
	T2 = 0.88*np.exp(-((R-4500)**2-(RS-4500)**2)/2900**2)*a(R)**-1.65*np.exp(-abs(z)/(1500*a(R)**1.65))
	return 4.8*10**-4 * (T1+T2)

def n_A_HII(R,z): 		# Ionisierter Wasserstoff insg.
	return n_A_HII_warm(R,z)+ n_A_HII_hot(R,z)

""" Innen"""

# Konstanten fuer CMZ
xc =-50			# Position Mitte in allg Koordinaten
yc = 50
TettaC = 70

#Konstanten fuer DISK
alpha = 13.5
beta = 20.
TettaD = 48.5

# Abmessungen in CMZ Koordinaten
XMAX=250		
XC = XMAX/2
LC = XMAX/(2*np.log(2)**0.25)
HC = 18.
HC2 = 54.

# Abmessungen in DISK Koordinaten
XD = 1200
LD = 438.
HD = 42.
HD2 = 120.

#Konstanten fuer HII -WIM-
y3 = -10
z3= -20
L3 = 145.
H3 = 26.
L2 = 3700.
H2 = 140.
L1 = 17000
H1=950.

#Konstanen fuer HII VHIM
alphaVH = 21
LVH=162
HVH = 90


def Bogenmass(x):			# Trafo ins Bogenmass fuer Winkel zur Berechnung
	return x*np.pi/180

def cos(x):					# Cos FKT fuer Gradmass
	x=Bogenmass(x)
	return np.cos(x)
def sin(x): 				# Sin FKT fuer Gradmass
	x=Bogenmass(x)
	return np.sin(x)
def sech2(x):
	return np.cosh(x)**-2

def CMZ_X_Trafo(x,y):
	return (x-xc)*cos(TettaC) +(y-yc)*sin(TettaC)
def CMZ_Y_Trafo(x,y):
	return -(x-xc)*sin(TettaC) +(y-yc)*cos(TettaC)

def DISK_X_Trafo(x,y,z):
	return x*cos(beta)*cos(TettaD) - y*(sin(alpha)*sin(beta)*cos(TettaD) -cos(alpha)*sin(TettaD))-z*(cos(alpha)*sin(beta)*cos(TettaD) +sin(alpha)*sin(TettaD))
def DISK_Y_Trafo(x,y,z):
	xT= x*cos(beta)*sin(TettaD)
	yT = y*(sin(alpha)*sin(beta)*sin(TettaD) +cos(alpha)*cos(TettaD))
	zT = z*(cos(alpha)*sin(beta)*sin(TettaD) -sin(alpha)*sin(TettaD))
	return -xT+yT+zT
def DISK_Z_Trafo(x,y,z):
	xT = x*sin(beta)
	yT = y*sin(alpha)*cos(beta)
	zT = z*cos(alpha)*cos(beta)
	return xT+yT+zT
	
#Mollekularer Wasserstoff im CMZ,
def n_I_H2_CMZ(x0,y0,z0): 			# Eingabe in Urspruenglichen koordinaten
	x = CMZ_X_Trafo(x0,y0)
	y = CMZ_Y_Trafo(x0,y0)
	XY_Help = ((np.sqrt(x**2+(2.5*y)**2)-XC)/LC)**4
	return 150*np.exp(-XY_Help)*np.exp(-(z0/HC)**2)

#Atomarer Wasserstoff im CMZ 
def n_I_HI_CMZ(x0,y0,z0):			#Eingabe in Urspruenglichen Koordinaten
	x=CMZ_X_Trafo(x0,y0)
	y=CMZ_Y_Trafo(x0,y0)
	A=np.sqrt(x**2 +(2.5*y)**2)
	B= (A-XC)/LC
	XY_Help=B**4
	Z = (z0/HC2)**2
	return 8.8*np.exp(-XY_Help)*np.exp(-Z)
	
#Mollekularer Wasserstoff in der DISK
def n_I_H2_DISK(x0,y0,z0):
	x= DISK_X_Trafo(x0,y0,z0)
	y= DISK_Y_Trafo(x0,y0,z0)
	z=DISK_Z_Trafo(x0,y0,z0)
	return 4.8*np.exp(-((np.sqrt(x**2 + (3.1*y)**2) - XD)/LD)**4)*np.exp(-(z/HD)**2)

#Atomarer Wasserstoff in der DISK
def n_I_HI_DISK(x0,y0,z0):
	x= DISK_X_Trafo(x0,y0,z0)
	y= DISK_Y_Trafo(x0,y0,z0)
	z=DISK_Z_Trafo(x0,y0,z0)
	return 0.34*np.exp(-((np.sqrt(x**2 + (3.1*y)**2) - XD)/LD)**4)*np.exp(-(z/HD2)**2)

#Ioniesierter Wasserstoff (Warm) 
def n_I_HII_WIM(x0,y0,z0):
	r=np.sqrt(x0**2+y0**2)
	P1 = np.exp(-(x0**2+(y0-y3)**2)/L3**2)*np.exp(-(z0-z3)**2/H3**2)
	P2 = np.exp(-((r-L2)/(0.5*L2))**2)*sech2(z0/H2)
	P3 = np.cos(np.pi*r*0.5/L1)*sech2(z0/H1)
	return 8.0*(P1+0.009*P2+0.005*P3)

#Ionieserter Wasserstoff (very hot)
def n_I_HII_VHIM(x0,y0,z0): 
	e = y0*cos(alphaVH)+z0*sin(alphaVH)
	s = -y0*sin(alphaVH) + z0*cos(alphaVH)
	return 0.29*np.exp(-((x0**2+e**2)/LVH**2 + s**2/HVH**2))

def n_I_HII(x0,y0,z0):	# Ioniesert insg.
	return n_I_HII_VHIM(x0,y0,z0) +n_I_HII_WIM(x0,y0,z0)
def n_I_HI(x,y,z):		#Atomar insg.
	return n_I_HI_DISK(x,y,z) + n_I_HI_CMZ(x,y,z)
def n_I_H2(x,y,z):		# Molekular insg.
	return n_I_H2_CMZ(x,y,z) + n_I_H2_DISK(x,y,z)

"""Zusammenfassen"""

def n_HI(x,y,z):
	R = np.sqrt(x**2+y**2)
	if R<=3000:
		return n_I_HI(x,y,z)
	else:
		return n_A_HI(R,z)
def n_HII(x,y,z):
    R=np.sqrt(x**2+y**2)
    if R<=3000:
        return n_I_HII(x,y,z)
    else:
        return n_A_HII(R,z)
def n_H2(x,y,z):
	R=np.sqrt(x**2+y**2)
	if R<=3000:
		return n_I_H2(x,y,z)
	else:
		return n_A_H2(R,z)

