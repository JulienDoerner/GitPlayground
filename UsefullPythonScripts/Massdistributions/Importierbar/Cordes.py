import numpy as np

#Konstanten 
ne1 = 0.025
ne2 = 0.2
H1 = 1.0
H2 = 0.15
A1 = 20
A2 = 2
R2 = 4.

def n_HII(R,z):
    R=abs(R)
    P1 = ne1 * np.exp(-abs(z)/H1)*np.exp(-(R/A1)**2)
    P2 = ne2 * np.exp(-abs(z)/H2)*np.exp(-((R-R2)/A2)**2)
    return P1+P2



