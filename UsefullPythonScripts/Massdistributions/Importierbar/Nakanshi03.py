import numpy as np

def n_HI(x,y,z):         #Eingabe in kpc
    R=np.sqrt(x**2+y**2)
    exp1 = np.exp(-R/2.4)
    exp2 = np.exp(-((R-9.5)/4.8)**2)
    Plane = 0.94*(0.6*exp1+0.24*exp2)
    scaleheight = 1.06*(116.3+19.3*R+4.1*R**2-0.05*R**3)
    return Plane*np.exp(np.log(2)*(-(z/scaleheight)**2))
