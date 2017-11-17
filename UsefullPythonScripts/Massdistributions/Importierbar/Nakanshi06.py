import numpy as np

def n_H2(x,y,z):     #Eingabe in kpc
    r=np.sqrt(x**2+y**2)/0.85
    E1=11.2*np.exp(-r**2/0.874)
    E2=0.83*np.exp(-((r-4)/3.2)**2)
    h=1.06*10**-3*(10.8*np.exp(0.28*r)+42.78)
    return 0.94*(E1+E2)*np.exp(-np.log(2)*(z/h)**2)
