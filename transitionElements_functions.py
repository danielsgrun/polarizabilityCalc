# Calculate the Einstein A coefficient and the dipole matrix element 
# Written by: Daniel S. Grun #
# IQOQI Innsbruck, 2023 #

import numpy as np
from scipy import constants as const

def squareDipole(J2, nu12, Aik):
    g2 = 2*J2+1
    u12 = g2/(2*np.pi*nu12)**3 * 3*np.pi*const.epsilon_0*const.hbar*const.c**3 * Aik
    
    return u12

if __name__ == "__main__":
    pass
