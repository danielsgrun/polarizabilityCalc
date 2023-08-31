# Polarizabilities functions #
# Written by: Daniel S. Grun #
# IQOQI Innsbruck, 2023 #

import numpy as np
from scipy import constants as const
from transitionElements_functions import squareDipole
from sympy.physics.wigner import wigner_6j


def alpha_JJp(u12,E1, E2, wl, J1, J2, Aik, K=0):
    omega12 = abs(E1-E2)/const.hbar
    omega = 2*np.pi*const.c/wl
    gamma1 = 0 # if ground-state polarizability: 0
    gamma2 = Aik
    gammaComp = gamma1+gamma2
    # gammaComp = gammaComp
    
    alpha_single = ( (-1)**(J1+J2+K)*np.sqrt(2*K+1)
                    *float(wigner_6j(1,K,1,J1,J2,J1))
                    *u12/const.hbar
                    *(1/(omega12-omega-1j*gammaComp/2)
                      +(-1)**K/(omega12+omega-1j*gammaComp/2)) )
    
    # print(alpha_single)
    
    return alpha_single
    

def alphaK(E1,E2,wl,J1,J2,Aik,K=0):
    E1, E2 = E1*100*const.h*const.c, E2*100*const.h*const.c
    omega12 = abs(E1-E2)/const.hbar
    nu12 = omega12/2/np.pi
    u12 = squareDipole(J2, nu12, Aik)
    alphaK = []
    for i in range(len(J2)):
        alphaK_single = alpha_JJp(u12[i], E1[i], E2[i], wl, J1[i], J2[i], Aik[i], K) 
        alphaK.append(alphaK_single)
    
    return sum(alphaK)         
        