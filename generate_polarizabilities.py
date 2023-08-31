# Generate individual polarizabilities
# Written by: Daniel S. Grun #
# IQOQI Innsbruck, 2023 #

import numpy as np
import matplotlib.pyplot as plt

from reading_functions import extract_lines
from polarizability_functions import alphaK
from scipy import constants as const
from tqdm import tqdm

a0 = 5.29e-11 #Bohr radius in m
atomic_unit = const.e**2 * a0**4 * const.m_e / (const.hbar**2)

def constants_svt(J):
    constant_s = -1/np.sqrt(3*(2*J+1)) / atomic_unit
    constant_v = np.sqrt(2*J/((J+1)*(2*J+1))) / atomic_unit
    constant_t = np.sqrt(2*J*(2*J-1)/(3*(J+1)*(2*J+1)*(2*J+3))) / atomic_unit
    
    return [constant_s, constant_v, constant_t]


def generate_polarizabilities(E1,E2,trapWavelengths,J1,J2,Aik,constants, verbose=True):
    constant_s, constant_v, constant_t = constants
    if verbose:
        loop_range = tqdm(range(len(trapWavelengths)))
    else:
        loop_range = range(len(trapWavelengths))
        
    alpha_s, alpha_v, alpha_t = [],[],[]    
        
    for i in loop_range:    
        wl = trapWavelengths[i]
        alpha_s.append(alphaK(E1, E2, wl, J1, J2, Aik, K=0)*constant_s)
        alpha_v.append(alphaK(E1, E2, wl, J1, J2, Aik, K=1)*constant_v)
        alpha_t.append(alphaK(E1, E2, wl, J1, J2, Aik, K=2)*constant_t)
        
    alpha_s, alpha_v, alpha_t = np.array(alpha_s), np.array(alpha_v), np.array(alpha_t)
    
    return np.array([alpha_s, alpha_v, alpha_t])
  
  
def organized_polarizabilities(trapWavelengths, stateLabels=['GS'], atom='Er', 
                               use_latest_GS=False, savefiles=True):
    wl = trapWavelengths
    if use_latest_GS:
        alpha_real, alpha_imag = [],[]
        wls, Aik, [E1, J1], [E2, J2] = extract_lines(atom, stateLabels[0], use_latest_GS=use_latest_GS)
        constants = constants_svt(J1[0])
        alphas = generate_polarizabilities(E1, E2, wl, J1, J2, Aik, constants)
        alpha_real = np.real(alphas)
        alpha_imag = np.imag(alphas)
        if savefiles:
            print("Saving files")
            directory = 'calculated_polarizabilities/'
            np.savetxt(directory+'GS_real_latestGS.txt', np.vstack([wl, alpha_real]))
            np.savetxt(directory+'GS_imag_latestGS.txt', np.vstack([wl, alpha_imag]))
    else:
        alpha_real, alpha_imag = [],[]
        for i in range(len(stateLabels)):    
            wls, Aik, [E1, J1], [E2, J2] = extract_lines(atom, stateLabels[i])
            constants = constants_svt(J1[0])
            alphas = generate_polarizabilities(E1, E2, wl, J1, J2, Aik, constants)
            print("Generating total polarizability vector")
            alpha_real.append(np.real(alphas))
            alpha_imag.append(np.imag(alphas))
            if savefiles:
                print("Saving files for {}".format(stateLabels[i]))
                directory = 'calculated_polarizabilities/'
                name = stateLabels[i]
                # re_alpha = np.reshape(alpha_real[name], (len(wl),3))
                # im_alpha = np.reshape(alpha_imag[name], (len(wl),3))
                # re_alpha = np.append(re_alpha, wl, axis=0)
                # im_alpha = np.append(im_alpha, wl, axis=0)
            
                # np.savetxt(directory+name+'_real.txt', np.c_[wl, re_alpha])
                # np.savetxt(directory+name+'_imag.txt', np.c_[wl, im_alpha])
                np.savetxt(directory+name+'_real.txt', np.vstack([wl, alpha_real[i]]))
                np.savetxt(directory+name+'_imag.txt', np.vstack([wl, alpha_imag[i]]))

    return [alpha_real, alpha_imag]
    
        

if __name__ == "__main__":
    
    atom = "Er"
    wl = np.linspace(300,1300,10000)*1e-9
    stateLabels = ['GS', '583nm', '841nm']
    
    alphas = organized_polarizabilities(wl, stateLabels, use_latest_GS = False, savefiles=False)
    
    plt.plot(1e9*wl, alphas['GS'][0], color='b', alpha=0.8, label='GS')
    plt.plot(1e9*wl, alphas['583nm'][0], color='r', alpha=0.8, label='583nm')
    plt.plot(1e9*wl, alphas['841nm'][0], color='g', alpha=0.8, label='841nm')
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Scalar polarizability (a.u.)")
    # plt.ylim(-200,600)
    plt.legend(loc=0)