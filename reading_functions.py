## Reading and preparing the file containing spectral data
# Written by: Daniel S. Grun #
# IQOQI Innsbruck, 2023 #

import numpy as np
import pandas as pd

def extract_lines(which_atom, stateLabel='GS', use_latest_GS = False):
    
    if use_latest_GS:
        print("Hi")
        data = np.loadtxt('lineEr_latest.txt', unpack=True)
        Aki = data[6]
        wavelengths = data[5]*1e-9
        E1, E2 = data[0], data[2]
        J1, J2 = data[1], data[3]
  
    else:
        if which_atom == 'Er':
            df = pd.read_excel('Erbium_lines.xlsx') # Read the lines data from file
        elif which_atom == 'Dy':
            df = pd.read_excel('Dysprosium_lines.xlsx')
        
        if stateLabel=='GS':
            df = df[df['En_low'] == 0]
        else:
            df_label = df[df['Label']== stateLabel]
            Energy = df_label['En_up']
            df = df[df['En_low'] == float(Energy)]
  
        wavelengths = df['Wavelength']
        Aki = df['Aki']
        E1, J1 = df['En_low'], df['J_low']
        E2, J2 = df['En_up'], df['J_up']
    
        wavelengths = wavelengths.to_numpy() * 1e-9
        Aki = np.float32(Aki.to_numpy())
        E1, J1 = E1.to_numpy(), J1.to_numpy()
        E2, J2 = E2.to_numpy(), J2.to_numpy()  
  
    return [wavelengths, Aki, [E1, J1], [E2, J2]]


def read_polarizabilities(use_Maxence=False, use_latest_GS=False):
    directory='calculated_polarizabilities/'
    alphas_real, alphas_imag = {},{}

    
    if use_Maxence:  
        stateLabels_real = ['GS', '401nm', '583nm', '841nm', '1299nm']
        stateLabels_imag = ['GS', '1299nm']
        
        for name in stateLabels_real:
            alphas_real[name] = np.loadtxt(directory+'Max_'+name+'_real.txt', unpack=True)
            alphas_real[name] = alphas_real[name][1:4]
            
        for name in stateLabels_imag:
            alphas_imag[name] = np.loadtxt(directory+'Max_'+name+'_imag.txt', unpack=True)
            alphas_imag[name] = alphas_imag[name][1:4]

        Energies = np.loadtxt(directory+'Max_'+name+'_real.txt', unpack=True, usecols=0)
        wavelengths = 1e-2/Energies
        print(Energies)
        
    elif use_latest_GS:
        name='GS'
        alphas_real[name] = np.loadtxt(directory+name+'_real_latestGS.txt', unpack=False)
        alphas_imag[name] = np.loadtxt(directory+name+'_imag_latestGS.txt', unpack=False)
        
        alphas_real[name] = alphas_real[name][1:]
        alphas_imag[name] = alphas_imag[name][1:]        
   
        wavelengths = np.loadtxt(directory+name+'_real_latestGS.txt', unpack=False)
        wavelengths = wavelengths[0]
        
    else:
        stateLabels = ['GS', '583nm', '841nm']
        for name in stateLabels:
            alphas_real[name] = np.loadtxt(directory+name+'_real.txt', unpack=False)
            alphas_imag[name] = np.loadtxt(directory+name+'_imag.txt', unpack=False)
            
            alphas_real[name] = alphas_real[name][1:]
            alphas_imag[name] = alphas_imag[name][1:]
            
        wavelengths = np.loadtxt(directory+name+'_real.txt', unpack=False)
        wavelengths = wavelengths[0]
            
    return [wavelengths, alphas_real, alphas_imag]


# if __name__ == "__main__":
#     which_atom = "Er"
#     stateLabel = '841nm'
#     wls, Aki, [E1, J1], [E2, J2] = extract_lines(which_atom, stateLabel)
  

