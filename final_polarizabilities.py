# Obtaining and plotting the final results for polarizabilities #
# Written by: Daniel S. Grun #
# IQOQI Innsbruck, 2023 #

from generate_polarizabilities import organized_polarizabilities
from final_polarizabilities_functions import organized_combined_alphas, generate_polarizations
from reading_functions import read_polarizabilities

import numpy as np
import matplotlib.pyplot as plt

def lambda_unique(wavelengths, lambda_0):
  
  diff = abs(wavelengths-lambda_0)  
  min_diff = min(abs(wavelengths-lambda_0))
  loc = np.where(diff == min_diff)[0][0]

  return loc

if __name__ == "__main__":
    atom = "Er"
    
    plot_versus_wavelength = True
    
    if plot_versus_wavelength:
    
        wl, alphas_real, alphas_imag = read_polarizabilities()
        wl_Max, alphas_real_Max, alphas_imag_Max = read_polarizabilities(use_Maxence=True)
        wl_late, alphas_real_late, alphas_imag_late = read_polarizabilities(use_latest_GS=True)
        
        stateLabels = ['GS']
        polarization=[[1,0]]
        B_field = [0,0,1]
        fine_structure_dict = {'GS':[6,-6], '583nm':[7,-7], '841nm':[7,-7]}
    
        alpha_T_real = organized_combined_alphas(polarization, wl, alphas_real, 
                                                fine_structure_dict, stateLabels=stateLabels, 
                                                B_field=B_field)

        alpha_T_real_Max = organized_combined_alphas(polarization, wl_Max, alphas_real_Max, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                B_field=B_field)
        
        alpha_T_real_late = organized_combined_alphas(polarization, wl_late, alphas_real_late, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                B_field=B_field)
        
        alpha_T_imag = organized_combined_alphas(polarization, wl, alphas_imag, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                B_field=B_field)

        alpha_T_imag_Max = organized_combined_alphas(polarization, wl_Max, alphas_imag_Max, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                B_field=B_field)
        
        alpha_T_imag_late = organized_combined_alphas(polarization, wl_late, alphas_imag_late, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                B_field=B_field)
    
        plt.figure()
        plt.plot(1e9*wl, alpha_T_real['GS'][0], lw=2, alpha=0.6, color='k', label='Daniel')
        plt.plot(1e9*wl_Max, alpha_T_real_Max['GS'][0], lw=2, alpha=0.6, color='g', label='Maxence')
        plt.plot(1e9*wl_late, alpha_T_real_late['GS'][0], lw=2, alpha=0.6, color='r', label='Latest')
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Re$(\\alpha)$ (a.u.)")
        # plt.xlim(500,540)
        # plt.ylim(0, 500)
        # plt.yticks(np.arange(0,501,50))
        # plt.yscale('log')
        plt.legend(loc=0)

        plt.figure()
        plt.plot(1e9*wl, alpha_T_imag['GS'][0], lw=2, alpha=0.6, color='k', label='Daniel')
        plt.plot(1e9*wl_Max, alpha_T_imag_Max['GS'][0], lw=2, alpha=0.6, color='g', label='Maxence')
        plt.plot(1e9*wl_late, alpha_T_imag_late['GS'][0], lw=2, alpha=0.6, color='r', label='Latest')
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Im$(\\alpha)$ (a.u.)")
        # plt.xlim(480,500)
        # plt.ylim(1e-5,5e-5) 
        plt.yscale('log')               
        plt.legend(loc=0)
    
    else:
        
        wl, alphas_real, alphas_imag = read_polarizabilities()
        wl_Max, alphas_real_Max, alphas_imag_Max = read_polarizabilities(use_Maxence=True)
        wl_late, alphas_real_late, alphas_imag_late = read_polarizabilities(use_latest_GS=True)
        
        wl_0 = 488e-9 # fixed trap wavelength
        pos, pos_Max, pos_late = lambda_unique(wl, wl_0), lambda_unique(wl_Max, wl_0), lambda_unique(wl_late, wl_0)
        
        theta, polarization = generate_polarizations()
        
        stateLabels = ['GS']
        
        fine_structure_dict = {'GS':[6,-6], '583nm':[7,-7], '841nm':[7,-7]}
    
        alpha_T_real = organized_combined_alphas(polarization, [1], alphas_real, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                index=[pos])

        alpha_T_real_Max = organized_combined_alphas(polarization, [1], alphas_real_Max, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                index=[pos_Max])
        
        alpha_T_real_late = organized_combined_alphas(polarization, [1], alphas_real_late, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                index=[pos_late])
        
        alpha_T_imag = organized_combined_alphas(polarization, [1], alphas_imag, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                index=[pos])

        alpha_T_imag_Max = organized_combined_alphas(polarization, [1], alphas_imag_Max, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                index=[pos_Max])
        
        alpha_T_imag_late = organized_combined_alphas(polarization, [1], alphas_imag_late, 
                                                fine_structure_dict, stateLabels=stateLabels,
                                                index=[pos_late])
        
        plt.figure()
        plt.plot(theta, alpha_T_real['GS'][:,0], lw=2, alpha=0.6, color='k', label='Daniel')
        plt.plot(theta, alpha_T_real_Max['GS'][:,0], lw=2, alpha=0.6, color='g', label='Max')
        plt.plot(theta, alpha_T_real_late['GS'][:,0], lw=2, alpha=0.6, color='r', label='Latest lines')
        plt.xticks(ticks=np.arange(-np.pi, np.pi+0.1, np.pi/2),
                   labels=['-$\pi$', '-$\pi/2$', '$0$', '$\pi/2$', '$\pi$'])
        plt.xlabel("Ellipticity angle (rad)")
        plt.ylabel("Re($\\alpha$) (a.u.)")
        plt.legend(loc=0)
        
        plt.figure()
        plt.plot(theta, alpha_T_imag['GS'][:,0], lw=2, alpha=0.6, color='k', label='Daniel')
        plt.plot(theta, alpha_T_imag_Max['GS'][:,0], lw=2, alpha=0.6, color='g', label='Max')
        plt.plot(theta, alpha_T_imag_late['GS'][:,0], lw=2, alpha=0.6, color='r', label='Latest lines')
        plt.xticks(ticks=np.arange(-np.pi, np.pi+0.1, np.pi/2),
                   labels=['-$\pi$', '-$\pi/2$', '$0$', '$\pi/2$', '$\pi$'])
        plt.xlabel("Ellipticity angle (rad)")
        plt.ylabel("Im($\\alpha$) (a.u.)")
        plt.legend(loc=0)