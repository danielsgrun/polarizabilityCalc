# Total polarizability including field polarization and B-field orientation #
# Written by: Daniel S. Grun #
# IQOQI Innsbruck, 2023 #

import numpy as np  

def waveplate(wavep, theta):
    if wavep=='L2':
        Jones_waveplate = np.array([[np.cos(2*theta), np.sin(2*theta)],[np.sin(2*theta), -np.cos(2*theta)]])
    elif wavep=='L4':
        Jones_waveplate = np.array([[np.cos(theta)**2 + 1j*np.sin(theta)**2, (1-1j)/2*np.sin(2*theta)],
                                [(1-1j)/2*np.sin(2*theta), np.sin(theta)**2 + 1j*np.cos(theta)**2]])
    
    return Jones_waveplate  

def mirror(phi, theta_z):
    coord_rot = np.array([[np.cos(theta_z), np.sin(theta_z)],
                        [-np.sin(theta_z), np.cos(theta_z)]])  
    Jones_mirror = 1/np.sqrt(2) * np.array([[1,0],[0,np.exp(1j*phi)]])  

    return Jones_mirror @ coord_rot


def generate_polarizations():
    polarization_0 = [1,0]
    theta = np.linspace(-np.pi, np.pi, 100)
    polarizations = np.array([waveplate('L4', theta[i]) @ polarization_0 for i in range(len(theta))])
  
    return [theta,polarizations]

def combined_alpha(polarization,
                   pol_vector,
                   fine_structure,
                   B_field=[0,0,1]):
    
    (alpha_scalar,
     alpha_vectorial,
     alpha_tensorial) = (pol_vector[0],
                      pol_vector[1],
                      pol_vector[2])  
  
    E_x, E_y = polarization
    E_x, E_y = E_x * np.conj(E_x), E_y * np.conj(E_x)
    phase_diff = np.angle(E_y)
    E_x, E_y = abs(E_x), abs(E_y)
  
    sin_2chi = 2*E_x*E_y/(E_x**2+E_y**2)*np.sin(phase_diff)
  
    ellipticity = np.tan(0.5 * np.arcsin(sin_2chi))
                      
    B_field = np.array(B_field)
    B_field = B_field/np.sqrt(sum(abs(B_field)**2))
                    
    cos_plane = np.sqrt(sum(abs(B_field[:2])**2))/np.sqrt(sum(abs(B_field)**2))
    
    #cos_p = (B_field[0]*E_x + B_field[1]*E_y) * (1-abs(ellipticity)) * cos_plane
    cos_p_sq = (E_x*B_field[0])**2 + (E_y*B_field[1])**2 + 2*(E_x*E_y*B_field[0]*B_field[1])*np.sin(phase_diff)
    cos_p = np.sqrt(cos_p_sq) * cos_plane
    cos_k = B_field[-1]/np.sqrt(sum(abs(B_field)**2))
  
    print(ellipticity, cos_k, cos_p)
                      
    J, mJ = fine_structure
    # A = np.outer(polarization, np.conj(polarization))
    A = 2*np.imag(np.conj(polarization[0])*polarization[1])
  
    alpha_full = (alpha_scalar + 
                  cos_k*A*mJ/(2*J)*alpha_vectorial +
                  (3*mJ**2 - J*(J+1))/(J*(2*J+1))*(3*cos_p**2-1)/2*alpha_tensorial)
  
    return alpha_full

def organized_combined_alphas(polarization,
                              waveLength,
                              alpha_dict,
                              fine_structure_dict,
                              B_field=[0,0,1],
                              stateLabels=['GS'],
                              index='all'):
    alpha_full = {}
    
    if index=='all':
        for name in stateLabels:
            fine_structure = fine_structure_dict[name]
            print(fine_structure)
            alpha_full[name] = np.array([[combined_alpha(polarization[j], alpha_dict[name][:,i], fine_structure, B_field)
                                          for i in range(len(waveLength))]
                                         for j in range(len(polarization))])
    else:
        for name in stateLabels:
            fine_structure = fine_structure_dict[name]
            print(fine_structure)
            alpha_full[name] = np.array([[combined_alpha(polarization[j], alpha_dict[name][:,index[i]], fine_structure, B_field)
                                          for i in range(len(index))]
                                         for j in range(len(polarization))])
    return alpha_full
    
    