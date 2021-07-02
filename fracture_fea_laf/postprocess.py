import numpy as np

def celsius2kelvin(temperature):
    return temperature + 273.15

def k2j(k, E, nu, plane_stress=False):
    """
    Convert fracture 
    
    Parameters
    ----------
    k: float
    E: float
        Young's modulus in GPa.
    nu: float
        Poisson's ratio
    plane_stress: bool
        True for plane stress (default) or False for plane strain condition.
    
    Returns
    -------
    float
    
    """
    
    if plane_stress:
        E = E / (1 - nu ** 2)
        
    return k ** 2 / E


def j2k(j, E, nu, plane_stress=True):
    """
    Convert fracture 
    
    Parameters
    ----------
    j: float (in N/mm)
    E: float
        Young's modulus in GPa.
    nu: float
        Poisson's ratio
    plane_stress: bool
        True for plane stress (default) or False for plane strain condition.
    
    Returns
    -------
    K : float
        Units are MPa m^0.5.
    """
    
    if plane_stress:
        E = E / (1 - nu ** 2)
        
    return (j * E) ** 0.5

def calc_cpf(data):
    data_ordered = np.sort(data)
    u =  (np.arange(len(data)) + 1 - 0.3) / (len(data) + 0.4)
    
    return data_ordered, u

def find_closest(a, b, x, tol=0.01):
    """
    For two arrays `a` and `b`, find the elements in `b` closest to the `x` element in `a`.
    
    Parameters
    ----------
    a: narray of shape (N,)
    b: narray of shape (N,)
    x: float
    tol: float
        Tolerance value. 
    Returns
    -------
    """
    a_idx = np.isclose(a, x, atol=tol)
    a_i = a[a_idx]
    b_i = b[a_idx]

    return a_i, b_i
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def Fk(shear_stress, temperature, Tath=700, Uk=0.33 * 1.60217657e-19, sigmap=900*1e6):
    """
    Parameters
    ----------
    shear_stress : ndarray of shape (N,)
        CRSS or maximum shear stress (in Pa or Nm-2)
    temperature : float
        Temperature (in Kelvin)
    Tath : float
        The athermal temperature in Kelvin. Default set for iron = 700 K.
    Uk : float
        The activation energy (in Joules). Default set for iron.
    sigmap : float
        critical Peierls stress (in Pascals). Default set for iron.

    Returns
    -------
    float
        Kink-formation free energy (in Joules).

    Ref. T. D. Swinburne, S. L. Dudarev, Kink-limited Orowan strengthening explains the brittle to ductile transition of irradiated and unirradiated bcc metals. Phys. Rev. Mater. 2, 73608 (2018).
    """

    fk = Uk *(1 - temperature / Tath - (shear_stress / sigmap) / (1 - temperature / Tath))
    fkneg_idx = np.where(fk<0)[0]

    fk[fkneg_idx] = 0.0
    
    return fk

def lstar(shear_stress, temperature, b=0.248e-9):
    """
    Parameters
    -------
    shear_stress : ndarray of shape (N,)
        CRSS or maximum shear stress (in Pa or Nm-2)
    temperature : float
        Temperature (in Kelvin)
    b : float
        Burger's vector (in m). Default set for iron. 
        
    Returns
    -------
    ndarray of shape (N,)
        Critical lenght L* (in m).
    
    """
    # Boltzmann constant
    kb =  1.3806498e-23 # J/K

    return b * np.exp( Fk(shear_stress, temperature) / (kb * temperature) )


def thinning_function(obstacle_distance, shear_stress, temperature):
    """
    Thinning function applied in a modified local approach to cleavage fracture predictions.
    
    Parameters
    ----------
    obstacle_distance: array_like
        Average distance between dislocation obstacles [1].
    shear_stress: array_like
        Maximum shear stress (in Pa or Nm-2)
    temperature: float 
        Temperature in degree Kelvin.
    Returns
    -------
    thinf : 

    References
    ----------
    [1] T. D. Swinburne and S. L. Dudarev, “Kink-limited Orowan strengthening explains the brittle to ductile transition of irradiated and unirradiated bcc metals,” Phys. Rev. Mater., vol. 2, p. 73608, 2018.
    
    """
    # Boltzmann constant
    kb = 1.3806498e-23 # Joules/Kelvin
    

    
    idx_more = lstar(shear_stress, temperature) > np.ones(len(shear_stress)) * obstacle_distance
    idx_less = ~idx_more
    ref_temp = - 200 + 273.15
    norm_factor = 1 / (
        2 * Fk(np.array([0]), ref_temp) / (kb * ref_temp)
    )[0]

    fk = np.zeros(len(shear_stress))
    fk[idx_more] = 2 * Fk(shear_stress, temperature)[idx_more] / (kb * temperature)
    fk[idx_less] = Fk(shear_stress, temperature)[idx_less] / (kb * temperature)


    thinf = norm_factor * fk
        
    return thinf


def calc_sw_diff(sw_data, n0, ri, mi, gdata_all, skipped_lines_all):
    
    sw0 = sw_data[n0][mi][gdata_all[n0]['j0_idx'] - np.array(skipped_lines_all[n0])]
    swi = sw_data[ri][mi][gdata_all[ri]['j0_idx'] - np.array(skipped_lines_all[ri])]
    
    return (swi - sw0) * 100 / swi

def pred_j0(pred_i, sw_data, idxm, jskipped_lines_all, gdata_all, temps, tol=0.002):
    """
    Predict the characteristic fracture toughness J0 based on the Weibull stress

    Parameters:
    -----------
    pred_i : integer
        Index basis for prediction (zero-based)
    sw_data :  a list of length i of lists of length Ni
        Weibull stress data for i temperatures and Ni time increments for each temperature 
    idxm : integer
        Index of shape parameter m to be used
    jskipped_lines_all : list
        Lines to be skipped in the fracture toughness data for each temperature due to no Weibull data, e.g. when J is very small and there is limited plasticity.
    gdata_all :  list of dicts
       Fracture toughness data for each temperature. Here, the 'j0' data is used from the dictionary.
    temps : ndarray of size (i,)
        Temperatures (in degrees Celsius)
    """
    sw_pred = sw_data[pred_i][idxm][gdata_all[pred_i]['j0_idx']-jskipped_lines_all[pred_i]]
 
    print('Weibull stress: ', (sw_pred))
    print('gdata_all[pred_i]', gdata_all[pred_i]['j0'])
    for ri, sw_r in enumerate(sw_data):
        sw = np.array(sw_r[idxm])
        print('Temperature: ', str(temps[ri]))
        if ri != pred_i:
            j0_pred = np.where(np.isclose(sw, sw_pred, tol))
            shift_idx = len(gdata_all[ri]['gdata']['jintc']) - len(sw)
            print('shift_idx: ',shift_idx)
#             print('   j0_pred: ', j0_pred)
#             print(' sw: ', np.array(sw).shape)
#             print(type(j0_pred))
            jpred = gdata_all[ri]['gdata']['jintc'][j0_pred[0]+shift_idx]

            if len(jpred)>0:
                print('   jpred: {:4.2f}'.format(jpred.values[0]))
        else:
            print('basis for prediction')

        print('###')