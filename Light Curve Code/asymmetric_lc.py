import numpy as np
from numpy.linalg import inv
from astropy.io import ascii
import json
import os
import time

#-------------------------------

#Define functions

def files_from_directory(directory):
    
    '''Takes the files from an input directory
       and returns the file names as strings 
       in a list'''
       
    for filename in os.listdir(directory):
        if filename.endswith('.dat'):
            files=[filename for filename in os.listdir(directory)]            
    return files     


def read_ascii_file(filepath): 
    data = ascii.read(filepath)
    try:    
        return np.array(data['JD']), np.array(data['g_transit_mag']),np.array(data['mag_error'])
    except ValueError: 
        print('Incompatible Field. Table must have columns containing'
        'JD, g_transit_mag, mag_error')


def generate_params(filepath):

    """Extracts necessary params from ASCII table"""

    _jd, _mag, _mag_err = read_ascii_file(filepath)
    _jd = _jd-2450000
    _nobs = len(_jd)
    _weights = 1/(_mag_err**2)
    _jd_min = np.min(_jd)
    _jd = _jd - _jd_min
    return _nobs, _weights, _jd, _mag, _jd_min, _mag_err


def generate_arrays(n_obs,n_points): 
    """Generates the arrays used by the chi_squared function """
    _chi = np.zeros(n_points,dtype=np.int64)
    _per = np.zeros(n_points,dtype=np.int64)
    _func = np.zeros(n_obs,dtype=np.int64)

    return _chi, _per, _func

def extract_source_id(filepath, delim = '/'):
    filename = filepath.split(delim)[-1]
    return filename.split('.')[0]


def write_to_json(dictionary,filepath):
    
    """Dumps the results of a curve fit to a JSON file"""
    
    folder = r"C:\Users\cuand\ESAC Internship\OHIR Sample\OHIR_sources_json"
    source = extract_source_id(filepath)
    new_file = f"{folder}\{source}.json"
    
    with open(new_file, "w") as outfile:
        json.dump(dictionary, outfile, indent=4) 


def list_to_file(_list, _directory): 
    return [f"{_directory}/{i}.dat" for i in _list]  

def read_list_from_file(file): 
    data = np.loadtxt(file,dtype='str')
    return data[:,0], [float(i) for i in data[:,1]]


#------------------------------------------
# Chi squared function

def run_chi_squared(filepath,period=None):

    TWO_PI = 2*np.pi

    if period is not None: 
        P_MIN=period-99
        P_MAX=period+99
        DELTA_P = int(P_MAX-P_MIN)
        print("Short Run")
    
    #else: 
    P_MIN = 100
    P_MAX = 2500
    DELTA_P = P_MAX - P_MIN
    PERIODS = np.arange(P_MIN,P_MAX)
    print("Long Run")

    nobs, weights, jd, mag, jd_min, mag_err = generate_params(filepath)

    chi, per, func = generate_arrays(nobs, DELTA_P)
    chi_min = 10000000.

    for i in PERIODS:
        
        chiminperiod = 10000000.

        for j in np.arange(1,100, dtype=np.float64):
            f=j/100.
            jd_int=np.int32((jd/per[i])+f)
            jdr=((jd/per[i])+f)-jd_int

            for k in np.arange(10,90, dtype=np.float64): 
                fact=0.01*k
                for m in np.arange(nobs):
                    if jdr[m] < fact:
                        func[m]=jdr[m]/(2*fact)
                        
                    else:
                        func[m]=((jdr[m]-1)/(2*(1-fact)))+1
                        
                    

                cc=np.sum(((np.cos(TWO_PI*func))**2)*weights)
                co=np.sum((np.cos(TWO_PI*func))*weights)
                err=np.sum(weights)
                com=np.sum(np.cos(TWO_PI*func)*mag*weights)
                m = np.sum(mag*weights)

                    #coefficients

                det = cc*err - co*co 
                b=(com*err-co*m)/det
                a=(cc*m-com*co)/det

                    #fit error 

                chi_phase = np.sum((weights*((mag-a-b*np.cos(TWO_PI*func))**2))/(nobs-2))

                if chi_phase <= chiminperiod:
                    chi[i]=chi_phase
                    chiminperiod=chi_phase

                if chi_phase <= chi_min: 
                    chi_min = chi_phase
                    ffactor=fact
                    error_ffactor=0.01
                    period = per[i]
                    error_period=(period**2)/(2*(np.max(jd)-np.min(jd)))
                    phase=f
                    error_phase=.01
                    matrix=np.array([[cc,co],[co,err]])                
                    matrix_inv=inv(matrix)
                    ampli=b*2
                    if ampli <= 0:
                        ampli=-ampli
                        ffactor=1-ffactor
                        phase=phase+ffactor

                    if phase > 1.:
                        phase=phase-1

                    error_ampli=2*np.sqrt(matrix_inv[0,0])
                    magmed=a
                    error_magmed=np.sqrt(matrix_inv[1,1])
                    #print(ffactor)

    jd = jd+jd_min

    results = {
        'time': jd.tolist(),
        'jd_min':jd_min,
        'chi': chi.tolist(),
        'mag':(mag.tolist(), mag_err.tolist()),
        'period':per.tolist(),
        'final_period':(period, error_period),
        'phase':(phase, error_phase),
        'amplitude':(ampli, error_ampli),
        'ffactor':(ffactor,error_ffactor), 
        'mag_med': (magmed, error_magmed)
    }

    return results



#-----------------------------
#Testing the script

#DIR = <Insert path of ascii files here>


def asymmetric_fit(filepath, period=None): 
    #print(period)   
    data = run_chi_squared(filepath, period)
    write_to_json(data, filepath)
    return data

    
      
if __name__=='__main__': 
    #print('Start script...')
    file_list = files_from_directory(DIR)
    for file in file_list:
        asymmetric_fit(file)

    

   

