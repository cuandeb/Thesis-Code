import matplotlib.pyplot as plt 
import json 
import numpy as np 
import os

#-----------------------------

#sets the style and fontsize for the light curve plots
font = {'family' : 'serif',
        'size'   : 12}
plt.rc('font', **font)


TWO_PI = 2*np.pi
#----------------------------

#data extraction stpdf 

def extract_json(filepath):
    with open(filepath, 'r') as f:
        data = json.load(f)
    return data


def extract_source_id(filepath, delim = '/'):
    filename = filepath.split(delim)[-1]
    return filename.split('.')[0]


def check_directory(path): 
    isExist = os.path.exists(path)
    if not isExist:  
  # Create a new directory because it does not exist 
        os.makedirs(path)
        print("New directory is created!")


def source_data(condition, filepath, input_data): 
    if condition: 
        data = input_data

    else: 
        data = extract_json(filepath)
    return data
#---------------------------

#Plotting functions

def generate_plots(filepath, input_data = None, no_file = True, mission = 'DR3'):

    data = source_data(no_file, filepath, input_data)

    _name = extract_source_id(filepath) 
    
    folder = r"C:\Users\cuand\ESAC Internship\OHIR sample\Source Plots"

    _jd = np.array(data['time'])
    path = os.path.join(folder, _name)

    check_directory(path)
    
    _period = data['final_period'][0]
    _err_period = data['final_period'][1]
    _phase = data['phase'][0]
    _ffac = data['ffactor'][0]
    _mag_med = data['mag_med'][0]
    _amp = data['amplitude'][0]
    _jd_min = data['jd_min']  
    
    
    chi_squared_plot(data, path, _name, _period)
    mag_jd_plot(data, path, _name, mission)
    _x,_y,_mi,_ma = plot_alc_model(data, path,
                   _name, 
                   _period, 
                   _err_period, 
                   _phase,
                   _amp,
                   _mag_med, 
                   _ffac,
                   _jd_min,
                   mission)
    folded_alc(data, path, 
                _name, 
                _jd, 
                _period, 
                _phase, 
                _mag_med,
                _jd_min, 
                _amp,
                _ffac,
                _x, _y, _mi, _ma)
    
    




def chi_squared_plot(data, directory, name, period):
    fig, ax = plt.subplots(1, figsize=(9,7))

    ax.plot(data['period'], data['chi'],'k-')
    ax.set_xlabel('Period')
    ax.set_ylabel(r'$\chi^2$')
    ax.set_xlim(0,2500)
    ax.axvline(period,c='r', linestyle='--')
    ax.set_title(f"Gaia-{name} Periodogram (P(d)={period})")
    plt.savefig(f'{directory}/{name}_periodogram.pdf', format='pdf')
    plt.clf()
    
    
    
def mag_jd_plot(data, directory, name, mission):
                
    fig, ax = plt.subplots(1, figsize=(9,7))

    ax.plot(data['time'], data['mag'][0],'gs',
            markerfacecolor='none',label=f'Gaia {mission}')
    ax.errorbar(data['time'], data['mag'][0],
                yerr=data['mag'][1],fmt='none',
                color='k', capsize=2)
    ax.set_xlabel('Julian Day (+2450000)')
    ax.set_ylabel(r'Mag (G-Band)')
    
    ax.set_title(f"Gaia-{name} G-Mag")
    ax.set_xlim(min(data['time'])-500, max(data['time'])+500)
    ax.invert_yaxis()
    
    ax.legend(loc='best') 
    
    plt.savefig(f'{directory}/{name}_mag.pdf', format='pdf')
    plt.clf()
    

def plot_alc_model(data, directory,
                   name, 
                   period, 
                   err_period, 
                   phase,
                   amp,
                   mag_med, 
                   ffac,
                   jd_min,
                   mission
                  ):
    
    jd_range = int(np.max(data['time'])-np.min(data['time'])+data['final_period'][0])
    
    x = np.ndarray(jd_range).astype(float)
    y = np.ndarray(jd_range).astype(float)
    
    mi = np.ndarray(2).astype(float)
    ma = np.ndarray(2).astype(float)
    
    
    for i in range(0,jd_range): 
        x[i]=i
        xf=((x[i]/period)+phase+0.5)-int((x[i]/period)+phase+0.5)
        
        if xf < ffac: 
            omega= xf/(2*ffac)
            y[i]=mag_med + (amp/2)*np.cos(TWO_PI*omega)
        else: 
            omega = ((xf-1)/(2*(1-ffac)))+1
            y[i] = mag_med + (amp/2)*np.cos(TWO_PI*omega)
            
    mi[0]= np.min(y)
    mi[1]= np.min(np.array(data['mag'][0])-np.array(data['mag'][1]))
    
    ma[0]= np.max(y)
    ma[1]= np.max(np.array(data['mag'][0])-np.array(data['mag'][1]))
    
    x = x + jd_min-data['final_period'][0]/2
    
    
    #estimation of the JD for the minimum of the LC
    
    minimum = 150
    for k in range(int(data['final_period'][0])): 
        if y[k] <= minimum: 
            minimum = y[k]
            jd_maximum = x[k]
            
            if jd_maximum <= jd_min: 
                jd_maximum = jd_maximum+data['final_period'][0]
                err_jdmax = np.abs((jd_min-jd_maximum)/(period*err_period))
                    
    fig, ax = plt.subplots(1, figsize=(10,7))


    ax.errorbar(data['time'], data['mag'][0],yerr=data['mag'][1],fmt='none',color='k', capsize=2)
    ax.plot(data['time'], data['mag'][0],'gs',markerfacecolor='none', label=f'Gaia {mission}')
    ax.plot(x,y,'r-', label = 'Fitted Model')
    ax.set_xlabel('Julian Day (+2450000)')
    ax.set_ylabel(r'Mag (G-Band)')
    ax.set_title(f"Gaia-{name} G-Mag")
    ax.set_xlim(min(data['time'])-500, max(data['time'])+500)
    ax.set_ylim(np.max(ma)+(np.max(ma)-np.min(mi))*.1,np.min(mi)-(np.max(ma)-np.min(mi))*.1)
    ax.legend(loc='best')

    plt.savefig(f'{directory}/{name}_alc_model.pdf', format='pdf')
    plt.clf()

    return x, y, mi, ma



def folded_alc(data, directory, 
                name, 
                jd, 
                period, 
                phase, 
                mag_med,
                jd_min, 
                amp,
                ffac,
                x, y, mi, ma):


    x=np.zeros(2000)
    y=np.zeros(2000)
    yz=np.zeros(2000)

    for i in range(2000):
        x[i]=i/1000.
        if x[i] < ffac:
            ffactor=x[i]/(2*ffac)
            yz[i]=mag_med+(amp/2)*np.cos(ffactor*TWO_PI)
        elif x[i] >= ffac and x[i] < 1:
            ffactor=((x[i]-1)/(2*(1-ffac)))+1
            yz[i]=mag_med+(amp/2)*np.cos(ffactor*TWO_PI)

        else:
            yz[i]=yz[i-1000]


    for z in range(1999):
        if z <= 1749:
            y[z]=yz[z+250]
        else:
            y[z]=yz[z-1750]


    jdf=((jd-jd_min)/period+phase-.25)-((jd-jd_min)/period+phase-.25).astype(int) 

    fig, ax = plt.subplots(1, figsize=(9,7))



    ax.plot(jdf, data['mag'][0],'gs',markerfacecolor='none')
    ax.plot(jdf+1, data['mag'][0],'gs',markerfacecolor='none')
    ax.plot(x,y,'r-')
    ax.errorbar(jdf, data['mag'][0],yerr=data['mag'][1],fmt='none',color='k', capsize=2)
    ax.errorbar(jdf+1, data['mag'][0],yerr=data['mag'][1],fmt='none',color='k', capsize=2)
    ax.set_xlabel('Phase')
    ax.set_ylabel(r'Mag (G-Band)')
    ax.set_title(f"{name} G-Mag Phase")
    ax.set_xlim(0,2)
    ax.set_ylim(np.max(ma)+(np.max(ma)-np.min(mi))*.1,np.min(mi)-(np.max(ma)-np.min(mi))*.1)


    plt.savefig(f'{directory}/{name}_folded_phase.pdf', format='pdf')
    plt.clf()
    

#------------------------

TEST_FILE = r"C:\Users\cuand\ESAC Internship\IDL_Code\OHIR_sources_json\4596654564903027456.json"


if __name__=="__main__":
    generate_plots(TEST_FILE)