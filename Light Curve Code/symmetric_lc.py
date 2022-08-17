import numpy as np 
from asymmetric_lc import generate_params
from numpy.linalg import inv
from tqdm import trange



def symmetric_fit(filepath): 

    P_MIN=100
    P_MAX=2500
    DELTA_P = P_MAX-P_MIN
    two_pi = 2*np.pi




    nobs, weights, jd, mag, jd_min, mag_err = generate_params(filepath)

    chi_min = 100000
    chi=np.zeros(DELTA_P)
    per=np.zeros(DELTA_P)

    for j in trange(DELTA_P): 
        per[j]=P_MIN+j

        pr=two_pi/per[j]
        jdr=jd*pr

        #Perform the summations

        sc = np.sum((np.sin(jdr)**2)*weights)
        sco = np.sum(np.sin(jdr)*np.cos(jdr)*weights)
        s = np.sum(np.sin(jdr)*weights)
        cc = np.sum((np.cos(jdr)**2)*weights)
        co = np.sum(np.cos(jdr)*weights)
        err = np.sum(weights)
        sm = np.sum(np.sin(jdr)*mag*weights)
        com = np.sum(np.cos(jdr)*mag*weights)
        m = np.sum(mag*weights)

        #calculate the Cramer coefficients 

        det=sc*cc*err+sco*co*s+s*sco*co-sc*co*co-sco*sco*err-s*cc*s
        a=(sm*cc*err+sco*co*m+s*com*co-sm*co*co-sco*com*err-s*cc*m)/det
        b=(sc*com*err+sm*co*s+s*sco*m-sc*co*m-sm*sco*err-s*com*s)/det
        c=(sc*cc*m+sco*com*s+sm*sco*co-sc*com*co-sco*sco*m-sm*cc*s)/det

        chi[j] = np.sum(weights*((mag-a*np.sin(jdr)-b*np.cos(jdr)-c)**2))/(nobs-3)

        if chi[j] < chi_min: 
            chi_min = chi[j]
            desvstchi = np.sqrt(2*(nobs-3))/(nobs-3)
            period = per[j]
            period_err = (period**2)/2*(np.max(jd)-np.min(jd))

            matrix = np.array([[err,s,co],[s,sc,sco],[co,sco,cc]])
            invmat=inv(matrix)

            phase = np.arctan(b/a)/two_pi
            phase_err=np.sqrt(((1/(two_pi*((a**2)+(b**2))))**2)*(((b**2)*invmat[1,1])+((a**2)*invmat[2,2])))

            amp = 2*a/np.cos(two_pi*phase)
            amp_err = np.sqrt((((2/np.cos(phase*two_pi))**2)*invmat[1,1])+(((2*two_pi*np.sin(phase*two_pi)*a/(np.cos(phase*two_pi)**2))**2)*(phase_err**2)))

            med_mag = c
            medmag_err = np.sqrt(invmat[0,0])

            if amp < 0:
                amp = -amp
                phase = phase + 0.5

            if phase > 1: 
                phase = phase - 1

            if phase <= 0: 
                phase = phase + 1

            jd_maximum = (-phase + 0.75)*period + jd_min
            jd_max_err = np.sqrt((period*phase_err)**2+((-phase+0.75)*period_err)**2)

            if jd_maximum < jd_min: 
                jd_maximum=(-phase+1.75)*period+jd_min
                jd_max_err = np.sqrt((period*phase_err)**2+((-phase+1.75)*period_err)**2)

    jd = jd+jd_min

    return {
        'jd': jd,
        'mag': (mag,mag_err),
        'period':(period, period_err),
        'amplitude': (amp,amp_err),
        'phase':(phase, phase_err),
        'chi':chi, 
        'med_mag': (med_mag, medmag_err),
        'jd_min':jd_min,
        'jd_max':(jd_maximum,jd_max_err)
    }

TEST_FILE=r"C:\Users\cuand\ESAC_Internship\Arecibo Sample\Arecibo_sources_ascii\4596654564903027456.dat"

if __name__=="__main__": 
    data = symmetric_fit(TEST_FILE)
    print(data)