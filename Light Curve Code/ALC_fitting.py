
from symmetric_lc import symmetric_fit
from asymmetric_lc import asymmetric_fit, files_from_directory,read_list_from_file,list_to_file
from lc_plotting import generate_plots

DIR = r"C:\Users\cuand\ESAC Internship\OHIR sample\OHIR_sources_ascii"



TEST_DIR = r"C:\Users\cuand\ESAC Internship\OHIR sample\OHIR_sources_ascii"

TEST_FILE=r"C:\Users\cuand\ESAC_Internship\Arecibo Sample\Arecibo_sources_ascii\4596654564903027456.dat"

DATA_FILE = r"C:\Users\cuand\OneDrive\Documents\LVP Data\new_periods.csv"

#---------------------------------------------------------------------------

def fit_lightcurve(filepath):
    #sym_data = symmetric_fit(filepath)
    #print(sym_data['period'][0])
    #asym_data = asymmetric_fit(filepath, sym_data['period'][0])
    asym_data = asymmetric_fit(filepath)
    #print(asym_data)
    generate_plots(filepath, input_data = asym_data, mission='DR3')



def fit_curves_from_list():
    i=0
    names, periods = read_list_from_file(DATA_FILE)
    file_list = list_to_file(names, DIR)[43:]

    for file, period in zip(file_list, periods):
        print(i)
        asym_data = asymmetric_fit(file, period)
        generate_plots(file, input_data = asym_data, mission='DR3')
        i+=1


def main(): 
    i=0
    files = files_from_directory(TEST_DIR)
    for file in files[:]:
        fit_lightcurve(f"{DIR}\\{file}")
        i+=1


if __name__ == '__main__':
    try:
        main()

    except KeyboardInterrupt:
        pass