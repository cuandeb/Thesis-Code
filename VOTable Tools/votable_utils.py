#import the necessary packages
import pandas as pd 
import numpy as np 
from collections import namedtuple
from astropy.io.votable import parse_single_table


#Include a dictionary to convert bands to keys
BAND_DICT = {
    'G':'g_',
    'BP':'bp',
    'RP':'rp'
}


#establish a set of namedtuples representing each of the bands and their Gaia keys
Band =  namedtuple('Band', ['name','source_id','time', 'mag', 'flux','flux_error','phot_flag'])

G_BAND = Band(
    name= 'G Band',
    source_id='source_id',
    time= 'g_transit_time',
    mag='g_transit_mag',
    flux='g_transit_flux',
    flux_error='g_transit_flux_error',
    phot_flag='photometry_flag_sm_reject')
    
BP_BAND = Band(
    name='BP Band',
    source_id='source_id',
    time= 'bp_obs_time',
    mag='bp_mag',
    flux='bp_flux',
    flux_error='bp_flux_error',
    phot_flag='photometry_flag_bp_reject')

RP_BAND = Band(
    name='RP Band',
    source_id='source_id',
    time= 'rp_obs_time',
    mag='rp_mag',
    flux='rp_flux',
    flux_error='rp_flux_error',
    phot_flag='photometry_flag_rp_reject')


#---------------------------------

def check_extension(file,ext):  
    
    """Returns True value if input file is of type ext.
     Otherwise returns False."""    
    
    if file.split('.')[1]==ext:
        return True
    else:
        return False

def read_phot_file(file_path): 
    
    """Extracts tabluar data from a VOTable XML file into
       a pandas dataframe."""
    
    if check_extension(file_path, 'xml'):
        return (parse_single_table(file_path).to_table()).to_pandas()
    else: 
        print('Wrong file format. Please use XML file.') 
    
def find_keys(df,string):
    """Returns the columns of a dataframe that contain
       an input string"""
    return [i for i in df.keys() if string in i]


#---------------------------------

def filter_by_source(table,source):
    
    """Extract data on one specific Gaia source
       using Gaia's source_id label"""
    
    if type(source) != list:
        return table.loc[table['source_id']==source]
    else:
        return table.loc[table['source_id'].isin(source)]
        
    
def filter_by_band(table, band):
    
    """Extract data on one specific Gaia source
       using Gaia's source_id label"""
    
    keys = find_keys(table, BAND_DICT[band])
    return table[['source_id',*keys]]
    
def single_source_band(table, source, band): 
    
    """Extract data for one band for one source"""
    
    _source = filter_by_source(table,source)
    return filter_by_band(_source, band)
    
    
    