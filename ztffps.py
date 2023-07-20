from matplotlib.pyplot import cycler
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.cm
from glob import glob
from astropy.time import Time
from astropy.table import Table
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii, fits
import pandas as pd
from astropy.coordinates import SkyCoord



def galactic_latitude(ra, dec):
    # l_ref = 33.012 # deg
    ra_ref = 282.25 # deg
    g = 62.6 # deg
    b =  np.arcsin(np.sin(np.deg2rad(dec)) * np.cos(np.deg2rad(g)) - \
                   np.cos(np.deg2rad(dec)) * np.sin(np.deg2rad(g)) * np.sin(np.deg2rad(ra) - np.deg2rad(ra_ref)))
    return np.rad2deg(b)

def get_cycle(cmap, N=None, use_index="auto"):
    if isinstance(cmap, str):
        if use_index == "auto":
            if cmap in ['Pastel1', 'Pastel2', 'Paired', 'Accent',
                        'Dark2', 'Set1', 'Set2', 'Set3',
                        'tab10', 'tab20', 'tab20b', 'tab20c']:
                use_index=True
            else:
                use_index=False
        cmap = matplotlib.cm.get_cmap(cmap)
    if not N:
        N = cmap.N
    if use_index=="auto":
        if cmap.N > 100:
            use_index=False
        elif isinstance(cmap, LinearSegmentedColormap):
            use_index=False
        elif isinstance(cmap, ListedColormap):
            use_index=True
    if use_index:
        ind = np.arange(int(N)) % cmap.N
        return cycler("color",cmap(ind))
    else:
        colors = cmap(np.linspace(0,1,N))
        return cycler("color",colors)
    
def get_lightcurve(data_file):
    ztf_fp = pd.read_csv(data_file, comment='#', index_col='index,', delim_whitespace=True)
    ztf_fp.columns = [x.split(',')[0] for x in ztf_fp.columns]
    ztf_fp.index.name = 'index'
    ztf_fp['mjd'] = Time(ztf_fp['jd'],format='jd').mjd
    
    # check for nonzero procstatus
    wfp_bad = (ztf_fp.procstatus.apply(lambda x: str(x)) != '0')
    np.sum(wfp_bad)
    # print(ztf_fp.loc[wfp_bad,['field','filter','ccdid','procstatus']].groupby(['field','filter','procstatus']).agg(len))
    
    # print(ztf_fp.loc[~wfp_bad,['field','filter','ccdid','procstatus']].groupby(['field','filter','procstatus']).agg(len))
    

    # TODO: need to check that nearest ref source is coincident and stellar
    ztf_fp.loc[:,'nearestrefflux'] = 10**(0.4 * (ztf_fp['zpdiff'] - ztf_fp['nearestrefmag'])) 
    ztf_fp.loc[:,'nearestreffluxunc'] = ztf_fp['nearestrefmagunc'] * ztf_fp['nearestrefflux']/ 1.0857
    ztf_fp.loc[:,'dc_flux'] = ztf_fp['forcediffimflux'] + ztf_fp['nearestrefflux']
    ztf_fp.loc[:,'dc_flux_err'] = np.sqrt(ztf_fp['forcediffimfluxunc']**2. - ztf_fp['nearestreffluxunc']**2.)
    # fix any imaginary ones
    ztf_fp.loc[ztf_fp.loc[:,'dc_flux_err'].isna(),'dc_flux_err'] = np.sqrt(ztf_fp['forcediffimfluxunc']**2. + ztf_fp['nearestreffluxunc']**2.)
    ztf_fp.loc[:,'dc_snr'] = ztf_fp['dc_flux']/ztf_fp['dc_flux_err']
    SNT = 3 # 3-sigma detection
    SNU = 5 # 5-sigma upper limit

    ztf_fp.loc[:,'dc_mag'] = np.NaN
    ztf_fp.loc[:,'dc_mag_err'] = np.NaN
    ztf_fp.loc[:,'dc_maglim'] = np.NaN

    wzdet = ztf_fp['dc_snr'] > SNT
    ztf_fp.loc[wzdet,'dc_mag'] = ztf_fp.loc[wzdet,'zpdiff'] - 2.5 * np.log10(ztf_fp.loc[wzdet,'dc_flux'])
    ztf_fp.loc[wzdet,'dc_mag_err'] = 1.0857 / ztf_fp['dc_snr']

    ztf_fp.loc[~wzdet,'dc_maglim'] = ztf_fp.loc[~wzdet,'zpdiff'] - 2.5 * np.log10(SNU*ztf_fp.loc[~wzdet,'dc_flux_err'])


    # if we don't have reference sources the calculation is slightly different
    wnoref = ztf_fp['nearestrefmag'].isna()
    wnorefdet = wnoref & ((ztf_fp.loc[wnoref,'forcediffimflux'] / ztf_fp.loc[wnoref,'forcediffimfluxunc']) > SNT)
    wnorefnondet = wnoref & ~((ztf_fp.loc[wnoref,'forcediffimflux'] / ztf_fp.loc[wnoref,'forcediffimfluxunc']) > SNT)

    ztf_fp.loc[wnorefdet,'dc_mag'] = ztf_fp.loc[wnorefdet,'zpdiff'] - 2.5 * np.log10(ztf_fp.loc[wnorefdet,'forcediffimflux'])
    ztf_fp.loc[wnorefdet,'dc_mag_err'] = 1.0857 / (ztf_fp['forcediffimflux'] / ztf_fp.loc[wnorefdet,'forcediffimfluxunc'])
    ztf_fp.loc[wnorefnondet,'dc_maglim'] = ztf_fp.loc[wnorefnondet,'zpdiff'] - 2.5 * np.log10(SNU*ztf_fp.loc[wnorefnondet,'forcediffimfluxunc'])

    return ztf_fp