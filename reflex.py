from __future__ import division

import sys
from astropy.table import Table, Column
import numpy as np
import pyfits 

class FORSReflex():
    indexes ='abcdefghijklmnopqrstuvxyz'

    keys = {'OBJECT':None,
            'RA':None,
            'DEC':None,
            'EQUINOX':None,
            'MJD-OBS':None,
            'DATE-OBS':None,
            'LST':None,
            'UTC':None,
            'ARCFILE':None,
            'BUNIT':None,
            'NAXIS1':None,
            'NAXIS2':None,
            'CRVAL1':None,
            'CRVAL2':None,
            'CD1_1':None,
            'CD1_2':None,
            'CD2_1':None,
            'CD2_2':None,
            'HIERARCH ESO OBS AIRM':None,
            'HIERARCH ESO TEL GEOELEV':None,
            'HIERARCH ESO TEL GEOLAT':None,
            'HIERARCH ESO TEL GEOLON':None,
            'HIERARCH ESO INS PIXSCALE':None,
            'HIERARCH ESO TEL FWHMLINOBS':None,
            'HIERARCH ESO INS SLIT POSANG':None,
            'HIERARCH ESO INS SLIT LEN':None,
            'HIERARCH ESO INS SHUT EXP TIME':None}

    def __init__(self,flux=None,error=None,header=None):
        self.flux = flux
        self.error = error
        self.header = header
        self.create_data_table()


    def create_data_table(self):
        self.data = Table()
        for i in range(self.flux[0].data.shape[0]):
            f = Column(self.flux[0].data[i,:],name = 'flux_' + self.indexes[i])
            e = Column(self.error[0].data[i,:],name = 'error_' + self.indexes[i])
            self.data.add_column(f)
            self.data.add_column(e)

        self.update_keys()
        self.data.meta = self.keys        

        w = Column( np.arange(self.keys['CRVAL1'], self.keys['CRVAL1'] + self.keys['CD1_1']*self.keys['NAXIS1'], self.keys['CD1_1']),name = 'w')
        self.data.add_column(w)

    @classmethod
    def readfile(cls,flux_file,error_file):
        try:
            flux = pyfits.open(flux_file)
        except (TypeError, NameError):
            print("Flux file does not exist or not in the right format")
        try:
            error = pyfits.open(error_file)
        except (TypeError, NameError):
            print("Flux file does not exist or not in the right format")

        assert flux[0].header['MJD-OBS'] == error[0].header['MJD-OBS']

        header = flux[0].header
        
        fors2_object = cls(flux,error,header)
        return fors2_object


    @staticmethod
    def mid_obs(mjd_start,exptime):
        exp_time_in_days = exptime / 60. / 60. / 24.
        mjd_mid = mjd_start + exp_time_in_days / 2.
        return mjd_mid
        


    def update_keys(self):
        for key in self.keys:
            if key in self.header:
                self.keys[key] = self.header[key]

