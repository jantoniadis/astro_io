from __future__ import division

import sys
from astropy.table import Table, Column
import numpy as np
import pyfits
import os
import fnmatch
from WDSpec.fitdaspectrum import *
from WDSpec.findvelocity import *
import time
from  astropy.stats.funcs import sigma_clip




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
            'HIERARCH ESO INS SLIT WID':None,
            'HIERARCH ESO TEL GEOELEV':None,
            'HIERARCH ESO TEL GEOLAT':None,
            'HIERARCH ESO TEL GEOLON':None,
            'HIERARCH ESO INS PIXSCALE':None,
            'HIERARCH ESO TEL FWHMLINOBS':None,
            'HIERARCH ESO INS SLIT POSANG':None,
            'HIERARCH ESO INS SLIT LEN':None,
            'HIERARCH ESO INS SHUT EXPTIME':None,
            'HIERARCH ESO TEL IA FWHMLINOBS':None,
            'HIERARCH ESO INS PIXSCALE':None,
            'HIERARCH ESO PRO REC1 PARAM1 VALUE':None}   #dispersion in Ang/pix

    def __init__(self,flux=None,error=None,header=None,mask=False):
        self.flux = flux
        self.error = error
        self.header = header
        self.mask = mask
        self.create_data_table()



    def create_data_table(self):
        self.data = Table()

        #Mask zeros by looking at the first star in table
        if self.mask:
            good = np.where(self.error[0].data[0,:] > 0)
        else:
            good = np.where(self.error[0].data[0,:] > -1e99)

        for i in range(self.flux[0].data.shape[0]):

            #Mask out zeros
            f = Column(self.flux[0].data[i,:][good],name = 'flux_' + self.indexes[i])
            e = Column(self.error[0].data[i,:][good],name = 'error_' + self.indexes[i])
            self.data.add_column(f)
            self.data.add_column(e)

        self.update_keys()
        self.data.meta = self.keys        

        w = Column( np.arange(self.keys['CRVAL1'], self.keys['CRVAL1'] + self.keys['CD1_1']*self.keys['NAXIS1'], self.keys['CD1_1']),name = 'w')
        w = w[good]
        self.data.add_column(w)

    @classmethod
    def readfile(cls,flux_file,error_file,mask=False):
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
        
        fors2_object = cls(flux,error,header,mask)
        return fors2_object



    def update_keys(self):
        for key in self.keys:
            if key in self.header:
                self.keys[key] = self.header[key]




def mid_obs(mjd_start,exptime):
    exp_time_in_days = exptime / 60. / 60. / 24.
    mjd_mid = mjd_start + exp_time_in_days / 2.
    return mjd_mid




class FORSSpectraTable():

    indexes ='abcdefghijklmnopqrstuvxyz'


    def __init__(self,dir='.',n_objects=1):

        self.dir = dir
        self.n_objects = n_objects
        self._fluxes, self._errors = self.get_FORS_spectra()
        self.data = self.create_table()


    def get_FORS_spectra(self):
        fluxlist = []
        errorlist = []
        for root, dirs, files in os.walk(self.dir):
            for name in files:
                if name.endswith(("REDUCED_FLUX_SCI_LSS.fits")):
                    fluxlist.append(os.path.join(root,name))
                elif  name.endswith(("REDUCED_FLUX_ERROR_SCI_LSS.fits")):
                    errorlist.append(os.path.join(root,name))

        return fluxlist, errorlist


    def create_table(self):
        self.mjds = []
        self.positions = []
        self.exptime = []
        self.posangle = []
        self.slit = []
        self.pixscale = []
        self.seeing = []
        self.airmass= []
        self.dispersion = []

        nn = len(self._fluxes)
        for i in range(nn):
            flux_file = self._fluxes[i]
            error_file = self._errors[i]
            fitsobject = FORSReflex.readfile(flux_file,error_file)

            if  i == 0:
                disp_size = fitsobject.data['w'].shape[0]
                wave = np.zeros((nn,disp_size))
                f = np.zeros((nn,disp_size,self.n_objects))
                e = np.zeros((nn,disp_size,self.n_objects))


                self.elev = fitsobject.data.meta['HIERARCH ESO TEL GEOELEV']
                self.lat = fitsobject.data.meta['HIERARCH ESO TEL GEOLAT']
                self.long = fitsobject.data.meta['HIERARCH ESO TEL GEOLON']





            print i

            self.mjds.append(fitsobject.data.meta['MJD-OBS'])
            self.positions.append({fitsobject.data.meta['RA'],fitsobject.data.meta['DEC']})
            self.exptime.append(fitsobject.data.meta['HIERARCH ESO INS SHUT EXPTIME'])
            self.posangle.append(fitsobject.data.meta['HIERARCH ESO INS SLIT POSANG'])
            self.slit.append(fitsobject.data.meta['HIERARCH ESO INS SLIT WID'])
            self.seeing.append(fitsobject.data.meta['HIERARCH ESO TEL IA FWHMLINOBS'])
            self.dispersion.append(fitsobject.data.meta['HIERARCH ESO PRO REC1 PARAM1 VALUE'])
            self.pixscale.append(fitsobject.data.meta['HIERARCH ESO INS PIXSCALE'])
            self.airmass.append(fitsobject.data.meta['HIERARCH ESO OBS AIRM'])

            wave[i,:] = fitsobject.data['w']
            for oo in range(self.n_objects):
                try:
                    f[i,:,oo] = fitsobject.data['flux_' + self.indexes[oo]]
                    e[i,:,oo] = fitsobject.data['error_' + self.indexes[oo]]
                except:
                    pass

        #join flux and e arrays
        r = np.zeros([f.shape[0],wave.shape[-1],self.n_objects,2])
        r[:,:,:,0] = f
        r[:,:,:,1] = e
        r = r.swapaxes(-1,-2)

            #sort by object flux
        spectra_table = np.msort(r.T)[::-1,:,:,:].T

        #mask zeros
        mask=(spectra_table <= 0)
        self.table = np.ma.array(spectra_table,mask=mask)
        self.wave = np.ma.array(wave,mask=mask[:,:,0,0])



def find(pattern,path):
    result = []
    for root,dirs,files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name,pattern):
                result.append(os.path.join(root,name))
    return result


class FORSObservations(FORSSpectraTable):

#  ID MJD Barycentric_MJD EXP_TIME(s) SLIT WIDTH
        def create_obs_table(self):
            self.observations = Table((self.mjds,self.exptime,self.slit,self.dispersion,self.seeing,self.pixscale),
                                names=('mjd','exptime','slit','disp','seeing','pixscale'))
            self.observations['mjd'] = mid_obs(self.observations['mjd'],self.observations['exptime'])



def wd_radial_velocities(tab,teff,
                                logg,
                                vgrid=np.arange(-500.,500.,5.),
                                obj_id=1,
                                pold=3,
                                sigmarange=30,
                                plot=False,
                                clip=30.0):

        i = obj_id
        if plot:
            from matplotlib.pylab import subplots,close
            fig,ax = subplots(2,1)
            ax[0].hold(True)
            ax[1].hold(True)
            plt.ion()
            plt.show()

        try:
            del tab.observations['wdvel']
            del tab.observations['wderr']
        except:
            pass


        wdvel = []
        wderr = []



        for j in range(len(tab.observations['mjd'])):
            print j

            #first load the masked array that contains the flux and error columns
            # and clip it


            clipped = sigma_clip(tab.table[j,:,:,i],sig=clip,iters=3,copy=True,axis=0,cenfunc=np.median)
            s = clipped.data*~clipped.mask

            good = np.where((s[:,0] != 0) & (s[:,1] != 0))

            s = np.array([tab.wave[j,:][good],
                            s[:,0][good],
                            s[:,1][good]])



            s = s.T
            m = DAModel(teff,logg,
                        tab.observations['seeing'][j],
                        tab.observations['slit'][j],
                        tab.observations['pixscale'][j],
                        tab.observations['disp'][j].astype(float))

            if s.size != 0:

                chi2,fit,sol = findvel(s,m,vgrid,pold,sigma=sigmarange,plot=False)
                wdvel.append(chi2.meta['vbest'])
                wderr.append(chi2.meta['verr'])

                if plot:
                    ax[0].cla()
                    ax[1].cla()
                    ax[0].scatter(chi2['v'], chi2['chi2'])
                    ax[0].plot(chi2['v'],polyval(chi2['v'],sol),color='r')
                    ax[1].plot(s[:,0],s[:,1])
                    ax[1].plot(s[:,0],fit,color='g')

                    plt.draw()
                    time.sleep(0.1)
                    plt.pause(0.0001)

            else:
                wdvel.append(0.)
                wderr.append(0.)

        wdvel = Column(wdvel,name='wdvel')
        wderr = Column(wderr,name='wderr')

        tab.observations.add_column(wdvel)
        tab.observations.add_column(wderr)




def calculate_phase(tab,t0,pb):
        circle = Column(((tab.observations['mjd'] - t0)/pb),name='circle')
        phase =Column((circle %1),name='phase')

        tab.observations.add_column(circle)
        tab.observations.add_column(phase)
