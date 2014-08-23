# this file is part of the IR-XXX project
# Copyright 2014 Melissa Ness

import pyfits 
import glob 
import pickle
from scipy import interpolate 
from scipy import ndimage 
import numpy as np

def readinspectra(filein1):
  af = open(filein1,'r')
  al = af.readlines()
  bl = []
  for each in al:
    bl.append(each.split()[0]) 
  for jj, each in enumerate(bl):
    a = pyfits.open(each)
    b = pyfits.getheader(each)
    start_wl =  a[1].header['CRVAL1']
    diff_wl = a[1].header['CDELT1']
    print np.atleast_2d(a[1].data).shape
    if jj == 0:
      nlam = len(a[1].data)
    val = diff_wl*(nlam) + start_wl 
    wl_full_log = np.arange(start_wl,val, diff_wl) 
    ydata = (np.atleast_2d(a[1].data))[0] 
    ydata_err = (np.atleast_2d(a[2].data))[0] 
    ydata_flag = (np.atleast_2d(a[3].data))[0] 
    assert len(ydata) == nlam
    wl_full = [10**aval for aval in wl_full_log]
    xdata= wl_full
    xdata =np.array(xdata)
    ydata =np.array(ydata)
    if jj == 0:
      npix = len(xdata) 
      test_y = np.zeros((npix, len(bl), 3))
    if jj > 0:
      assert xdata[0] == test_y[0, 0,0]
    test_y[:, jj, 0] = xdata
    test_y[:, jj, 1] = ydata
    test_y[:, jj, 2] = ydata_err
  return test_y

# this is to return the parameters for all of the existing stars in the calibration   - this is not cross validation 
def getpar_one(coeffs,scatters,metaall, starnumber): 
  """
  return stellar parameters for a single or multiple input spectrum : this is not for cross validation this is just to see how it compares 
  inputs = all of the outputs from the fitspectra, the data, the coefficients, the scatters - going to be used to get params given a spectra 
  testspectrum = the spectra we want to determine the parameters of 
  # inputs go here
  """
  xdata = dataall[:,300,0] #testspectrum[:,0,0]
  ydata = dataall[:,300,1]  #testspectrum[:,0,1]
  ysigma = dataall[:,300,2] #testspectrum[:,0,2]
  #bad3 = dataall[:,300,1] == 0. 
  #dataall[:,300,2][bad3] == 100000. 
  #ysigma[bad3] = 100000. #testspectrum[:,0,2]
  coeffs_reshape = coeffs[:,-3:]
  ydata_norm = ydata  - coeffs[:,0] 
  coeffs_reshape = coeffs[:,-3:]
  Cinv = 1. / (ysigma ** 2 + scatters ** 2)
  #ind1 = logical_and(xdata > 16055, xdata < 16061) 
  xdata = array(xdata) 
  ind1 = logical_and(xdata > 14400, xdata < 19780) 
  #indb = logical_or(logical_and(xdata > 15894, xdata < 16389) , inda) 
  #ind1 = logical_or(indb, indc)
  #ind1 = logical_and(scatters < 0.08 , scatters > 0.0001) 
  #ind1 =  scatters > 0.00010 
  #MCM_rotate = np.dot(coeffs_reshape[ind1].T, Cinv[:,None][ind1] * coeffs_reshape[ind1])
  #MCy_vals = np.dot(coeffs_reshape[ind1].T, Cinv[ind1] * ydata_norm[ind1]) 
  MCM_rotate = np.dot(coeffs_reshape.T, Cinv[:,None] * coeffs_reshape)
  MCy_vals = np.dot(coeffs_reshape.T, Cinv * ydata_norm) 
  Params = linalg.solve(MCM_rotate, MCy_vals)
  Params = Params + [mean(metaall[:,0]) , mean(metaall[:,1]), mean(metaall[:,2])] # must be synchronized with fitspectra.py 
  print MCM_rotate, MCy_vals, Params 
  return Params,Cinv

params_one,Cinv_Params = getpar_one(coeffs, scatters, metaall,"testspectrum.txt") 
  

#
#
#
#dataall, metaall, predictors, count = get_data()
#blob = do_regressions(dataall, predictors)
#coeffs = np.array([b[0] for b in blob])
#invcovs = np.array([b[1] for b in blob])
#covs = np.array([np.linalg.inv(b[1]) for b in blob])
#chis = np.array([b[2] for b in blob])
#chisqs = np.array([np.dot(b[2],b[2]) - b[3] for b in blob]) # holy crap be careful
#scatters = np.array([b[4] for b in blob])
