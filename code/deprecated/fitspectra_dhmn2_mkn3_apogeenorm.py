# this file is part of the XXX project
# Copyright 2014 Melissa Ness

import scipy
from scipy import interpolate 
from scipy import ndimage 
# import numpy as np

import numpy as np

def weighted_median(values, weights,quantile):
    """
    """
    sindx = np.argsort(values)
    cvalues = 1. * np.cumsum(weights[sindx])
    cvalues = cvalues / cvalues[-1]
    indx = (sindx[cvalues > quantile])[0]
    return values[indx]

def get_continuum(dataall, delta_lambda=50):
    """
    ## inputs:
    dataall:       (Nlambda, Nstar, 3) wavelengths, flux densities, errors
    delta_lambda:  half-width of meadian region in angstroms

    ## output:
    continuum:     (Nlambda, Nstar) continuum level

    ## bugs:
    * for loops!
    """
    Nlambda, Nstar, foo = dataall.shape
    continuum = np.zeros((Nlambda, Nstar))
    assert foo == 3
    for star in range(Nstar):
        print "get_continuum(): working on star" ,star
        for ll, lam in enumerate(dataall[:, 0, 0]):
            assert dataall[ll, star, 0] == lam
            indx = (np.where(abs(dataall[:, star, 0] - lam) < delta_lambda))[0]
            ivar = 1. / (dataall[indx, star, 2] ** 2)
            #continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar)
            continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar,0.80)
    return continuum

def get_data():
  count = 0 
  #T_est,g_est,feh_est = loadtxt("starsin_new.txt", usecols = (5,7,9), unpack =1) 
  #T_est,g_est,feh_est = loadtxt("starsin_M53_N6819.txt", usecols = (4,6,8), unpack =1) 
  T_est,g_est,feh_est = loadtxt("starsin_new_all.txt", usecols = (4,6,8), unpack =1) 
  thismeta = array([T_est, feh_est, g_est])
  thismeta = [T_est, feh_est, g_est]
  #thismeta = zip(T_est, feh_est, g_est)
  #a = open('starsin_new.txt','r')
  #a = open("starsin_M53_N6819.txt", 'r')
  a = open("starsin_new_all.txt", 'r')
  al = a.readlines() 
  bl = []
  for each in al:
    bl.append(each.split()[0]) 
  for jj,each in enumerate(bl):
    count = count+1
    each = each.strip('\n')
    a = pyfits.open(each) 
    b = pyfits.getheader(each) 
    start_wl =  a[1].header['CRVAL1']
    diff_wl = a[1].header['CDELT1']
    print atleast_2d(a[1].data).shape
    if jj == 0:
      nmeta = len(thismeta)
      #nlam = len(a[1].data[0])
      nlam = len(a[1].data)
    val = diff_wl*(nlam) + start_wl 
    wl_full_log = arange(start_wl,val, diff_wl) 
    ydata = (atleast_2d(a[1].data))[0] 
    ydata_err = (atleast_2d(a[2].data))[0] 
    ydata_flag = (atleast_2d(a[3].data))[0] 
    assert len(ydata) == nlam
    wl_full = [10**aval for aval in wl_full_log]
    testdata = scipy.ndimage.gaussian_filter(ydata, 20 ) 
    xdata= wl_full
    xdata =array(xdata)
    ydata =array(ydata)
    ydata_err =array(ydata_err)
    ydata_flag =array(ydata_err)
    #takevalues = logical_and(ydata_err > 10, logical_and(ydata_flag !=1, ydata_flag <= 4500)  
    takevalues = logical_and(ydata  >20, logical_and(ydata_flag >1, ydata_flag <= 3000) )  
    takevalues = logical_and(ydata > 20, logical_and(ydata_flag >1, ydata_flag <= 3000) ) 
    a1 = xdata#[takevalues]
    b1 = ydata#[takevalues]
    b2 = scipy.ndimage.gaussian_filter(b1,1) 
    b3 = scipy.ndimage.gaussian_filter(b1[takevalues],100) 
    a1 = array(a1) 
    b1 = array(b1) 
   # fit1 = pylab.polyfit(array(a1[takevalues]), array(b3), 3) 
   # y1 = polyval(array(fit1), a1) 
    ynew = b2#/y1
    y2new = b2#/y1
    xgrid1 = a1 
    starname2 = each.split('.fits')[0]+'.txt'
    #ivarnew2 = (atleast_2d(a[2].data))[0] /b2
    sigma = (atleast_2d(a[2].data))[0]# /y1
    if jj == 0:
      npix = len(xgrid1) 
      dataall = zeros((npix, len(bl), 3))
     # metaall = ones((len(bl), nmeta))
      metaall = ones((len(bl), nmeta))
      #metaall2 = ones((len(bl), npix, nmeta+1))
    if jj > 0:
      assert xgrid1[0] == dataall[0, 0, 0]
    dataall[:, jj, 0] = xgrid1
    dataall[:, jj, 1] = y2new
    dataall[:, jj, 2] = sigma
    for k in range(0,len(bl)): 
        metaall[k,0] = T_est[k] 
        metaall[k,1] = g_est[k] 
        metaall[k,2] = feh_est[k] 
  predictors = hstack((ones((len(bl), 1)), metaall - mean(metaall, axis=0)[None, :]))
  #continuum = get_continuum(dataall)
  #dataall[:, :, 1] /= continuum
  #dataall[:, :, 2] /= continuum
  return dataall, metaall, predictors, count

def do_one_regression(data, predictors):
  """
  ## inputs
  - data [nobjs, 3] wavelengths, fluxes, invvars
  - meta [nobjs, nmeta] Teff, Feh, etc, etc

  ## outputs
  - chi-squared at best fit
  - coefficients of the fit
  - inverse covariance matrix for fit coefficients
  """
  print "do_one_regression(): working on wavelength", data[0, 0]
  nobj, nmeta = metaall.shape
  assert data.shape == (nobj, 3)
  # least square fit
  Cinv = 1. / (data[:, 2] ** 2) # invvar slice of data
  M = predictors
  MTCinvM = dot(M.T, Cinv[:, None] * M) # craziness b/c Cinv isnt a matrix
  x = data[:, 1] # intensity slice of data
  MTCinvx = dot(M.T, Cinv * x)
  coeff = linalg.solve(MTCinvM, MTCinvx)
  chisq = sum(Cinv * (x - dot(M, coeff)) ** 2)
  chisq_wl = Cinv * (x - dot(M, coeff)) ** 2
  return chisq, coeff, MTCinvM, chisq_wl

def do_regressions(dataall, predictors):
  nlam, nobj, ndata = dataall.shape
  nobj, npred = predictors.shape
  predictorsall = np.zeros((nlam,nobj,npred))
  predictorsall[:, :, :] = predictors[None, :, :]
  return map(do_one_regression, dataall, predictorsall)

def plot_one_fit(dataall, metaall, chisqs, coeffs, invcovs, index):
  return None

dataall, metaall, predictors, count = get_data()
#do_one_regression(dataall[2804], predictors)
blob = do_regressions(dataall, predictors)
chisqs = np.array([b[0] for b in blob])
chisqs_wl = np.array([b[3] for b in blob])
coeffs = np.array([b[1] for b in blob])
invcovs = np.array([b[2] for b in blob])
#predictions = [np.dot(cc, pp) for cc, pp in zip(coeffs, predictors)]
