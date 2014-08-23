# this file is part of the IR analysis project
# Copyright 2014 Melissa Ness
# use this one, unless you copy it Aug 18 2014 

import pyfits 
import glob 
import pickle
from scipy import interpolate 
from scipy import ndimage 
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
            ivar = array(ivar) 
            bad = isinf(ivar)
            ivar[bad] =1.
            continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar, 0.90)
    return continuum

def get_data():
  #if glob.glob('normed_data.pickle'): # MOVE THIS WAY UP AND PICKLE metaall TOO
  #  file_in2 = open('normed_data.pickle', 'r') 
  #  dataall,metaall,predictors,count = pickle.load(file_in2)
  #  file_in2.close()
  #  return dataall, metaall, predictors, count
  count = 0 
  T_est,g_est,feh_est = np.loadtxt("starsin_test.txt", usecols = (4,6,8), unpack =1) 
  T_est,g_est,feh_est = np.loadtxt("starsin_test3.txt", usecols = (4,6,8), unpack =1) 
  T_est,g_est,feh_est = np.loadtxt("starsin_test2.txt", usecols = (4,6,8), unpack =1) 
  T_est,g_est,feh_est = np.loadtxt("starsin_test2.txt", usecols = (4,6,8), unpack =1) 
  T_est,g_est,feh_est = np.loadtxt("starsin_new_all_ordered.txt", usecols = (4,6,8), unpack =1) 
  #T_est = loadtxt("Temperature_Alonso.txt", usecols = (0,), unpack =1) 
  thismeta = np.array([T_est, feh_est, g_est])
  thismeta = [T_est, feh_est, g_est]
  #a = open("starsin_new_all.txt", 'r')
  #a = open("starsin_M53_N6819.txt", 'r')
  a = open("starsin_test.txt", 'r')
  a = open("starsin_test3.txt", 'r')
  a = open("starsin_test2.txt", 'r')
  a = open("starsin_test2.txt", 'r')
  a = open("starsin_new_all_ordered.txt", 'r')
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
    print np.atleast_2d(a[1].data).shape
    if jj == 0:
      nmeta = len(thismeta)
      nlam = len(a[1].data)
    val = diff_wl*(nlam) + start_wl 
    wl_full_log = np.arange(start_wl,val, diff_wl) 
    ydata = (np.atleast_2d(a[1].data))[0] 
    ydata_err = (np.atleast_2d(a[2].data))[0] 
    ydata_flag = (np.atleast_2d(a[3].data))[0] 
    assert len(ydata) == nlam
    wl_full = [10**aval for aval in wl_full_log]
    testdata = scipy.ndimage.gaussian_filter(ydata, 20 ) 
    xdata= wl_full
    xdata =np.array(xdata)
    ydata =np.array(ydata)
    ydata_err =np.array(ydata_err)
    ydata_flag =np.array(ydata_err)
    a1 = xdata
    b1 = ydata
    b2 = scipy.ndimage.gaussian_filter(b1,1) 
    a1 = np.array(a1) 
    b1 = np.array(b1) 
    ynew = b2
    y2new = b2
    xgrid1 = a1 
    starname2 = each.split('.fits')[0]+'.txt'
    sigma = (np.atleast_2d(a[2].data))[0]# /y1
    if jj == 0:
      npix = len(xgrid1) 
      dataall = np.zeros((npix, len(bl), 3))
      metaall = np.ones((len(bl), nmeta))
    if jj > 0:
      assert xgrid1[0] == dataall[0, 0, 0]
    dataall[:, jj, 0] = xgrid1
    dataall[:, jj, 1] = y2new
    dataall[:, jj, 2] = sigma
    bad = np.logical_or(np.isnan(sigma), np.isnan(y2new)) 
    print "get_data(): fixing %d bad values" % np.sum(bad)
    dataall[bad, jj, 1] = 0.
    dataall[bad, jj, 2] = np.Inf
    bad = np.isinf(dataall[:, jj, 1])

    for k in range(0,len(bl)): 
      metaall[k,0] = T_est[k] 
      metaall[k,1] = g_est[k] 
      metaall[k,2] = feh_est[k] 
  predictors = np.hstack((np.ones((len(bl), 1)), metaall - np.mean(metaall, axis=0)[None, :]))
  continuum = get_continuum(dataall)
  dataall[:, :, 1] /= continuum
  dataall[:, :, 2] /= continuum
  bad1 = np.isnan(dataall[:,:,1])
  # commented out just for now the badness
  dataall[bad1,1] = 0.
  ##bad2 = np.isnan(dataall[:,:,2])
  ##dataall[bad2,2] = 100000.0#np.Inf
  ##bad3 = dataall[:,:,1] == 0.
  ##dataall[bad3,2] = 100000.0#np.Inf# trying this - may be wrong but trying 
  #bad3 = dataall[:,:,2] == 0.
  #dataall[bad3,2] = 100000.0#np.Inf
  #bad = np.logical_or(np.isnan(dataall[:, jj, 2]), np.isnan(dataall[:, jj, 1])) 
  print "get_data(): fixing %d bad values" % np.sum(bad)

  #file_in = open('normed_data.pickle', 'w')  
  #pickle.dump((dataall,metaall, predictors, count),  file_in)
  #file_in.close()

  return dataall, metaall, predictors, count

#def do_one_regression_at_fixed_scatter(data, predictors, scatter=0):
def do_one_regression_at_fixed_scatter(data, predictors, scatter):
  """
  ## inputs
  - data [nobjs, 3] wavelengths, fluxes, invvars
  - meta [nobjs, nmeta] Teff, Feh, etc, etc
  - scatter
  ## outputs
  - chi-squared at best fit
  - coefficients of the fit
  - inverse covariance matrix for fit coefficients
  """
  nobj, nmeta = metaall.shape
  assert data.shape == (nobj, 3)
  # least square fit
  Cinv = 1. / (data[:, 2] ** 2 +scatter**2)  # invvar slice of data
  M = predictors
  MTCinvM = np.dot(M.T, Cinv[:, None] * M) # craziness b/c Cinv isnt a matrix
  x = data[:, 1] # intensity slice of data
  MTCinvx = np.dot(M.T, Cinv * x)
  coeff = np.linalg.solve(MTCinvM, MTCinvx)
  # below is all new 
  #gp = GaussianProcess(corr='cubic', theta0=1e-2, thetaL=1e-4, thetaU=1e-1,
  #                   random_start=100)
  # I removed this below - thought we should be checking for nan not assert isfinite
  #assert np.all(np.isfinite(coeff)) 
  #chisq = sum(Cinv * (x - dot(M, coeff)) ** 2)
  #chiq_wl =Cinv * (x - dot(M, coeff)) ** 2
  chi = np.sqrt(Cinv) * (x - np.dot(M, coeff)) 
  logdet_Cinv = np.sum(np.log(Cinv)) 
  # resid2 = sum( (x - dot(M, coeff)) ** 2)
  bad1 = isnan(data) 
  bad2 = isnan(coeff) 
  bad3 = isnan(Cinv) 
  return (coeff, MTCinvM, chi, logdet_Cinv )

def do_one_regression(data, metadata):
    """
    blah blah blah.
    # inputs:
    """
    ln_s_values = np.arange(np.log(0.0001), 0., 0.5)
    chis_eval = np.zeros_like(ln_s_values)
    for ii, ln_s in enumerate(ln_s_values):
        foo, bar, chi, logdet_Cinv = do_one_regression_at_fixed_scatter(data, metadata, scatter = np.exp(ln_s))
        chis_eval[ii] = np.sum(chi * chi) - logdet_Cinv
    if np.any(np.isnan(chis_eval)):
        s_best = np.exp(ln_s_values[-1])
        return do_one_regression_at_fixed_scatter(data, metadata, scatter = s_best) + (s_best, )
    lowest = argmin(chis_eval)
    if lowest == 0 or lowest == len(ln_s_values) + 1:
        s_best = np.exp(ln_s_values[lowest])
        return do_one_regression_at_fixed_scatter(data, metadata, scatter = s_best) + (s_best, )
    ln_s_values_short = ln_s_values[np.array([lowest-1, lowest, lowest+1])]
    chis_eval_short = chis_eval[np.array([lowest-1, lowest, lowest+1])]
    z = np.polyfit(ln_s_values_short, chis_eval_short, 2)
    f = np.poly1d(z)
    fit_pder = np.polyder(z)
    fit_pder2 = pylab.polyder(fit_pder)
    s_best = np.exp(np.roots(fit_pder)[0])
    return do_one_regression_at_fixed_scatter(data, metadata, scatter = s_best) + (s_best, )

def do_regressions(dataall, predictors):
  nlam, nobj, ndata = dataall.shape
  nobj, npred = predictors.shape
  predictorsall = np.zeros((nlam,nobj,npred))
  predictorsall[:, :, :] = predictors[None, :, :]
  return map(do_one_regression, dataall, predictorsall)

def plot_one_fit(dataall, metaall, chisqs, coeffs, invcovs, index):
  return None

#covars = np.array([linalg.inv(cinv) for cinv in invcovs]) 

dataall, metaall, predictors, count = get_data()
blob = do_regressions(dataall, predictors)
coeffs = np.array([b[0] for b in blob])
invcovs = np.array([b[1] for b in blob])
covs = np.array([np.linalg.inv(b[1]) for b in blob])
chis = np.array([b[2] for b in blob])
chisqs = np.array([np.dot(b[2],b[2]) - b[3] for b in blob]) # holy crap be careful
scatters = np.array([b[4] for b in blob])
bad1 = np.isnan(dataall[:,:,2])
bad2 = np.isnan(dataall[:,:,1])
bad3 = np.isnan(coeffs) 
dataall[bad1,2] = np.Inf 
dataall[bad2,1] = 0. 
coeffs[bad3] = 0. 

#predictions = [np.dot(cc, pp) for cc, pp in zip(coeffs, predictors)]
#savetxt("data_test.txt", zip(dataall[:,0,0], coeffs[:,0], coeffs[:,1], coeffs[:,2], chisqs), fmt = "%s")
#savetxt("data_test.txt", zip(coeffs[:,0], coeffs[:,1], coeffs[:,2], coeffs[:,3], chis, chisqs), fmt = "%s")


