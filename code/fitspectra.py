"""
This file is part of the IR analysis project.
Copyright 2014 Melissa Ness.
# http://iopscience.iop.org/1538-3881/146/5/133/suppdata/aj485195t4_mrt.txt for calibration stars 
# http://data.sdss3.org/irSpectrumDetail?locid=4330&commiss=0&apogeeid=2M17411636-2903150&show_aspcap=True object explorer 
# http://data.sdss3.org/basicIRSpectra/searchStarA
# http://data.sdss3.org/sas/dr10/apogee/spectro/redux/r3/s3/a3/ for the data files 
# to-do
- need to add a test that the wavelength range is the same - and if it isn't interpolate to the same range 
- format PEP8-ish (four-space tabs, for example)
- take logg_cut as an input
- extend to perform quadratic fitting
"""

from astropy.io import fits as pyfits 
import scipy 
import glob 
import pickle
import pylab 
from scipy import interpolate 
from scipy import ndimage 
import numpy as np
LARGE = 1e2 # sigma value to use for bad continuum-normalized data; MAGIC

def weighted_median(values, weights, quantile):
    """
    """
    sindx = np.argsort(values)
    cvalues = 1. * np.cumsum(weights[sindx])
    cvalues = cvalues / cvalues[-1]
    foo = sindx[cvalues > quantile]
    if len(foo) == 0:
        return values[0]
    indx = foo[0]
    return values[indx]

def continuum_normalize(dataall, delta_lambda=50):
    """
    ## inputs:
    dataall:       (Nlambda, Nstar, 3) wavelengths, flux densities, errors
    delta_lambda:  half-width of meadian region in angstroms

    ## output:
    continuum:     (Nlambda, Nstar) continuum level

    ## comments:
    * does a lot of stuff *other* than continuum normalization

    ## bugs:
    * for loops!
    """
    Nlambda, Nstar, foo = dataall.shape
    continuum = np.zeros((Nlambda, Nstar))
    # sanitize inputs
    for jj in range(Nstar):
    #    #BROKEN
        bad_a = np.logical_or(np.isnan(dataall[:, jj, 1]) ,np.isinf(dataall[:,jj, 1]))
        bad_b = np.logical_or(dataall[:, jj, 2] <= 0. , np.isnan(dataall[:, jj, 2]))
        bad = np.logical_or( np.logical_or(bad_a, bad_b) , np.isinf(dataall[:, jj, 2]))
        dataall[bad, jj, 1] = 0.
        dataall[bad, jj, 2] = np.Inf #LARGE#np.Inf #100. #np.Inf
        continuum = np.zeros((Nlambda, Nstar))
    assert foo == 3
    for star in range(Nstar):
        print "get_continuum(): working on star" ,star
        for ll, lam in enumerate(dataall[:, 0, 0]):
            if dataall[ll, star, 0] != lam:
                print dataall[ll,star,0], lam , dataall[ll,0,0] 
                print ll, star 
                print ll+1, star+1, dataall[ll+1, star+1, 0], dataall[ll+1,0,0] 
                print ll+2, star+2, dataall[ll+2, star+2, 0], dataall[ll+2,0,0] 
                assert False
            indx = (np.where(abs(dataall[:, star, 0] - lam) < delta_lambda))[0]
            ivar = 1. / (dataall[indx, star, 2] ** 2)
            ivar = np.array(ivar)
            continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar, 0.90)
    # sanitize outputs
    #dataall[:, :, 1] /= continuum
    #dataall[:, :, 2] /= continuum
    #bad1 = np.isnan(dataall[:,:,1])
    #dataall[bad1,1] = 0.
    #bad2 = np.isnan(dataall[:,:,2])
    #dataall[bad2,2] = 1000000. #note if infinity falls over later 
    for jj in range(Nstar):
        bad = np.where(continuum[:,jj] <= 0) 
        continuum[bad,jj] = 1.
        dataall[:, jj, 1] /= continuum[:,jj]
        dataall[:, jj, 2] /= continuum[:,jj]
        dataall[bad,jj, 1] = 1. 
        dataall[bad,jj, 2] = LARGE 
        bad = np.where(dataall[:, jj, 2] > LARGE) 
        dataall[bad,jj, 1] = 1. 
        dataall[bad,jj, 2] = LARGE 

    #    #BROKEN
    #    bad = np.where(continuum[:, jj] <= 0) 
    #    continuum[bad] = 1.
    #    dataall[:, :, 1] /= continuum[:,jj]
    #    dataall[:, :, 2] /= continuum[:,jj]
    #    dataall[bad, 1] = 1. 
    #    dataall[bad, 2] = LARGE 
    #    bad = np.where(dataall[:, :, 2] > LARGE) 
    #    dataall[bad, 1] = 1. 
    #    dataall[bad, 2] = LARGE 
    return dataall 

def get_normalized_test_data(testfile): 
  """
  ## inputs
  the file in with the list of fits files want to test - if normed, move on, if not normed, norm it 
  """
  name = testfile.split('/')[-2]
  testdir = testfile.split('stars')[0]
  #if glob.glob(name+'.pickle'):
  #    file_in2 = open(name+'.pickle', 'r') 
  #    testdata = pickle.load(file_in2)
  #    file_in2.close()
  #    return testdata 

  a = open(testfile, 'r')
  al2 = a.readlines()
  bl2 = []
  for each in al2:
    bl2.append(testdir+each.strip())
  for jj,each in enumerate(bl2):
    a = pyfits.open(each) 
    ydata = a[1].data
    ysigma = a[2].data
    start_wl =  a[1].header['CRVAL1']
    diff_wl = a[1].header['CDELT1']
    if jj == 0:
        nlam = len(a[1].data)
        testdata = np.zeros((nlam, len(bl2), 3))
    val = diff_wl*(nlam) + start_wl 
    wl_full_log = np.arange(start_wl,val, diff_wl) 
    wl_full = [10**aval for aval in wl_full_log]
    xdata = wl_full
    testdata[:, jj, 0] = xdata
    testdata[:, jj, 1] = ydata
    testdata[:, jj, 2] = ysigma
  testdata = continuum_normalize(testdata) # testdata
  file_in = open(name+'.pickle', 'w')  
  pickle.dump(testdata,  file_in)
  file_in.close()
  return testdata 

def get_normalized_training_data():
  if glob.glob('normed_data.pickle'): 
        file_in2 = open('normed_data.pickle', 'r') 
        dataall, metaall, labels = pickle.load(file_in2)
        file_in2.close()
        return dataall, metaall, labels
  fn = "starsin_test2.txt"
  fn = "starsin_test.txt"
  fn = "starsin_new_all_ordered.txt"
  T_est,g_est,feh_est = np.loadtxt(fn, usecols = (4,6,8), unpack =1) 
  labels = ["teff", "logg", "feh"]
  a = open(fn, 'r') 
  al = a.readlines() 
  bl = []
  for each in al:
    bl.append(each.split()[0]) 
  for jj,each in enumerate(bl):
    each = each.strip('\n')
    a = pyfits.open(each) 
    b = pyfits.getheader(each) 
    start_wl =  a[1].header['CRVAL1']
    diff_wl = a[1].header['CDELT1']
    print np.atleast_2d(a[1].data).shape
    if jj == 0:
      nmeta = len(labels)
      nlam = len(a[1].data)
    val = diff_wl*(nlam) + start_wl 
    wl_full_log = np.arange(start_wl,val, diff_wl) 
    ydata = (np.atleast_2d(a[1].data))[0] 
    ydata_err = (np.atleast_2d(a[2].data))[0] 
    ydata_flag = (np.atleast_2d(a[3].data))[0] 
    assert len(ydata) == nlam
    wl_full = [10**aval for aval in wl_full_log]
    xdata= np.array(wl_full)
    ydata = np.array(ydata)
    ydata_err = np.array(ydata_err)
    starname2 = each.split('.fits')[0]+'.txt'
    sigma = (np.atleast_2d(a[2].data))[0]# /y1
    if jj == 0:
      npix = len(xdata) 
      dataall = np.zeros((npix, len(bl), 3))
      metaall = np.ones((len(bl), nmeta))
    if jj > 0:
      assert xdata[0] == dataall[0, 0, 0]

    dataall[:, jj, 0] = xdata
    dataall[:, jj, 1] = ydata
    dataall[:, jj, 2] = sigma

    for k in range(0,len(bl)): 
        # must be synchronised with labels 
      metaall[k,0] = T_est[k] 
      metaall[k,1] = g_est[k] 
      metaall[k,2] = feh_est[k] 
  dataall = continuum_normalize(dataall) #dataall
  offsets = np.mean(metaall, axis=0)
  predictors = np.hstack((np.ones((len(metaall), 1)), metaall - offsets[None, :]))

  file_in = open('normed_data.pickle', 'w')  
  pickle.dump((dataall, metaall, labels),  file_in)
  file_in.close()
  return dataall, metaall, labels 

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
    # least square fit
    Cinv = 1. / (data[:, 2] ** 2 +scatter**2)  # invvar slice of data
    M = predictors
    MTCinvM = np.dot(M.T, Cinv[:, None] * M) # craziness b/c Cinv isnt a matrix
    x = data[:, 1] # intensity slice of data
    MTCinvx = np.dot(M.T, Cinv * x)
    coeff = np.linalg.solve(MTCinvM, MTCinvx)
    assert np.all(np.isfinite(coeff)) 
    chi = np.sqrt(Cinv) * (x - np.dot(M, coeff)) 
    logdet_Cinv = np.sum(np.log(Cinv)) 
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
    lowest = np.argmin(chis_eval)
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
    """
    """
    nlam, nobj, ndata = dataall.shape
    nobj, npred = predictors.shape
    predictorsall = np.zeros((nlam,nobj,npred))
    predictorsall[:, :, :] = predictors[None, :, :]
    return map(do_one_regression, dataall, predictorsall)

def read_and_train(dataall, metaall, order, fn, logg_cut=100., teff_cut=0., leave_out=None):
    """
    - `leave out` must be in the correct form to be an input to `np.delete`
    """
    #good = np.logical_and((metaall[:, 1] < logg_cut), (metaall[:,0] > teff_cut) ) 
    dataall = dataall#[:, good]
    metaall = metaall#[good]

    if leave_out is not None:
        dataall = np.delete(dataall, leave_out, axis = 1) 
        metaall = np.delete(metaall, leave_out, axis = 0) 

    assert order == 1 # because we suck
    if order == 1:
        offsets = np.mean(metaall, axis=0)
        #predictors = np.hstack((np.ones((len(metaall), 1)), metaall - offsets[None, :]))
        predictors = np.hstack((np.ones((len(metaall), 1)), metaall - offsets)) 

    blob = do_regressions(dataall, predictors)
    coeffs = np.array([b[0] for b in blob])
    #invcovs = np.array([b[1] for b in blob])
    covs = np.array([np.linalg.inv(b[1]) for b in blob])
    #chis = np.array([b[2] for b in blob])
    #chisqs = np.array([np.dot(b[2],b[2]) - b[3] for b in blob]) # holy crap be careful
    scatters = np.array([b[4] for b in blob])

    fd = open(fn, "w")
    pickle.dump((dataall, metaall, labels, offsets, coeffs, covs, scatters), fd)
    fd.close()


def infer_tags(fn_pickle,testdata, fout_pickle, weak_lower,weak_upper):
#def infer_tags(fn_pickle,testdata, fout_pickle, weak_lower=0.935,weak_upper=0.98):
    """
    best log g = weak_lower = 0.95, weak_upper = 0.98
    best teff = weak_lower = 0.95, weak_upper = 0.99
    best_feh = weak_lower = 0.935, weak_upper = 0.98 
    this returns the parameters for a field of data  - and normalises if it is not already normalised 
    this is slow because it reads a pickle file 
    """
    file_in = open(fn_pickle, 'r') 
    dataall, metaall, labels, offsets, coeffs, covs, scatters = pickle.load(file_in)
    file_in.close()
    nstars = (testdata.shape)[1]
    ntags = len(labels)
    Params_all = np.zeros((nstars, ntags))
    MCM_rotate_all = np.zeros((nstars, ntags, ntags))
    for jj in range(0,nstars):
      if np.any(testdata[:,jj,0] != dataall[:, 0, 0]):
          print testdata[range(5),jj,0], dataall[range(5),0,0]
          assert False
      xdata = testdata[:,jj,0]
      ydata = testdata[:,jj,1]
      ysigma = testdata[:,jj,2]
      ydata_norm = ydata  - coeffs[:,0] # subtract the mean 
      coeffs_slice = coeffs[:,-3:]
      ind1 = np.logical_and(ydata > weak_lower , ydata < weak_upper)
      Cinv = 1. / (ysigma ** 2 + scatters ** 2)
      MCM_rotate = np.dot(coeffs_slice[ind1].T, Cinv[:,None][ind1] * coeffs_slice[ind1])
      MCy_vals = np.dot(coeffs_slice[ind1].T, Cinv[ind1] * ydata_norm[ind1]) 
      Params = np.linalg.solve(MCM_rotate, MCy_vals)
      Params = Params + offsets 
      print Params
      Params_all[jj,:] = Params 
      MCM_rotate_all[jj,:,:] = MCM_rotate 
    file_in = open(fout_pickle, 'w')  
    pickle.dump((Params_all, MCM_rotate_all),  file_in)
    file_in.close()
    return Params_all , MCM_rotate_all

if __name__ == "__main__":
    dataall, metaall, labels = get_normalized_training_data()
    fpickle = "coeffs.pickle" 
    if not glob.glob(fpickle):
        #read_and_train(dataall, metaall, 1,  fpickle, logg_cut= 4.)
        read_and_train(dataall, metaall, 1,  fpickle, logg_cut= 40.,teff_cut = 0.)# - did this before and makes pleaides out 
    testfile = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/stars_list_all.txt"
    self_flag = 0
    if self_flag != 1:
      field = "4332_"
      testdataall = get_normalized_test_data(testfile) # if flag is one, do on self 
      testmetaall, inv_covars = infer_tags("coeffs.pickle", testdataall, field+"tags.pickle",-10.94,10.99) 
    if self_flag == 1:
      field = "self_"
      file_in = open('normed_data.pickle', 'r') 
      testdataall, metaall, labels = pickle.load(file_in)
      file_in.close() 
      testmetaall, inv_covars = infer_tags("coeffs.pickle", testdataall, field+"tags.pickle",-10.94,10.99) 
    
    #def labels_on_apogee_data(fin,fsave):
        #fsave = "labels.pickle"
        #dataall, metaall, labels = get_normalized_training_data()


if False: 
    testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
    file2 = '4332_data_all_more.txt'
    file2in = testdir+file2
    t,g,feh,feh_err = loadtxt(file2in, usecols = (1,3,5,6), unpack =1) 
