# this is to test a field to get stellar parameters 
import pyfits 
import glob 
import pickle
from scipy import interpolate 
from scipy import ndimage 
import numpy as np
from glob import glob 

# this is the file with the list of continuumn normalised APOGEE fits files in a given field 
testfile = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/stars_list_all.txt"
#testfile = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/stars_list_some.txt"

def test_fits_list(filein): 
  """
  ## inputs
  the file in with the list of fits files want to test - if normed, move on, if not normed, norm it 
  """
  name = filein.split('/')[-2]
  if glob(name+'.pickle'):
    file_in2 = open(name+'.pickle', 'r') 
    testdata = pickle.load(file_in2)
    file_in2.close()
    return testdata 
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
     # if jj > 0:
     #     assert xdata[0] == testdata[0, 0, 0]
      testdata[:, jj, 0] = xdata
      testdata[:, jj, 1] = ydata
      testdata[:, jj, 2] = ysigma
  bad3 = testdata[:,:,1] == 0
  testdata[bad3,2] = 1000000.0 #np.Inf
  continuum = get_continuum(testdata) 
  testdata[:, :, 1] /= continuum
  testdata[:, :, 2] /= continuum
  bad1 = np.isnan(testdata[:,:,1])
  bad2 = np.isnan(testdata[:,:,2])
  testdata[bad1,1] = 0.
  testdata[bad2,2] = 1000000.0 #np.Inf
  file_in = open(name+'.pickle', 'w')  
  pickle.dump(testdata,  file_in)
  file_in.close()

  return testdata

#def return_test_params(filein,scatters,coeffs,weak_lower=0.935,weak_upper=0.99):
#def return_test_params(filein,scatters,coeffs,weak_lower=0.935,weak_upper=0.99):
def return_test_params(filein,scatters,coeffs,weak_lower=0.960,weak_upper=.99):
#def return_test_params(filein,scatters,coeffs,weak_lower=0.0,weak_upper=1.99):
    """
    best log g = weak_lower = 0.95, weak_upper = 0.98
    best teff = weak_lower = 0.95, weak_upper = 0.99
    best_feh = weak_lower = 0.935, weak_upper = 0.98 
    this returns the parameters for a field of data  - and normalises if it is not already normalised 
    """
    testdata = test_fits_list(testfile) 
    nstars = shape(testdata)[1]
    Params_all = np.zeros((nstars,  3))
    MCM_rotate_all = np.zeros((nstars, 3,3))
    for jj in range(0,nstars):
      xdata = testdata[:,jj,0]
      ydata = testdata[:,jj,1]
      ysigma = testdata[:,jj,2]
      coeffs_reshape = coeffs[:,-3:]
      ydata_norm = ydata  - coeffs[:,0] 
      coeffs_reshape = coeffs[:,-3:]
      ind1 = logical_and(ydata > weak_lower , ydata < weak_upper)
      Cinv = 1. / (ysigma ** 2 + scatters ** 2)
      MCM_rotate = np.dot(coeffs_reshape[ind1].T, Cinv[:,None][ind1] * coeffs_reshape[ind1])
      MCy_vals = np.dot(coeffs_reshape[ind1].T, Cinv[ind1] * ydata_norm[ind1]) 
     # MCM_rotate = np.dot(coeffs_reshape.T, Cinv[:,None] * coeffs_reshape)
     # MCy_vals = np.dot(coeffs_reshape.T, Cinv * ydata_norm) 
      Params = linalg.solve(MCM_rotate, MCy_vals)
      Params = Params + [mean(metaall[:,0]) , mean(metaall[:,1]), mean(metaall[:,2])] # must be synchronized with fitspectra.py 
      print Params
      Params_all[jj,:] = Params 
      MCM_rotate_all[jj,:,:] = MCM_rotate 
    return Params_all , MCM_rotate_all

params,icovs_params = return_test_params(testfile,scatters,coeffs)
covs_params = np.linalg.inv(icovs_params) 
params = array(params) 
tme = params[:,0]
gme = params[:,1]
fehme = params[:,2]
# below, this reads in ASPCAP values for comparison for plotting 
testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
file2 = '4332_data_all_more.txt'
file2in = testdir+file2
t,g,feh,feh_err = loadtxt(file2in, usecols = (1,3,5,6), unpack =1) 
