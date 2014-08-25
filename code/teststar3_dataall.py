# this is to test a field to get stellar parameters but for the data not for a field 
import pyfits 
import glob 
import pickle
from scipy import interpolate 
from scipy import ndimage 
import numpy as np
from glob import glob 

testfile = "normed_data.pickle"

def test_fits_list(filein): 
  """
  ## inputs
  the file in with the list of fits files want to test - if normed, move on, if not normed, norm it 
  """
  file_in2 = open('normed_data.pickle', 'r') 
  testdata_in,metaall_in,labels, offsets, coeffs, covs, scatters = pickle.load(file_in2)
  file_in2.close()
  nlen = shape(dataall)[1]
  for jj in nlen: 
      testdata[:, jj, 0] = testdata_in[:,jj,0]
      testdata[:, jj, 1] = testdata_in[:,jj,1]
      testdata[:, jj, 2] = testdata_in[:,jj,2] 
  #continuum = get_continuum(testdata) 
  testdata[:, :, 1]# /= continuum
  testdata[:, :, 2]# /= continuum
  bad1 = np.isnan(testdata[:,:,1])
  bad2 = np.isnan(testdata[:,:,2])
  testdata[bad1,1] = 0.
  testdata[bad2,2] = 1000000.0 #np.Inf
  name = 'dataall_self'
  file_in = open(name+'.pickle', 'w')  
  pickle.dump(testdata,  file_in)
  file_in.close()

  return testdata

#def return_test_params(filein,scatters,coeffs,weak_lower=0.935,weak_upper=0.99):
#def return_test_params(filein,scatters,coeffs,weak_lower=0.935,weak_upper=0.99):
#def return_test_params(filein,scatters,coeffs,weak_lower=0.960,weak_upper=.99):
#def return_test_params(filein,scatters,coeffs,weak_lower=0.0,weak_upper=1.99):
def return_test_params(filein,scatters,coeffs,weak_lower=0.97,weak_upper=0.99):
    """
    best log g = weak_lower = 0.95, weak_upper = 0.98
    best teff = weak_lower = 0.95, weak_upper = 0.99
    best_feh = weak_lower = 0.935, weak_upper = 0.98 
    this returns the parameters for a field of data  - and normalises if it is not already normalised 
    """
    #testdata = test_fits_list(testfile) 
    file_in = open(filein, 'r') 
    dataall, metaall, labels, offsets, coeffs, covs, scatters = pickle.load(file_in) 
    file_in.close()
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

testfile = 'coeffs.pickle'
params,icovs_params = return_test_params(testfile,scatters,coeffs)
covs_params = np.linalg.inv(icovs_params) 
params = array(params) 
tme = params[:,0]
gme = params[:,1]
fehme = params[:,2]
all_params = [ params[:,0], params[:,1], params[:,2], (covs_params**0.5)[:,0,0], (covs_params**0.5)[:,1,1], (covs_params**0.5)[:,2,2]] 
file_in = open('datain_self.pickle', 'w')  
file_in = open('coeffs.pickle', 'w')  
pickle.dump((all_params), file_in) 
file_in.close()
# below, this reads in ASPCAP values for comparison for plotting 
testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
file2 = '4332_data_all_more.txt'
file2in = testdir+file2
t,g,feh,feh_err = loadtxt(file2in, usecols = (1,3,5,6), unpack =1) 
