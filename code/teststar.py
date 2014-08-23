testspectrum = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/aspcapStar-v304-2M18024372-2921285.fits" 
testfile = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/stars_list_all.txt"
testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
a = open(testfile, 'r')
al = a.readlines()
bl = []
Params_all = [] 
for each in al:
  bl.append(testdir+each.strip())
for each in bl: 
    a = pyfits.open(each)
    ydata = a[1].data
    ysigma = a[2].data
    testdata = np.zeros((npix, 1, 3))
    testdata[:, 0, 0] = xgrid1
    testdata[:, 0, 1] = ydata
    testdata[:, 0, 2] = ysigma
    continuum = get_continuum(testdata) 
    testdata[:, :, 1] /= continuum
    testdata[:, :, 2] /= continuum
    bad1 = np.isnan(testdata[:,:,1])
    testdata[bad1] = 0.
    bad2 = np.isnan(testdata[:,:,2])
    testdata[bad2] = np.Inf
    xdata = testdata[:,0,0]
    ydata = testdata[:,0,1]
    ysigma = testdata[:,0,2]
    coeffs_reshape = coeffs[:,-3:]
    ydata_norm = ydata  - coeffs[:,0] 
    coeffs_reshape = coeffs[:,-3:]
    Cinv = 1. / (ysigma ** 2 + scatters ** 2)
    MCM_rotate = np.dot(coeffs_reshape.T, Cinv[:,None] * coeffs_reshape)
    MCy_vals = np.dot(coeffs_reshape.T, Cinv * ydata_norm) 
    Params = linalg.solve(MCM_rotate, MCy_vals)
    Params = Params + [mean(metaall[:,0]) , mean(metaall[:,1]), mean(metaall[:,2])] # must be synchronized with fitspectra.py 
    print MCM_rotate, MCy_vals, Params 
    Params_all.append(Params) 

t_me = Params_all[:,0]
g_me = Params_all[:,1]
feh_me = Params_all[:,2]
testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
file2 = '4332_data_all_more.txt'
file2in = testdir+file2
t,g,feh,feh_err = loadtxt(file2in, usecols = (1,3,5,6), unpack =1) 

a = open('starsin_new_all_ordered.txt', 'r' )
al = a.readlines()
names = []
for each in al:
  names.append(each.split()[1])
unames = unique(names)
starind = arange(0,len(names), 1)
name_ind = []
names = array(names)
for each in unames:
  takeit = each == names
  name_ind.append(starind[takeit][-1]+1. )

name_ind_sort = sort(name_ind)
T_meta = metaall[:,0]
g_meta = metaall[:,1]
feh_meta = metaall[:,2]

