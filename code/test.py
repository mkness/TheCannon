a = open("starsin_new_all_ordered.txt", 'r')
a = open("starsin_test.txt", 'r')
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
  if jj > 0:
    assert xgrid1[0] == dataall[0, 0, 0]
  dataall[:, jj, 0] = xgrid1
  dataall[:, jj, 1] = y2new
  dataall[:, jj, 2] = sigma
test = dataall[:,1,:] 
indx = (np.where(abs(dataall[:, star, 0] - lam) < delta_lambda))[0]
ivar = 1. / (dataall[indx, star, 2] ** 2)
ivar = array(ivar) 
bad = isinf(ivar)
ivar[bad] = 1.0
continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar, 0.90)

