def weighted_median2(values, weights,quantile):
    """
    """
    sindx = np.argsort(values)
    cvalues = 1. * np.cumsum(weights[sindx])
    cvalues = cvalues / cvalues[-1]
    indx = (sindx[cvalues > quantile])[0]
    return values[indx]

def get_continuum2(dataall, delta_lambda=50):
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
    test = dataall[:,:,1] == 0
    dataall[test,1] = 1.
    continuum = np.zeros((Nlambda, Nstar))
    assert foo == 3
    for star in range(Nstar):
        print "get_continuum(): working on star" ,star
        for ll, lam in enumerate(dataall[:, 0, 0]):
            assert dataall[ll, star, 0] == lam
            indx = (np.where(abs(dataall[:, star, 0] - lam) < delta_lambda))[0]
            dataall_cut = dataall[:,star,2]
            ivar = 1. / (dataall[indx, star, 2] ** 2)
            #test = ivar == 0.
            #ivar[test] = 100.
            #ivar = array(ivar) 
            bad = isinf(ivar)
            ivar[bad] =1.
            #test = dataall == 0.
            #dataall[test] = 0.001
            continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar, 0.90)
    return continuum


each = bl2[0]
a = pyfits.open(each)
ydata = a[1].data
ysigma = a[2].data
testdata = np.zeros((npix, 1, 3))
testdata[:, 0, 0] = xgrid1
testdata[:, 0, 1] = ydata
testdata[:, 0, 2] = ysigma
continuum = get_continuum2(testdata) 
testdata[:, :, 1] /= continuum
testdata[:, :, 2] /= continuum
bad1 = np.isnan(testdata[:,:,1])
testdata[bad1,1] = 0.
#testdata[bad1] = 0.
bad2 = np.isnan(testdata[:,:,2])
testdata[bad2,2] = 100000.#np.Inf
#testdata[bad2] = 100000.#np.Inf
#bad3 = testdata[:,:,2] == 0.
#testdata[bad3,1] = 100000.
#testdata[bad1,1] = 0.
#testdata[bad2,1] = np.Inf
coeffs_reshape = coeffs[:,-3:]
ydata = testdata[:,0,1]
ysigma = testdata[:,0,2]
ydata_norm = ydata  - coeffs[:,0] 
coeffs_reshape = coeffs[:,-3:]
Cinv = 1. / (ysigma ** 2 + scatters ** 2)
MCM_rotate = np.dot(coeffs_reshape.T, Cinv[:,None] * coeffs_reshape)
MCy_vals = np.dot(coeffs_reshape.T, Cinv * ydata_norm) 
Params = linalg.solve(MCM_rotate, MCy_vals)
print Params 
