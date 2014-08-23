            #continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar, 0.90)
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
            test = dataall[indx, star,2] == 0. 
            dataall[test] = 1. 
            ivar = 1. / (dataall[indx, star, 2] ** 2)
            test = ivar == 0.
            ivar[test] = 100.
            test = dataall == 0.
            dataall[test] = 0.0001
            ivar = array(ivar) 
            #bad = isinf(ivar)
            #ivar[bad] =1.
            continuum[ll, star] = weighted_median2(dataall[indx, star, 1], ivar, 0.90)
    return continuum

each = bl[1]
a = pyfits.open(each)
ydata = a[1].data
ysigma = a[2].data
testdata = np.zeros((npix, 1, 3))
testdata[:, 0, 0] = xgrid1
testdata[:, 0, 1] = ydata
testdata[:, 0, 2] = ysigma
continuum = get_continuum2(testdata) 
print continuum , median(continuum) 
