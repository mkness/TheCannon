delta_lambda=50
q = 0.5
print "get_continuum(): working on star" ,star
Nlambda, Nstar, foo = testdata2b.shape
continuum = np.zeros((Nlambda, 1))
bad = logical_or(np.isnan(testdata2b[:,star,1]), np.isinf(testdata2b[:,star,2]))
ratio_norm[bad] =1.
ratio_norm[bad] =0.
#testdata2b[bad,star,2] =200#np.Inf
#testdata2b[bad,star,1] 0#np.Inf
testdata2b[bad,star,2] =200#np.Inf
testdata2b[bad,star,1] = 0#np.Inf
for ll, lam in enumerate(testdata2b[:, 0, 0]):
    indx = (np.where(abs(testdata2b[:, star, 0] - lam) < delta_lambda))[0]
    #ivar = (1. / ((testdata2[:, star, 2] ** 2) )  + 1. / ((testdata2b[:, star, 2] ** 2) )  ) 
    ivar = 1./((testdata2b[:, star, 2] ** 2) )
    continuum[ll] = weighted_median(testdata2b[indx,star,1], ivar, q)
