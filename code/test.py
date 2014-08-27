import numpy as np
import scipy.optimize as opt
#name = 'self_tags'
#file_in2 = open(name+'.pickle', 'r') 
#params, icovs_params = pickle.load(file_in2)
#file_in2.close()
#
#fn = 'coeffs.pickle'
#fd = open(fn, "r")
#dataall, metaall, labels, offsets, coeffs, covs, scatters = pickle.load(fd) 
#fd.close()
#fn = 'coeffs_2nd_order.pickle'
#fd = open(fn, "r")
#dataall, metaall, labels, offsets, coeffs, covs, scatters = pickle.load(fd) 


import numpy as np
import scipy.optimize as opt
#ytestdata = dataall[:,300,1] - coeffs[:,0] 
#y0 = 0.9 #whatever the star is 
nobservations = 554 
def main():
    nobservations = 554
    t,g,feh = metaall[:,0], metaall[:,1], metaall[:,2] 
    offsets = np.mean(metaall, axis=0)
    a, b, c = t[0]-offsets[0], g[0]-offsets[0], feh[0]-offsets[0]
    a, b, c = t[20]-offsets[0], g[20]-offsets[0], feh[20]-offsets[0]
    a, b, c = t[300], g[300], feh[300]
    #a, b, c = t[300], g[300], feh[300]
    #a, b, c = t[300], g[300], feh[300]
    #f, x1, x2, x3, x4, x5, x6, x7, x8, x9 = generate_data(nobservations, a, b, c)
    ytestdata = dataall[:,300,1] - coeffs[:,0] 
    sigmavals = dataall[:,300,2]
    #ytestdata = dataall[:,0,1] - coeffs[:,0] 
    #sigmavals = dataall[:,0,2]
    f = ytestdata 
    x0,x1,x2,x3,x4,x5,x6,x7,x8,x9 = coeffs[:,0], coeffs[:,1], coeffs[:,2], coeffs[:,3], coeffs[:,4], coeffs[:,5], coeffs[:,6] ,coeffs[:,7], coeffs[:,8], coeffs[:,9] 
    a, b, c = t[300], g[300], feh[300]
    #print 'Non-linear results (should be {}, {}, {}):'.format(a+offsets[0], b+offsets[0], c+offsets[0])
    print 'Non-linear results (should be {}, {}, {}):'.format(a,b, c)
    print nonlinear_invert(f, x1, x2, x3, x4, x5, x6, x7, x8, x9,sigmavals ) + offsets 


def func(x1, x2, x3, x4, x5, x6, x7, x8, x9, a, b, c):
    f = (0 
         + x1*a #x1 * a
         + x2*b #x2 * b
         + x3*c #x3 * c 
         + x4* a**2# #x4 * a**2 
         + x5 * a * b
         + x6 * a * c 
         + x7*b**2
         + x8  * b * c 
         + x9*c**2 )
         #+ x5 * a * b
         #+ x6 * a * c 
         #+ x7 * b**2  ) 
         #+ x8  * b * c 
         #+ x9 + c**2 )
         #-0)
    return f

def nonlinear_invert(f, x1, x2, x3, x4, x5, x6, x7, x8, x9,sigmavals ):
    # "curve_fit" expects the function to take a slightly different form...
    def wrapped_func(observation_points, a, b, c):
        x1, x2, x3, x4, x5, x6, x7, x8, x9  = observation_points
        return func(x1, x2, x3, x4, x5, x6, x7, x8, x9,  a, b, c)

    xdata = np.vstack([x1, x2, x3, x4, x5, x6, x7, x8, x9 ])
    model, cov = opt.curve_fit(wrapped_func, xdata, f, sigma = sigmavals)#,sigmavals)
    return model

main()

def generate_data(nobservations, a, b, c, noise_level=0.01):
    #noise = noise_level * np.random.normal(0, noise_level, nobservations)
    f = func(x1, x2, x3,x4, x5, x6, x7, x8, x9, a, b, c) #+ noise
    return f, x1, x2, x3,x4, x5, x6, x7, x8, x9


if False:
    def func(x1, x2, x3,x4,x5,x6,x7,x8,x9, a, b, c):
        f = (x1*a**2
             + x2 * b**2
             + x3 * c**2 
             + x4 * a * b 
             + x5 * a * c 
             + x6 * b * c 
             + x7 * a 
             + x8 * b 
             + x9 * c ) 
        return f

    x0,x1,x2,x3,x4,x5,x6,x7,x8,x9= coeffs[1,:]
    def nonlinear_invert(f, x1, x2,x3,x4,x5,x6,x7,x8,x9):
        def wrapped_func(observation_points, t, g, feh):
            x1, x2, x3, x4, x5, x6, x7, x8, x9 = observation_points
            return func(x1,x2,x3,x4,x5,x6,x7,x8,x9, t, g, feh)

        xdata = np.vstack([x1,x2,x3,x4,x5,x6,x7,x8,x9])
        model, cov = opt.curve_fit(wrapped_func, xdata, f)
        return model

    test = nonlinear_invert(0.98, x1, x2,x3,x4,x5,x6,x7,x8,x9)
