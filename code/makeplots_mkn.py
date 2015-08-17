import numpy as np
import matplotlib.pyplot as plt
import triangle

def hogg_savefig(fn, figure=None):
    print "hogg_savefig(): Writing %s ..." % fn
    if figure is None:
        plt.savefig(fn)
    else:
        figure.savefig(fn)
    print "hogg_savefig(): Wrote %s" % fn


def make_grid(xs, xmin, xmax, dx, ys, ymin, ymax, dy, data):
    young = 4. # Gyr
    # end magic number time
    xgrid = np.arange(xmin + 0.5 * dx, xmax, dx)
    nx = len(xgrid)
    ygrid = np.arange(ymin + 0.5 * dy, ymax, dy)
    ny = len(ygrid)
    xgrid = xgrid[:, None] * (np.ones(ny))[None, :]
    ygrid = (np.ones(nx))[:, None] * ygrid[None, :]
    Nstars = np.zeros_like(xgrid)
    young_fracs = np.zeros_like(xgrid) + np.NaN
    median_age = np.zeros_like(xgrid) + np.NaN
    median_age_mono = np.zeros_like(xgrid) + np.NaN
    median_age_low_alpha = np.zeros_like(xgrid) + np.NaN
    median_age_high_alpha = np.zeros_like(xgrid) + np.NaN
    median_feh = np.zeros_like(xgrid) + np.NaN
    median_alpha = np.zeros_like(xgrid) + np.NaN
    low_alpha_fracs = np.zeros_like(xgrid) + np.NaN
    std_age = np.zeros_like(xgrid) + np.NaN
    for xxi in range(nx):
        for yyi in range(ny):
            x, y = xgrid[xxi, yyi], ygrid[xxi, yyi]
            inside = ((xs > x - dx) *
                      (xs < x + dx) *
                      (ys > y - dy) *
                      (ys < y + dy))
            N = np.float(np.sum(inside))
            Nstars[xxi, yyi] = N
            median_age[xxi, yyi] = np.median(data[inside, 9]) 
            median_age_low_alpha_pick = data[inside, 8] < 0.08 - 0.15 * data[inside,8] 
            median_age_mid_feh_pick = logical_and(data[inside, 7] < 0, data[inside,7] > -0.2)
            median_age_mono_pick = logical_and(median_age_low_alpha_pick, median_age_mid_feh_pick)
            median_age_high_alpha_pick = data[inside, 8] >= 0.08 - 0.15 * data[inside,8] 
            median_age_low_alpha[xxi, yyi] = np.median(data[inside, 9][median_age_low_alpha_pick])
            median_age_high_alpha[xxi, yyi] = np.median(data[inside, 9][median_age_high_alpha_pick])
            median_age_mono[xxi, yyi] = np.median(data[inside, 9][median_age_mono_pick])
            std_age[xxi, yyi] = np.std(data[inside, 9]) 
            median_feh[xxi, yyi] = np.median(data[inside, 7]) 
            median_alpha[xxi, yyi] = np.median(data[inside, 8]) 
            if N > 8: # magic 8
                young_fracs[xxi, yyi] = np.sum(data[inside, 9] < young) / N
                low_alpha_fracs[xxi, yyi] = np.sum(data[inside, 8] <
                    0.05 - 0.2 * data[inside, 7]) / N # magic number insanity
    return xgrid, ygrid, Nstars, young_fracs, low_alpha_fracs, median_age,std_age, median_feh, median_alpha, median_age_low_alpha, median_age_high_alpha, median_age_mono

if __name__ == "__main__":
    print "Hello World!"

    fn = "../HWR_redclump_sample.txt"
    fn = "../Jon/redclump_sample_A_updatedvalues_only.txt"
    print "Reading %s ..." % fn
    data = np.genfromtxt(fn)
    names = ["ID", "distance", "Radius_gal", "Phi_gal", "z", "Teff",
             "logg", "[Fe/H]", "[alpha/Fe]", "age", "mass"]
    print data[2]
    print data.shape
    print "Read %s" % fn

    print "Making age histogram ..."
    plt.clf()
    plt.hist(data[:, -1], bins=32, histtype="step")
    plt.xlabel("age (Gyr)")
    hogg_savefig("age_hist.png")

    print "Making (R,z) grid..."
    # start magic number time
    Rmin, Rmax = 2., 16. # kpc
    Rmin, Rmax = 4., 15. # kpc
    zmin, zmax = -5.0, 5.0 # kpc
    zmin, zmax = -2.0, 2.0 # kpc
    dR = 0.1 # kpc
    (Rgrid, zgrid, Nstars, young_fracs,
     low_alpha_fracs, median_age,std_age, median_feh, median_alpha, median_age_low_alpha, median_age_high_alpha, median_age_mono) = make_grid(data[:, 2], Rmin, Rmax, dR,
                                  (data[:, 4]), zmin, zmax, dR, data)

    plt.figure(figsize=(12,6))
    plt.clf()
    imshow_kwargs = {"interpolation": "nearest",
                     "aspect": 1.0,
                     "origin": "lower",
                     "extent": (Rmin, Rmax, zmin, zmax)}
    plt.imshow(young_fracs.T, cmap=plt.cm.YlGnBu, vmin=0., vmax=1.,
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label("fraction of stars with low age")
    hogg_savefig("young_fracs_abs.png")

    plt.clf()
    plt.imshow(low_alpha_fracs.T, cmap=plt.cm.YlGnBu, vmin=0., vmax=1.,
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label("fraction of stars with low chemical age")
    hogg_savefig("low_alpha_fracs_abs.png")

    plt.clf()
    plt.imshow(Nstars.T, cmap=plt.cm.afmhot, vmin=0.,
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(r"number of stars in a $%.1f\times%.1f$ kpc box" %
                  (2. * dR, 2 * dR))
    hogg_savefig("nstars_abs.png")
    plt.clf()
    plt.imshow(median_age.T, cmap=plt.cm.YlGnBu, vmin=0.0, vmax = 13., 
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label("median age (Gyr)") 
    hogg_savefig("median_age_abs.png")
    plt.clf()
    plt.imshow(std_age.T, cmap=plt.cm.YlGnBu, vmin=0, vmax = 6., 
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(" age dispersion (Gyr)") 
    hogg_savefig("std_age_abs.png")
    plt.clf()
    plt.imshow(median_feh.T, cmap=plt.cm.YlGnBu, vmin=-0.5, vmax = 0.3, 
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(" median [Fe/H] (dex)") 
    hogg_savefig("median_feh_abs.png")
    plt.clf()
    plt.imshow(median_alpha.T, cmap=plt.cm.YlGnBu, vmin=-0.020, vmax = 0.20, 
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(" median [alpha/Fe] (dex)") 
    hogg_savefig("median_alpha_abs.png")
    plt.clf()
    plt.imshow(median_age_low_alpha.T, cmap=plt.cm.YlGnBu, vmin=0, vmax = 12, 
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(" median age (Gyr)") 
    hogg_savefig("median_age_low_alpha_abs.png")
    plt.clf()
    plt.imshow(median_age_high_alpha.T, cmap=plt.cm.YlGnBu, vmin=0, vmax = 12, 
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(" median age (Gyr)") 
    hogg_savefig("median_age_high_alpha_abs.png")
    plt.clf()
    plt.clf()
    plt.imshow(median_age_mono.T, cmap=plt.cm.YlGnBu, vmin=0, vmax = 12, 
               **imshow_kwargs)
    plt.xlabel("Galactocentric radius $R$ (kpc)")
    plt.ylabel("Galactic height $z$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(" median age (Gyr)") 
    hogg_savefig("median_age_mono_abs.png")
    plt.clf()




    print "Making (x,y) grid..."
    # start magic number time
    inplane = (data[:, 4] < 0.30) * (data[:, 4] > -0.30) # kpc
    xs = data[inplane, 2] * np.cos(data[inplane, 3])
    ys = data[inplane, 2] * np.sin(data[inplane, 3])
    xmin, xmax = 5., 14.
    ymin, ymax = -4.5, 4.5
    dR = 0.1 # kpc
    (xgrid, ygrid, Nstars, young_fracs,
     low_alpha_fracs,median_age, std_age, median_feh, median_alpha, median_age_low_alpha, median_age_high_alpha, median_age_mono) = make_grid(xs, xmin, xmax, dR,
                                  ys, ymin, ymax, dR, data[inplane])
    plt.figure(figsize=(5,4))
    plt.clf()
    imshow_kwargs = {"interpolation": "nearest",
                     "aspect": "equal",
                     "origin": "lower",
                     "extent": (xmin, xmax, ymin, ymax)}
    plt.imshow(young_fracs.T, cmap=plt.cm.YlGnBu, vmin=0., vmax=1.,
               **imshow_kwargs)
    plt.xlabel(r"Galactocentric $R\,\cos\phi$ (kpc)")
    plt.ylabel(r"Galactocentric $R\,\sin\phi$ (kpc)")
    foo = plt.colorbar()
    foo.set_label("fraction of stars with low age")
    hogg_savefig("fo_young_fracs_abs.png")

    plt.clf()
    plt.imshow(low_alpha_fracs.T, cmap=plt.cm.YlGnBu, vmin=0., vmax=1.,
               **imshow_kwargs)
    plt.xlabel(r"Galactocentric $R\,\cos\phi$ (kpc)")
    plt.ylabel(r"Galactocentric $R\,\sin\phi$ (kpc)")
    foo = plt.colorbar()
    foo.set_label("fraction of stars with low chemical age")
    hogg_savefig("fo_low_alpha_fracs_abs.png")

    plt.clf()
    plt.imshow(Nstars.T, cmap=plt.cm.afmhot, vmin=0.,
               **imshow_kwargs)
    plt.xlabel(r"Galactocentric $R\,\cos\phi$ (kpc)")
    plt.ylabel(r"Galactocentric $R\,\sin\phi$ (kpc)")
    foo = plt.colorbar()
    foo.set_label(r"number of stars in a $%.1f\times%.1f$ kpc box" %
                  (2. * dR, 2 * dR))
    hogg_savefig("fo_nstars_abs.png")
    plt.clf()
    imshow_kwargs = {"interpolation": "nearest",
                     "aspect": "equal",
                     "origin": "lower",
                     "extent": (xmin, xmax, ymin, ymax)}
    plt.imshow(median_age.T, cmap=plt.cm.YlGnBu, vmin=4., vmax=12.,
               **imshow_kwargs)
    plt.xlabel(r"Galactocentric $R\,\cos\phi$ (kpc)")
    plt.ylabel(r"Galactocentric $R\,\sin\phi$ (kpc)")
    foo = plt.colorbar()
    foo.set_label("")
    hogg_savefig("fo_median_age_abs.png")

    

#    print "Making triangle plot ..."
#    figure = triangle.corner(data[:, 1:], labels=names[1:])
#    hogg_savefig("triangle.png", figure=figure)

    print "Goodbye World!"
