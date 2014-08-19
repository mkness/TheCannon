#!/usr/bin/python
import numpy
from numpy import savetxt
import matplotlib
from matplotlib import pyplot
import scipy
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
s.set_family('serif')
s.set_size(14)
from matplotlib import rc
rc('text', usetex=False)
rc('font', family='serif')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
s.set_family('serif')
rcParams["xtick.labelsize"] = 14
rcParams["ytick.labelsize"] = 14
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
majorLocator   = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(5)
yminorLocator   = MultipleLocator(10)
yminorLocator2   = MultipleLocator(25)
xminorLocator   = MultipleLocator(5)
yminorLocator   = MultipleLocator(5)
ymajorLocator   = MultipleLocator(50)
xmajorLocator   = MultipleLocator(10)
rcParams['figure.figsize'] = 15.0, 10.0

#x, median_y, t_y, g_y,feh_y,chi_y = loadtxt('data_test.txt', usecols = (0,1,2,3,4,5), unpack =1) 
#fig1 = pyplot.figure()
#ax0 = fig1.add_subplot(111)
fig, ax = plt.subplots()

sortindx = 2
sortname = ["Teff", "logg", "Fe/H"]
index_use = argsort(metaall[:,sortindx])
ax.set_title("Per-pixel scaled residuals ($\chi$); spectra ordered by %s" % (sortname[sortindx]),fontsize = 20 ) 
ax.set_xlabel("Wavelength, $\AA$",fontsize = 20,labelpad = 10 ) 
ax.set_ylabel("Star Number",fontsize = 20) 
print "Ordered by %s" % (sortname[sortindx]) 

wl = dataall[:,0,0] 
image = np.arcsinh(chis) 
#image2 = np.insert(image[index_use].T, name_ind, values=-10, axis =1)
#test = ax.imshow(image[:,index_use].t, cmap=plt.cm.pink_r, interpolation="nearest", vmin = -5, vmax = 5 ,aspect = 'auto',origin = 'lower', extent = (wl.min(), wl.max(), 0, len(image.t))) 
test = ax.imshow(image[:,index_use].T, cmap=plt.cm.pink_r, interpolation="nearest", vmin = -5, vmax = 5 ,aspect = 'auto',origin = 'lower', extent = (wl.min(), wl.max(), 0, len(image.T))) 
cb = fig.colorbar(test) 
cb.set_label("arcsinh($\chi$)", fontsize = 20 ) 
