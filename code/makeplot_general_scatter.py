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
from matplotlib import pyplot
from matplotlib.pyplot import * 
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

wl1,wl2,wl3,wl4,wl5,wl6 = 15392, 15697, 15958.8, 16208.6, 16120.4, 16169.5 
def plotdata(wl0, bw): 
    #x, median_y, t_y, g_y,feh_y,chi_y = loadtxt('data_test.txt', usecols = (0,1,2,3,4,5), unpack =1) 
    fig, temp = pyplot.subplots(3,1, sharex=True, sharey=False)
    ax1 = temp[0]
    ax2 = temp[1]
    ax3 = temp[2]

    def _plot_something(ax, wl, val, var, color, lw=2, label=""):
      factor = 1.
      if label == "Teff": factor = 1000. # yes, I feel dirty; MAGIC
      sig = np.sqrt(var)
      ax.plot(wl, factor*(val+sig), color=color, lw=lw, label=label)
      ax.plot(wl, factor*(val-sig), color=color, lw=lw, label=label)
      return None

    _plot_something(ax1, dataall[:, 0, 0], scatters, 0., 'k')
    for ax, indx, color, label in [(ax2, 0, "k", "continuum"),
                                   (ax3, 1, "g", "Teff"),
                                   (ax3, 2, "b", "logg"),
                                   (ax3, 3, "r", "FeH")]:
      _plot_something(ax, dataall[:, 0, 0], coeffs[:, indx], covs[:, indx, indx], color, label=label)

    ax3.legend(fontsize = 11,numpoints=1)
    ax1.set_ylim(0,median(scatters)*10.) 
    ax2.set_ylim(0.7,1.2) 
    ax3.set_ylim(np.median(coeffs[:,2])-15.0*np.std(coeffs[:,2][coeffs[:,2] < 1000]) ,np.median(coeffs[:,2])+15.0*np.std(coeffs[:,2][coeffs[:,2] < 1000])) 
    #ax1.text(wl0-bw/2.+2., np.median(chisqs)*10. , "chi2  " , fontsize = 12) 
    ax2.text(wl0-bw/2.+2., np.median(coeffs[:,0])*1.1, "mean spectra" , fontsize = 12) 
    ax3.text(wl0-bw/2.+2., np.median(coeffs[:,2])*10. , "[Fe/H]  coeff, log g coeff, Teff coeff*1000" , fontsize = 12) 
    
    # attach lines to plots
    axlist = [ax1,ax2,ax3]
    line_kwargs = {"color": "k", "alpha": 0.25}
    for each in axlist:
      each.axvline(wl0, **line_kwargs)
      each.axhline(0., **line_kwargs)
      each.axvspan(wl0-bw/2., wl0+bw/2., facecolor='c', alpha=0.1)
      each.set_xlim(wl0-bw/2.,wl0+bw/2.) 
    ax2.axhline(1., **line_kwargs)
    ax1.axhline(dataall.shape[1], **line_kwargs)
    ax1.set_ylabel("scatter", fontsize = 20) 
    ax2.set_ylabel("coeff a0", fontsize = 20) 
    ax3.set_ylabel("coeff a1,a2,a3", fontsize = 20) 
    #ax1.semilogy()

    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    return 

