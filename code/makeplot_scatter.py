#!/usr/bin/python
import numpy
import pickle 
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
rcParams['figure.figsize'] = 12.0, 10.0

def plotdata(file_in, wl0, bw,prefix): 
    file_in2 = open(file_in, 'r') 
    dataall, metall, labels, offsets, coeffs, covs, scatters, chis, chisqs = pickle.load(file_in2)
    file_in2.close()

    rcParams['figure.figsize'] = 13.0, 10.0
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
      ax.plot(wl, factor*(val-sig), color=color, lw=lw) 
      ax.fill_between(wl, factor*(val+sig), factor*(val-sig), color = color, alpha = 0.2) 

      return None
  
    _plot_something(ax1, dataall[:, 0, 0], scatters, 0., 'k')
    for ax, indx, color, label in [(ax2, 0, "k", ""),
                                   (ax3, 1, "g", "Teff"),
                                   (ax3, 2, "b", "logg"),
                                   (ax3, 3, "r", "FeH")]:
      _plot_something(ax, dataall[:, 0, 0], coeffs[:, indx], covs[:, indx, indx], color, label=label)
    
    legend(numpoints=1)

    #ax3.legend(numpoints=1)
    ax1.set_ylim(0,np.median(scatters)*4.5) 
    ax2.set_ylim(0.7,1.2) 
    ax3.set_ylim(-0.3,0.3) 
    ax1.text(wl0-bw/2.+2., np.median(scatters)*3.5 , "scatter  " , fontsize = 12) 
    ax2.text(wl0-bw/2.+2., np.median(coeffs[:,0])*1.1, "mean spectra" , fontsize = 12) 
    ax2.text(wl0-bw/2.+2., np.median(coeffs[:,0])*1.1, "mean spectra" , fontsize = 12) 
    ax3.text(wl0-bw/2.+2., 0.2 , "[Fe/H]  coeff, log g coeff, Teff coeff * 1000 K" , fontsize = 12) 
    
    # attach lines to plots
    axlist = [ax1,ax2,ax3]
    line_kwargs = {"color": "k", "alpha": 0.25}
    for each in axlist:
      each.set_xlim(wl0-bw/2.,wl0+bw/2.) 
    ax2.axhline(1., **line_kwargs)
    ax1.axhline(dataall.shape[1], **line_kwargs)
    ax1.set_ylabel("scatter", fontsize = 20) 
    ax2.set_ylabel("coeff a0", fontsize = 20) 
    ax3.set_ylabel("coeff a1,a2,a3", fontsize = 20) 
    #ax1.semilogy()

    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    savefig(fig, prefix, transparent=False, bbox_inches='tight', pad_inches=0.5)
    return 

def savefig(fig, prefix, **kwargs):
    for suffix in (".eps", ".png"):
        print "writing %s" % (prefix + suffix)
        fig.savefig(prefix + suffix, **kwargs)

if __name__ == "__main__": #args in command line 
    wl1,wl2,wl3,wl4,wl5,wl6 = 15392, 15697, 15958.8, 16208.6, 16120.4, 16169.5 
    #plotdata('coeffs.pickle', wl3,100, "/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/R1_example") 
    #plotdata('coeffs_2nd_order.pickle', wl3,100, "/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/R1_example_2nd_order") 
    plotdata('coeffs_2nd_order.pickle', wl1,100, "/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/R1_example_2nd_order") 

