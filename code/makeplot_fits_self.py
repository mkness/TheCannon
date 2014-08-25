#!/usr/bin/python 
import scipy 
import numpy 
import pickle
from numpy import * 
from scipy import ndimage
from scipy import interpolate 
from numpy import loadtxt
import os 
import numpy as np
from numpy import * 
import matplotlib 
from pylab import rcParams
from pylab import * 
from matplotlib import pyplot
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.pyplot import axes
from matplotlib.pyplot import colorbar
#from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
s.set_family('serif')
s.set_size(14)
from matplotlib import rc
rc('text', usetex=False)
rc('font', family='serif')
#rcParams["xtick.labelsize"] = 14
#rcParams["ytick.labelsize"] = 14
#rcParams['figure.figsize'] = 12.0, 10.0

def plotfits():
    #file_in = "4332_tags.pickle" #4332
    file_in = "self_tags.pickle"
    file_in2 = open(file_in, 'r') 
    params, icovs_params = pickle.load(file_in2)
    params = array(params)
    file_in2.close()

    #testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
    #file2 = '4332_data_all_more.txt'
    #file2in = testdir+file2
    #t,t_err,g,g_err,feh,feh_err = loadtxt(file2in, usecols = (1,2,3,4,5,6), unpack =1) 
     
    filein2 = 'starsin_test2.txt' # this is for self test this is dangerous - need to implement the same logg cut here, this is original data values or otherwise metadata 
    filein2 = 'starsin_new_all_ordered.txt' # this is for self test this is dangerous - need to implement the same logg cut here, this is original data values or otherwise metadata 
    t,g,feh,t_err,feh_err = loadtxt(filein2, usecols = (4,6,8,16,17), unpack =1) 
    g_err = [0]*len(g) 
    g_err = array(g_err)
    #g_cut, teff_cut = 40, 4000
    #vals_in = [t, g, feh, t_err, g_err, feh_err] 
    #cut = logical_and(g < g_cut, t > teff_cut) 
    #for each in vals_in:

    #file_in = "coeffs.pickle"
    #file_in2 = open(file_in, 'r') 
    #dataall, metall, labels, offsets, coeffs, covs, scatters = pickle.load(file_in2) - but need apogee errors on these so can only use the metadata for the 3 main t,g,feh - should also save apogee
    #errors 
    #file_in2.close()

    params = array(params) 
    covs_params = np.linalg.inv(icovs_params) 
    rcParams['figure.figsize'] = 12.0, 10.0
    fig, temp = pyplot.subplots(3,1, sharex=False, sharey=False)
    ax1 = temp[0]
    ax2 = temp[1]
    ax3 = temp[2]
    params_labels = [params[:,0], params[:,1], params[:,2] , covs_params[:,0,0]**0.5, covs_params[:,1,1]**0.5, covs_params[:,2,2]**0.5 ]  
    pick = logical_and(g > 0, logical_and(t_err < 300, feh > -4.0) ) 
    cval = ['k', 'b', 'r'] 
    input_ASPCAP = [t, g, feh, t_err, g_err, feh_err] 
    cind = array(input_ASPCAP[1][pick]) 
    listit_1 = [0,1,2]
    listit_2 = [1,0,0]
    axs = [ax1,ax2,ax3]
    labels = ["ASPCAP log g", "ASPCAP Teff", "ASPCAP Teff"]
    for ax, num,num2,label1,x1,y1 in zip(axs, listit_1,listit_2,labels, [4800,3.0,0.3], [3400,1,-1.5]): 
      cind = array(input_ASPCAP[num2][pick]).flatten() 
      s1 = ax.scatter(input_ASPCAP[num][pick], params_labels[num][pick], c = cind, s = 30,alpha = 1.0, linewidths = 0 ,cmap = cm.jet ) 
      c_T = fig.colorbar(s1,ax=ax) 
      c_T.set_label(label1,fontsize = 14,labelpad = 10 ) 
      a,b,c1 = ax.errorbar(input_ASPCAP[num][pick], params_labels[num][pick],yerr= params_labels[num+3][pick],marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
      a,b,c2 = ax.errorbar(input_ASPCAP[num][pick], params_labels[num][pick],xerr=input_ASPCAP[num+3][pick],marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
      g_color = c_T.to_rgba(cind) 
      c1[0].set_color(g_color)
      c2[0].set_color(g_color)
      ax.text(x1,y1,"y-axis, $<\sigma>$ = "+str(round(mean(params_labels[num+3][pick]),2)),fontsize = 14) 

    ax1.plot([0,6000], [0,6000], linewidth = 1.5, color = 'k' ) 
    ax2.plot([0,5], [0,5], linewidth = 1.5, color = 'k' ) 
    ax3.plot([-3,2], [-3,2], linewidth = 1.5, color = 'k' ) 
    ax1.set_xlim(3500, 5500) 
    ax2.set_xlim(0, 5) 
    ax3.set_xlim(-3, 2) 
    ax1.set_xlabel("ASPCAP Teff, [K]", fontsize = 14,labelpad = 5) 
    ax1.set_ylabel("NHR+ Teff, [K]", fontsize = 14,labelpad = 5) 
    ax2.set_xlabel("ASPCAP logg, [dex]", fontsize = 14,labelpad = 5) 
    ax2.set_ylabel("NHR+ logg, [dex]", fontsize = 14,labelpad = 5) 
    ax3.set_xlabel("ASPCAP [Fe/H], [dex]", fontsize = 14,labelpad = 5) 
    ax3.set_ylabel("NHR+ [Fe/H], [dex]", fontsize = 14,labelpad = 5) 
    ax1.set_ylim(min(params_labels[0][pick])-250, max(params_labels[0][pick])+250) 
    ax2.set_ylim(round(min(params_labels[1][pick])-0.2,1), round(max(params_labels[1][pick])+0.2,1)) 
    ax3.set_ylim(min(params_labels[2][pick])-0.4, max(params_labels[2][pick])+0.4) 
    # attach lines to plots
    fig.subplots_adjust(hspace=0.22)
    #prefix = "/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/fits_3_self_cut"
    prefix = "/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/test_self"
    savefig(fig, prefix, transparent=False, bbox_inches='tight', pad_inches=0.5)
    return 

def savefig(fig, prefix, **kwargs):
    for suffix in (".eps", ".png"):
        print "writing %s" % (prefix + suffix)
        fig.savefig(prefix + suffix, **kwargs)

if __name__ == "__main__": #args in command line 
    wl1,wl2,wl3,wl4,wl5,wl6 = 15392, 15697, 15958.8, 16208.6, 16120.4, 16169.5 
    plotfits() 

