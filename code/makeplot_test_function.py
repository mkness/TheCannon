#!/usr/bin/python 
import scipy 
import numpy 
from numpy import * 
from scipy import ndimage
from scipy import interpolate 
from numpy import loadtxt
import os 
import numpy as np
from numpy import * 
from matplotlib import pyplot
import matplotlib.pyplot as plt
from matplotlib.pyplot import axes
from matplotlib.pyplot import colorbar
from matplotlib.ticker import NullFormatter
import pyfits
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

rc('text', usetex=True)
rc('font', family='serif')
plt.rcParams['xtick.major.pad'] = 9
rcParams["xtick.labelsize"] = 16
rcParams["ytick.labelsize"] = 16
rcParams['xtick.major.pad']='8'
rcParams['ytick.major.pad']='8'
#rcParams['figure.figsize'] = 10, 10
rcParams['figure.figsize'] = 12.0, 14.0


def plot_compare_params(params_labels, input_ASPCAP): 
  fig, temp = pyplot.subplots(3,1, sharex=False, sharey=False)
  ax1 = temp[0]
  ax2 = temp[1]
  ax3 = temp[2]
  g = input_ASPCAP[1,:]
  
  pick = logical_and(g > 0, logical_and(t_err < 300, feh > -4.0) ) 
  cind = array(input_ASPCAP[1][pick]) 
  listit_1 = [0,1,2]
  listit_2 = [1,0,0]
  cval = ['k', 'b', 'r'] 
  axs = [ax1,ax2,ax3]
  labels = ["ASPCAP log g", "ASPCAP Teff", "ASPCAP Teff"]
  for ax, num,num2,label1,x1,y1 in zip(axs, listit_1,listit_2,labels, [4800,3.0,0.3], [3400,1,-1.5]): 
    cind = array(input_ASPCAP[num2][pick]).flatten() 
    s1 = ax.scatter(input_ASPCAP[num][pick], params_labels[num][pick], c = cind, s = 30,alpha = 1.0, linewidths = 0 ,cmap = cm.jet ) 
    c_T = fig.colorbar(s1,ax=ax) 
    c_T.set_label(label1,fontsize = 14,labelpad = 10 ) 
    a,b,c1 = ax.errorbar(input_ASPCAP[num][pick], params_labels[num][pick],yerr= params_labels[num+3][pick],marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
    a,b,c2 = ax.errorbar(input_ASPCAP[num][pick], params_labels[num][pick], xerr=input_ASPCAP[num+3][pick],marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
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
  ax1.set_ylim(min(tme)-250, max(tme)+250) 
  ax2.set_ylim(round(min(gme)-0.2,1), round(max(gme)+0.2,1)) 
  ax3.set_ylim(min(fehme)-0.4, max(fehme)+0.4) 
  # attach lines to plots
  fig.subplots_adjust(hspace=0.22)
  #fig.savefig('/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/fits_all3_self.eps', transparent=True, bbox_inches='tight', pad_inches=0.2)
  #fig.savefig('/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/fits_all3_continuumcut.eps', transparent=True, bbox_inches='tight', pad_inches=0.2)
  #fig.savefig('/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/fits_all3_continuumcut2.eps', transparent=True, bbox_inches='tight', pad_inches=0.2)
  #fig.savefig('/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/fits_all3.eps', transparent=True, bbox_inches='tight', pad_inches=0.2)
  #fig.subplots_adjust(wspace=0.3)
  return 

#plot_compare_params(params_labels, input_ASPCAP) 
file3 = 'starsin_new_all_ordered.txt'
t,g,feh,t_err,feh_err = loadtxt(file3, usecols = (4,6,8,16,17), unpack =1) 
g_err = [0]*len(g) 
input_ASPCAP_cal = [t, g, feh, t_err,g_err, feh_err] 
file_in2 = open('datain_self.pickle', 'r') 
params_labels = pickle.load(file_in2)
file_in2.close()
params_labels = array(params_labels) 
input_ASPCAP_cal = array(input_ASPCAP_cal) 
plot_compare_params(params_labels, input_ASPCAP_cal) 
