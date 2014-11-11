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

def returnscatter(x,y):
    xd = x
    yd = y
    # sort the data
    reorder = sorted(range(len(xd)), key = lambda ii: xd[ii])
    xd = [xd[ii] for ii in reorder]
    yd = [yd[ii] for ii in reorder]
    # determine best fit line
    slope=1.
    intercept= -4.
    xl = [min(xd), max(xd)]
    yl = [slope*xx + intercept  for xx in xl]
    # coefficient of determination, plot text
    variance = np.var(yd)
    residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(xd,yd)])
    #residuals_test = (np.sum([(slope*xx + intercept - yy)**2  for xx,yy in zip(xd,yd)])/len(xd))**0.5
    rms = (np.sum([ (xx - yy)**2  for xx,yy in zip(xd,yd)])/len(xd))**0.5
    bias = (np.mean([ (xx - yy)  for xx,yy in zip(xd,yd)]))
    Rsqr = np.round(1-residuals/variance, decimals=2)
    # error bounds
    yerr = [abs(slope*xx + intercept - yy)  for xx,yy in zip(xd,yd)]
    par = np.polyfit(xd, yerr, 2, full=True)
    yerrUpper = [(xx*slope+intercept)+(par[0][0]*xx**2 + par[0][1]*xx + par[0][2]) for xx,yy in zip(xd,yd)]
    yerrLower = [(xx*slope+intercept)-(par[0][0]*xx**2 + par[0][1]*xx + par[0][2]) for xx,yy in zip(xd,yd)]
    return bias, rms#**0.5

def plotfits():
#    file_in = "self_tags.pickle"
    file_in = "self_2nd_order_tags.pickle"
    file_in2 = open(file_in, 'r') 
    params, covs_params,variable1,ids = pickle.load(file_in2)
    #params, icovs_params,chi2,ids = pickle.load(file_in2)
    params = array(params)
    file_in2.close()

    filein2 = 'test18.txt' # this is for self test this is dangerous - need to implement the same logg cut here, this is original data values or otherwise metadata 
    #filein2 = 'mkn_labels_Atempfeh_edit.txt'  # this is for using all stars ejmk < 0.3 but with offest to aspcap values done in a consistent way to rest of labels 
    a = open(filein2) 
    al = a.readlines() 
    names = []
    for each in al:
      names.append(each.split()[1]) 
    unames = unique(names) 
    starind = arange(0,len(names), 1) 
    name_ind = [] 
    names = array(names) 
    for each in unames:
      takeit = each == names 
      name_ind.append(np.int(starind[takeit][-1]+1. ) )
    cluster_ind = [0] + list(sort(name_ind))# + [len(al)]
    plot_markers = ['ko', 'yo', 'ro', 'bo', 'co','k*', 'y*', 'r*', 'b*', 'c*', 'ks', 'rs', 'bs', 'cs', 'rd', 'kd', 'bd', 'cd', 'mo', 'ms' ]
    cluster_names = names[list([a-1 for a in cluster_ind])][1:]
    #plot_markers = ['k', 'y', 'r', 'b', 'c','k', 'y', 'r', 'b', 'c', 'k', 'r', 'b', 'c', 'r', 'k', 'b', 'c', 'm', 'm' ]
    #cv_ind = np.arange(395,469,1)
    #a  = open(filein2)
    #al = a.readlines() 
    #bl = []
    #for each in al:
    #    bl.append(each.strip())
    #bl = np.delete(bl, [cv_ind], axis = 0) 
    #savetxt("starsin_cut.txt", bl, fmt = "%s") 
    #filein3 = 'starsin_cut.txt'
    if filein2 == 'test18.txt':
      t,g,feh,t_err,feh_err = loadtxt(filein2, usecols = (4,6,8,16,17), unpack =1) 
    if filein2 == 'mkn_labels_Atempfeh_edit.txt':
      t,g,feh,t_err,feh_err = loadtxt(filein2, usecols = (3,5,7,3,3), unpack =1) 
    g_err = [0]*len(g) 
    feh_err = [0]*len(g) 
    t_err = [0]*len(g) 
    g_err = array(g_err)
    t_err = array(t_err)
    feh_err = array(feh_err)
    
    params = array(params) 
    #covs_params = np.linalg.inv(icovs_params) 
    rcParams['figure.figsize'] = 12.0, 10.0
    fig = plt.figure()
    ax = fig.add_subplot(111,frameon = 0 )
    ax.set_ylabel("The Cannon: output labels",labelpad = 40, fontsize = 20 ) 
    ax.set_xlabel("ASPCAP input labels",labelpad = 30, fontsize = 20 ) 
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax1 = fig.add_subplot(311)#,sharey = None)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    params_labels = [params[:,0], params[:,1], params[:,2] , covs_params[:,0,0]**0.5, covs_params[:,1,1]**0.5, covs_params[:,2,2]**0.5 ]  
    cval = ['k', 'b', 'r'] 
    input_ASPCAP = [t, g, feh, t_err, g_err, feh_err] 
    listit_1 = [0,1,2]
    listit_2 = [1,0,0]
    axs = [ax1,ax2,ax3]
    labels = ["ASPCAP log g", "ASPCAP Teff", "ASPCAP Teff"]
    names =array(names) 
    for i in range(0,len(cluster_ind)-1): 
      indc1 = cluster_ind[i]
      indc2 = cluster_ind[i+1]
      name = names[i] 
      for ax, num,num2,label1,x1,y1 in zip(axs, listit_1,listit_2,labels, [4800,3.0,0.3], [3400,1,-1.5]): 
        pick = logical_and(g[indc1:indc2] > 0, logical_and(t_err[indc1:indc2] < 300, feh[indc1:indc2] > -4.0) ) 
        cind = array(input_ASPCAP[1][indc1:indc2][pick]) 
        cind = array(input_ASPCAP[num2][indc1:indc2][pick]).flatten() 
        #ax.errorbar(input_ASPCAP[num][indc1:indc2][pick], params_labels[num][indc1:indc2][pick],yerr= params_labels[num+3][indc1:indc2][pick],marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
        #ax.errorbar(input_ASPCAP[num][indc1:indc2][pick], params_labels[num][indc1:indc2][pick],xerr=input_ASPCAP[num+3][indc1:indc2][pick],marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
        #ax.plot(input_ASPCAP[num][indc1:indc2][pick], params_labels[num][indc1:indc2][pick], plot_markers[i], label = str(name)) 
        ax.plot(input_ASPCAP[num][indc1:indc2][pick], params_labels[num][indc1:indc2][pick], plot_markers[i], label = str(cluster_names[np.int(i)]), ms = 6)
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.45),fancybox=True, shadow=True, numpoints = 1, ncol = 7,fontsize  = 12 ) 
    for ax, num in zip(axs, listit_1): 
        scatter1,scatter2 = returnscatter(input_ASPCAP[num][pick], params_labels[num][pick])
        if abs(scatter2) > 1:
          round2 = round(scatter2, 1)
        else: 
          round2 = round(scatter2, 2)
        if abs(scatter1) > 1:
          round1 = round(scatter1, 1)
        else: 
          round1 = round(scatter1, 2)
        if num == 0: 
          scatter3 = round(mean(params_labels[num+3][pick]),2) 
        else: 
          scatter3 = round(mean(params_labels[num+3][pick]),3) 
        if abs(scatter3) > 1:
          round3 = round(scatter3, 3)
        else: 
          round3 = round(scatter3, 3)
        ax.text(0.7, 0.1,"bias, rms, precision = "+str(round1)+", "+str(round2)+", "+str(round3), ha='center', va='center',
              transform=ax.transAxes,fontsize = 14)
          

    #plot_markers = ['ko', 'yo', 'ro', 'bo', 'co','k*', 'y*', 'r*', 'b*', 'c*', 'ks', 'rs', 'bs', 'cs', 'rd', 'kd', 'bd', 'cd', 'mo', 'ms' ]
    #for a,b in zip(plot_markers, cluster_names):
    #    ax1.plot(0,0, a,label = b) 
    #ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.35),fancybox=True, shadow=True, numpoints = 1, ncol = 10) 
    #ax1.legend() 
    
    #for ax,num in zip([ax1,ax2,ax3], [0,1,2]):
    #ax1.text(0.7, 0.1,"y-axis, $<\sigma>$ = "+str(round(mean(params_labels[0+3][pick]),2)), transform = ax1.transAxes, fontsize = 14) 
    #ax2.text(0.7, 0.1,"y-axis, $<\sigma>$ = "+str(round(mean(params_labels[1+3][pick]),3)), transform = ax2.transAxes, fontsize = 14) 
    #ax3.text(0.7, 0.1,"y-axis, $<\sigma>$ = "+str(round(mean(params_labels[2+3][pick]),3)), transform = ax3.transAxes, fontsize = 14) 
    #ax.text(0.7, 0.1,"bias, rms, precision = "+str(round1)+", "+str(round2)+", "+str(round3), ha='center', va='center',

    ax1.plot([0,6000], [0,6000], linewidth = 1.5, color = 'k' ) 
    ax2.plot([0,5], [0,5], linewidth = 1.5, color = 'k' ) 
    ax3.plot([-3,2], [-3,2], linewidth = 1.5, color = 'k' ) 
    ax1.set_xlim(3500, 6000) 
    ax1.set_ylim(3500, 6000) 
    ax2.set_xlim(0, 5) 
    ax2.set_xlim(0, 5) 
    ax1.set_xlabel(" Teff, [K]", fontsize = 14,labelpad = 2) 
    ax1.set_ylabel("Teff, [K]", fontsize = 14,labelpad = 5) 
    ax2.set_xlabel(" logg, [dex]", fontsize = 14,labelpad = 2) 
    ax2.set_ylabel("logg, [dex]", fontsize = 14,labelpad = 5) 
    ax3.set_xlabel(" [Fe/H], [dex]", fontsize = 14,labelpad = 2) 
    ax3.set_ylabel("[Fe/H], [dex]", fontsize = 14,labelpad = 5) 
    ax3.set_ylim(-3,2) 
    # attach lines to plots
    fig.subplots_adjust(hspace=0.22)
    #prefix = "/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/fits_3_self_cut"
#    prefix = "/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/test_self"
#    savefig(fig, prefix, transparent=False, bbox_inches='tight', pad_inches=0.5)
    return 

def savefig(fig, prefix, **kwargs):
    for suffix in (".eps", ".png"):
        print "writing %s" % (prefix + suffix)
        fig.savefig(prefix + suffix, **kwargs)

if __name__ == "__main__": #args in command line 
    wl1,wl2,wl3,wl4,wl5,wl6 = 15392, 15697, 15958.8, 16208.6, 16120.4, 16169.5 
    plotfits() 

