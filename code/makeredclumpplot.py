#!/usr/bin/python
import pyfits 
import numpy as np
import pickle 
from numpy import savetxt
import matplotlib
from matplotlib import pyplot
#a.close() 
import scipy
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
s = matplotlib.font_manager.FontProperties()
s.set_family('serif')
s.set_size(14)
from matplotlib import rc
rc('text', usetex=True)
rc('text', usetex=True)
rc('font', family='serif')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import pyplot
from matplotlib.pyplot import * 
import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
import matplotlib as mpl 
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
rcParams["xtick.labelsize"] = 14
rcParams["ytick.labelsize"] = 14
s = matplotlib.font_manager.FontProperties()
s.set_size(18)
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
majorLocator   = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(5)
ymajorLocator1   = MultipleLocator(0.005)
ymajorLocator2   = MultipleLocator(0.1)
ymajorLocator3   = MultipleLocator(0.1)
xminorLocator   = MultipleLocator(5)
yminorLocator   = MultipleLocator(5)
ymajorLocator   = MultipleLocator(50)
xmajorLocator   = MultipleLocator(10)
######rcParams['figure.figsize'] = 22.0, 10.0

import pickle 
a = open('coeffs_2nd_order_5HWR_2.pickle', 'r') # was originally just coeffs_2nd_order_5.pickle but that is corrupted 
bc = pickle.load(a)
coeffs = bc[4]
wl_star = bc[0][:,0,0] 
a.close() 

def makespec(t,g,feh,alpha,mass): 
    nlabels = 5
    features_data = np.ones((nlabels, 1))
    labels = [t,g,feh,alpha,log(mass)] # put all in 
    offsets = [4.85481344e+03,   2.58309295e+00,   8.87884779e-04 , 4.4598135e-02, 1]
    labels = array(labels)
    offsets = array(offsets) 
    features_data = np.hstack((1, labels - offsets))
    newfeatures_data = np.array([np.outer(labels-offsets, labels-offsets)[np.triu_indices(nlabels)] ]) 
    newfeatures2_data = array([( labels[1] - offsets[1])**3]).T
    features_data_final = np.hstack((features_data, newfeatures_data.flatten(), newfeatures2_data.flatten()))
    jj= 0 
    model_gen2 = np.dot(coeffs,features_data_final.T)
    return model_gen2 

#filein = './Jon/redclump_sample_A_updatedvalues_only.txt'
#a = open(filein, 'r')
#al = a.readlines() 
##tvals,gvals,fehvals,alphavals, agevals, massvals = genfromtxt(filein, usecols = (5,6,7,8,9,10) , unpack =1) 
#ids = []
#for each in al:
#    ids.append(each.split()[0]) 


#plotstars(file1, wl_star,params1, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/spectra_fits", cent_wl) 
def plotstars(filein, wl, params,prefix,cent_wl1, cent_wl2,cent_wl3, cent_wl4,cent_wl5, cent_wl6, cent_wl7, cent_wl8, bw):
  rcParams['figure.figsize'] = 14.0, 14.0
  fig, temp = pyplot.subplots(4,2, sharex=False, sharey=False)
  ax1 = temp[0,0]
  ax2 = temp[1,0]
  ax3 = temp[2,0]
  ax4 = temp[3,0]
  ax5 = temp[0,1]
  ax6 = temp[1,1]
  ax7 = temp[2,1]
  ax8 = temp[3,1]
  t1,g1,feh1,alpha1,mass1 = params[0], params[1],params[2], params[3],params[4]
  param_str = " Teff="+str(np.int(t1))+", logg="+str(np.round(g1,1))+", [Fe/H]="+ str(np.round(feh1,2))+r", [$\alpha$/Fe]="+str(np.round(alpha1,2))+", Mass="+str(np.round(mass1,2))
  def _plot_something(filein,params, ax, wl, indx, color, lw=1.0, label=""):
      a = pyfits.open(filein) 
      data = a[1].data
      sigma = a[2].data
      bad = sigma > 0.1
      data[bad] = None
      model_aspcap = a[3].data
      t1,g1,feh1,alpha1,mass1 = params[0], params[1],params[2], params[3],params[4]
      model_tc = makespec(t1,g1,feh1,alpha1,mass1 )
      #if indx == 0: 
      #  ax.plot(wl,model_aspcap, color='r', lw=lw,alpha = 0.5) 
      lw2 = 0.1
      lw1 = 2.
      lw2 = 1.
      if indx == 1: 
        ax.plot(wl,model_aspcap, color='gray', lw=lw1,alpha = 0.6, label= "ASPCAP model", linestyle = 'dashed') 
        ax.plot(wl,model_tc, color='r', lw=lw1,alpha = 0.6, label = 'The Cannon model', linestyle = '-')
      ax.plot(wl,data , color='k', lw=lw2,label = 'data') #, label=label)
      return None

  axes = [ ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8] 
  for ax in axes: 
    _plot_something(filein, params,ax, wl_star, 1,'k')
  #_plot_something(filein, params,ax1, wl_star, 1,'k')
  #_plot_something(filein, params,ax3, wl_star,1, 'k')
  #_plot_something(filein, params,ax2, wl_star, 1,'k')
  #_plot_something(filein, params,ax4, wl_star,1, 'k')
  leg = ax1.legend(numpoints=1, fontsize = 12,loc = 3, frameon = False) 
  ax5.text(cent_wl5 - 13, 0.72, param_str)  
  #for legobj in leg.legendHandles:
  #  legobj.set_linewidth(2.0)
#  ax3.plot(wl_star, model_tc, 'b',alpha=0.5) 
  #fig.subplots_adjust(hspace=0)
  fig.subplots_adjust(wspace=0.09)
  fig.subplots_adjust(hspace=0.30)

  axlist = [ax1,ax2,ax3,ax4, ax5, ax6, ax7, ax8]
  axlist1 = [ax1,ax2,ax3,ax4]
  axlist2 = [ax3,ax4]
  fs = 16
  for ax in axlist1:
      ax.set_ylabel("Normalized Flux", fontsize = fs) 
  for ax in [ax1,ax1]:
      ax.set_xlim(cent_wl1-bw,cent_wl1+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  #for ax in axlist2:
  for ax in [ax2,ax2]:
      ax.set_xlim(cent_wl2-bw,cent_wl2+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  for ax in [ax3,ax3]:
      ax.set_xlim(cent_wl3-bw,cent_wl3+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  for ax in [ax4,ax4]:
      ax.set_xlim(cent_wl4-bw,cent_wl4+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  for ax in [ax5,ax5]:
      ax.set_xlim(cent_wl5-bw,cent_wl5+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  for ax in [ax6,ax6]:
      ax.set_xlim(cent_wl6-bw,cent_wl6+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  for ax in [ax7,ax7]:
      ax.set_xlim(cent_wl7-bw,cent_wl7+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  for ax in [ax8,ax8]:
      ax.set_xlim(cent_wl8-bw,cent_wl8+bw) # logg1
      ax.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = fs,labelpad = 5) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
  for ax in axlist:
      ax.set_ylim(0.69,1.05) # logg1
  ax1.set_ylim(0.49,1.05) # logg1
  #ax1.set_xticklabels([])
  #ax3.set_xticklabels([])
  ax5.set_yticklabels([])
  ax6.set_yticklabels([])
  ax7.set_yticklabels([])
  ax8.set_yticklabels([])
  #savefig_mkn(fig, prefix, transparent=False, bbox_inches='tight', pad_inches=0.5)
  return 


def savefig_mkn(fig, prefix, **kwargs):
    suffix1 =  ".png"
    print "writing %s" % (prefix + suffix1)
    fig.savefig(prefix + suffix1, **kwargs)


fn = '/Users/ness/new_laptop/Apogee_ages/Jon/aspcapStar-r5-v603-'
nd = '.fits'
file1,file2,file3,file4,file5 = fn+'2M00000211+6327470'+nd, fn+'2M00000446+5854329'+nd , fn+'2M00001242+5524391'+nd , fn+'2M23591255+5724354'+nd , fn+'2M00001962+6502429'+nd
file6,file7 =fn+'2M23591931+6244094'+nd, fn+'2M00010744+6445170'+nd
t1,g1,feh1,alpha1,age1, mass1 = 4708.49983734, 2.17652168722, 0.0617198884138, 0.015425129037, 1.55488236009, 1.88141910392
t2,g2,feh2,alpha2,age2, mass2 = 4745.52159963, 2.35896450654, -0.0142214905547, 0.0355349312863, 1.09579636191, 2.10351779915
t3,g3,feh3,alpha3,age3,mass3 =  4634.56555856, 2.27026176286, 0.148285193614, 0.0544413463398, 3.25464515338, 1.48482119829
t4,g4,feh4,alpha4,age4,mass4 = 5065.9005193, 2.65343988317, -0.272045304817, 0.0425327062973, 0.463644013007, 2.73168129836
t5,g5, feh5, alpha5, age5, mass5 = 4547.8021233, 2.36446969168, 0.149600284589, 0.0625720823962, 1.47742666725, 1.94314195429
t6,g6,feh6,alpha6,age6,mass6 = 4592.01662402, 2.37847512277, -0.0995352783245, 0.0849492206682 ,5.81821837085, 1.17919262914
t7,g7,feh7,alpha7,age7,mass7 = 4843.4531731, 2.49841383705, -0.0586235260442, 0.0386181490951, 10.5928384269, 1.00701251325
params1 = [4708,2.17,0.06,0.015,1.88] 
params4= [t4,g4,feh4,alpha4,mass4]
params7= [t7,g7,feh7,alpha7,mass7]

#logg 
cent_wl1 = 15770 # log g max
cent_wl5 = 16810 # log g 2nd 
# t 
cent_wl2 = 15339 # teff 2nd max 
cent_wl6 = 15720 # teff max 
# feh,alpha 
cent_wl3 = 15221.5 # highest feh
cent_wl7 = 16369 # highest alpha 
# mass
cent_wl4 = 15241 # highest mass for _5 and _5 HWR
cent_wl8 = 15332 # second highest mass for _5
# V
cent_wla = 15929 # highest mass for _5 and _5 HWR
cent_wlb = 16410.7 # second highest mass for _5
#plotdata('coeffs_2nd_order_5.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_t_3", cent_g1, cent_g2,0) 
    # feh,alpha 
    #cent_g1 = 15221.5 # highest feh
    #cent_g2 = 16369 # highest alpha 
    #plotdata('coeffs_2nd_order_5.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_af_3", cent_g1, cent_g2,2) 
    # mass 
#plotstars(file1, wl_star,params1, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/spectra_fits_1", cent_wl1,cent_wl2,cent_wl3,cent_wl4, cent_wl5, cent_wl6, cent_wl7, cent_wl8, 14) 
#plotstars(file7, wl_star,params7, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/spectra_fits_7", cent_wl1,cent_wl2,cent_wl3,cent_wl4,cent_wl5, cent_wl6, cent_wl7,cent_wl8, 14) 
plotstars(file1, wl_star,params1, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/spectra_fits_1", cent_wl1,cent_wl2,cent_wl3,cent_wl4,cent_wl5, cent_wl6, cent_wl7,cent_wl8, 14) 
#plotstars(file2, wl_star,params2, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/spectra_fits_2", cent_wl1,cent_wl2, cent_wl3, cent_wl4, cent_wl5, cent_wl6, cent_wl7,cent_wl8,14) 
plotstars(file7, wl_star,params7, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/spectra_fits_7", cent_wl1,cent_wl2,cent_wl3,cent_wl4,cent_wl5, cent_wl6, cent_wl7,cent_wl8, 14) 
#plotstars(file7, wl_star,params7, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/spectra_fits_elem", cent_wl1,cent_wl2,cent_wl3,cent_wl4,cent_wl5, cent_wl6, cent_wla,cent_wlb, 14) 
show() 
