#!/usr/bin/python
import numpy
import pickle 
from numpy import savetxt
import matplotlib
from matplotlib import pyplot
import scipy
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
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

def plotcutlines2(filein,ax,wl_cent,indx):
    wl,elem,offset1,offset2 = loadtxt(filein, usecols = (0,2,-2,-1), unpack =1,delimiter = ',') 
    a = open(filein,'r') 
    al = a.readlines() 
    elem_name = []
    for each in al:
        elem_name.append(each.split(',')[-3].strip() ) 
    elem_name = array(elem_name) 
    pick = logical_and(wl > wl_cent -13, wl < wl_cent + 13) 
    for a1,b1,c1,d1 in  zip(wl[pick],elem_name[pick],offset1[pick],offset2[pick]):
      if indx==1:
        ax.vlines(a1, c1,d1,linewidth=1.0, linestyle = 'dotted', color = 'grey' ,alpha=1.0)
        ax.text(a1, c1-0.03,b1, rotation='vertical',color = 'k', horizontalalignment= 'center',alpha=1.0,fontsize = 10) 
      if indx==0:
        ax.vlines(a1, -2,2,linewidth=1.0, linestyle = 'dotted', color = 'grey' ,alpha=1.0)
      if indx==3:
        ax.vlines(15770.9, -2,2.8,linewidth=1.5, linestyle = 'dashed', color = 'blue' ,alpha=0.5)
        ax.vlines(15769.4, -2,2.8,linewidth=1.5, linestyle = 'dashed', color = 'blue' ,alpha=0.5)
      if indx==4:
        ax.vlines(15770.9, -2,0.8,linewidth=1.5, linestyle = 'dashed', color = 'blue' ,alpha=0.5)
        ax.vlines(15769.4, -2,0.8,linewidth=1.5, linestyle = 'dashed', color = 'blue' ,alpha=0.5)
      draw() 
    return

def plotcutlines(filein,ax,wl_cent):
    wl,elem = loadtxt(filein, usecols = (0,2), unpack =1) 
    wl_A = wl*10.
    pick = logical_and(wl_A > wl_cent -10, wl_A < wl_cent + 10) 
    for a1,b1 in  zip(wl_A[pick],elem[pick]):
      b1 = np.int(b1) 
      strval = mydict[b1]
      ax.vlines(a1, 1.02,1.1,linewidth=1.0, linestyle = '-', color = 'k' ,alpha=0.5)
      ax.vlines(a1, 0.8,0.9,linewidth=1.0, linestyle = '-', color = 'k' ,alpha=0.5)
      ypick = 0.65+0.03*mod(b1,5.33) 
      ax.text(a1, 0.75,strval, rotation='vertical',color = 'k', horizontalalignment= 'center',alpha=1.0,fontsize = 10) 
      draw() 
      ax.vlines(15771,0,2)
      ax.vlines(15769,0,2)
    return

def plotdata(file_in, wl0, bw,prefix, cent_g1, cent_g2,indx_coeff): 
    file_in2 = open(file_in, 'r') 
    dataall, metall, labels, offsets, coeffs, covs, scatters,chis,chisqs = pickle.load(file_in2)
    file_in2.close()

    rcParams['figure.figsize'] = 14.0, 8.0
    #x, median_y, t_y, g_y,feh_y,chi_y = loadtxt('data_test.txt', usecols = (0,1,2,3,4,5), unpack =1) 
    fig, temp = pyplot.subplots(3,2, sharex=False, sharey=False)
    ax2 = temp[0,0]
    ax3 = temp[1,0]
    ax1 = temp[2,0]
    ax5 = temp[0,1]
    ax6 = temp[1,1]
    ax4 = temp[2,1]
    #ax2 = temp[0,0]
    #ax3 = temp[1,0]
    #ax5 = temp[0,1]
    #ax6 = temp[1,1]

    def _plot_something(ax, wl, val, var, color, lw=0.2, alpha= 0.8, alpha2=0.5,label=""):
      factor = 1.
      if label == "Teff_test "+r"($\times$ 1000)": factor = 1000. # magic numbers
      if label == "log Mass_test "+r"($\times$ 10)": factor = 10. # magic numbers
      sig = np.sqrt(var)
      ax.plot(wl, factor*(val+sig), color=color, lw=lw, alpha=alpha,label=label)
      ax.plot(wl, factor*(val-sig), color=color, lw=lw,alpha=alpha) 
      ax.fill_between(wl, factor*(val+sig), factor*(val-sig), color = color, alpha = alpha2) 

      return None
    
    if indx_coeff == 1:  
      file_gval = "redclump1_1pt5_2pt1_3pt3.txt"
      wl, g1,g2,g3 = loadtxt(file_gval, usecols = (0,1,2,3), unpack =1) 
      axlist = [ax1,ax4]
      for ax in axlist:
        ax.plot(wl,g1, 'k',linewidth = 0.5,linestyle = '-', label = "log g = 1.5")
        ax.plot(wl,g2, 'k',linewidth = 1,linestyle = 'dashed', label = "log g = 2.1")
        ax.plot(wl,g3, 'k',linewidth = 1, linestyle = '-', label = "log g = 3.3")
      leg1 = ax4.legend(numpoints=1, fontsize = 10,loc = 3,ncol = 5, frameon = False)
    if indx_coeff == 0:  
      file_gval = "redclump1_4100_4700_5100.txt"
      wl, g1,g2,g3 = loadtxt(file_gval, usecols = (0,1,2,3), unpack =1) 
      axlist = [ax1,ax4]
      for ax in axlist:
        ax.plot(wl,g1, 'k',linewidth = 0.5,linestyle = '-', label = "Teff = 4100K")
        ax.plot(wl,g2, 'k',linewidth = 1,linestyle = 'dashed', label = "Teff = 4700K")
        ax.plot(wl,g3, 'k',linewidth = 1, linestyle = '-', label = "Teff = 5100K")
      leg1 = ax4.legend(numpoints=1, fontsize = 10,loc = 3,ncol = 5, frameon = False)
    if indx_coeff == 2:  
      file_1 = "redclump1_m8_p0_p2.txt"
      file_2 = "redclump1_a0_a15_a25.txt"
      wl, f1,f2,f3 = loadtxt(file_1, usecols = (0,1,2,3), unpack =1) 
      wl, a1,a2,a3 = loadtxt(file_2, usecols = (0,1,2,3), unpack =1) 
      ax1.plot(wl,f1, 'k',linewidth = 0.5,linestyle = '-', label = "[Fe/H] = -0.8")
      ax1.plot(wl,f2, 'k',linewidth = 1,linestyle = 'dashed', label = "[Fe/H] = 0")
      ax1.plot(wl,f3, 'k',linewidth = 1, linestyle = '-', label = "[Fe/H] = 0.2")
      leg1 = ax1.legend(numpoints=1, fontsize = 10,loc = 3,ncol = 5, frameon = False)
      ax4.plot(wl,a1, 'k',linewidth = 0.5,linestyle = '-',label = "[alpha/Fe] = 0")
      ax4.plot(wl,a2, 'k',linewidth = 1,linestyle = 'dashed',label = "[alpha/Fe] = 0.15")
      ax4.plot(wl,a3, 'k',linewidth = 1, linestyle = '-', label = "[alpha/Fe] = 0.25")
      leg1 = ax4.legend(numpoints=1, fontsize = 10,loc = 3,ncol = 5, frameon = False)
    if indx_coeff == 3:  
      file_gval = 'redclump1_mpt7_m2pt1_m3pt3.txt'
      wl, g1,g2,g3 = loadtxt(file_gval, usecols = (0,1,2,3), unpack =1) 
      axlist = [ax1,ax4]
      for ax in axlist:
        ax.plot(wl,g1, 'k',linewidth = 0.5,linestyle = '-',label = "Mass = 0.7 Msun")
        ax.plot(wl,g2, 'k',linewidth = 1,linestyle = 'dotted',label = "Mass = 2.1 Msun")
        ax.plot(wl,g3, 'k',linewidth = 1, linestyle = '-', label = "Mass = 3.3 Msun")
      leg1 = ax4.legend(numpoints=1, fontsize = 10,loc = 3,ncol = 5, frameon = False)

    _plot_something(ax2, dataall[:, 0, 0], scatters+coeffs[:,0], 0., '#B0B0B0')
    _plot_something(ax2, dataall[:, 0, 0], -1.*scatters+coeffs[:,0], 0., '#B0B0B0')
    _plot_something(ax5, dataall[:, 0, 0], scatters+coeffs[:,0], 0., '#B0B0B0')
    _plot_something(ax5, dataall[:, 0, 0], -1.*scatters+coeffs[:,0], 0., '#B0B0B0')
    for ax, indx, color, label in [(ax2, 0, "k", ""),
            (ax5, 0, "k", "")]:
      _plot_something(ax, dataall[:, 0, 0], coeffs[:, indx] , covs[:, indx, indx], color, label=label)
     
    a1,a2,a3,a4,a5  = [0.4,0.4,0.4,0.4,0.4]
    avals= [a1,a2,a3,a4,a5] 
    avals2= [0.4,0.4,0.4,0.4,0.4] 
    avals[indx_coeff] = 1.0
    avals2[indx_coeff] = 1.0
    if indx_coeff !=  2: 
      for ax, indx, color, label,alval,alval2 in [(ax3, 1, "g", "Teff",avals[0],avals2[0]),
                                     (ax6, 1, "g", "Teff",avals[0],avals2[0]),
                                     (ax3, 2, "b", "logg",avals[1],avals2[1]),
                                     (ax6, 2, "b", "logg",avals[1],avals2[1]),
                                     (ax3, 3, "m", "[Fe/H]",avals[2],avals2[2]), 
                                     (ax6, 3, "m", "[Fe/H]",avals[2],avals2[2]),
                                     (ax3, 4, "c", r"[$\alpha$/Fe]",avals[2],avals2[2]), 
                                     (ax6, 4, "c",r"[$\alpha$/Fe]",avals[2],avals2[2]),
                                     (ax3, 5, "r","log Mass",avals[3],avals2[3]),
                                     (ax6, 5, "r", "log Mass",avals[3],avals2[3])]:
        # note this is ln Mass and not log mass 
          if indx_coeff == 1:
            _plot_something(ax, dataall[:, 0, 0], coeffs[:, indx]/max(abs(coeffs[0:4000,indx])) , covs[:, indx, indx]/max(abs(coeffs[:,indx])), color, alpha=alval,alpha2=alval2,label=label)
          else: 
            _plot_something(ax, dataall[:, 0, 0], coeffs[:, indx]/max(abs(coeffs[:,indx])) , covs[:, indx, indx]/max(abs(coeffs[:,indx])), color, alpha=alval,alpha2=alval2,label=label)
    if indx_coeff ==  2: 
      for ax, indx, color, label,alval,alval2 in [(ax3, 1, "g", "Teff",avals[0],avals2[0]),
                                     (ax6, 1, "g", "Teff",avals[0],avals2[0]),
                                     (ax3, 2, "b", "logg",avals[1],avals2[1]),
                                     (ax6, 2, "b", "logg",avals[1],avals2[1]),
                                     (ax3, 3, "m", "[Fe/H]",avals[2],avals2[2]), 
                                     (ax6, 3, "m", "[Fe/H]",avals[0],avals2[0]),
                                     (ax3, 4, "c", r"[$\alpha$/Fe]",avals[0],avals2[0]), 
                                     (ax6, 4, "c",r"[$\alpha$/Fe]",avals[2],avals2[2]),
                                     (ax3, 5, "r","log Mass",avals[3],avals2[3]),
                                     (ax6, 5, "r", "log Mass",avals[3],avals2[3])]:
        # note this is ln Mass and not log mass 
          if indx_coeff == 1:
            _plot_something(ax, dataall[:, 0, 0], coeffs[:, indx]/max(abs(coeffs[0:4000,indx])) , covs[:, indx, indx]/max(abs(coeffs[:,indx])), color, alpha=alval,alpha2=alval2,label=label)
          else: 
            _plot_something(ax, dataall[:, 0, 0], coeffs[:, indx]/max(abs(coeffs[:,indx])) , covs[:, indx, indx]/max(abs(coeffs[:,indx])), color, alpha=alval,alpha2=alval2,label=label)
    
    pick2 = abs(coeffs[:,2]) > 5*std(coeffs[:,2]) 
    
    plotcutlines2('keeplines3.txt',ax2,cent_g1,1) 
    plotcutlines2('keeplines3.txt',ax5,cent_g2,1) 
    plotcutlines2('keeplines3.txt',ax1,cent_g1,0) 
    plotcutlines2('keeplines3.txt',ax4,cent_g2,0) 
    plotcutlines2('keeplines3.txt',ax3,cent_g1,0) 
    plotcutlines2('keeplines3.txt',ax6,cent_g2,0) 
    if indx_coeff == 1:
      plotcutlines2('keeplines3.txt',ax2,cent_g2,4) 
      plotcutlines2('keeplines3.txt',ax5,cent_g2,4) 
      plotcutlines2('keeplines3.txt',ax1,cent_g2,3) 
      plotcutlines2('keeplines3.txt',ax4,cent_g2,3) 
      plotcutlines2('keeplines3.txt',ax3,cent_g2,3) 
      plotcutlines2('keeplines3.txt',ax6,cent_g2,3) 

    leg = ax6.legend(numpoints=1, fontsize = 10,loc = 3,ncol = 5, frameon = False)
    for legobj in leg.legendHandles:
      legobj.set_linewidth(2.0)
    if indx_coeff == 2: 
      leg2 = ax3.legend(numpoints=1, fontsize = 10,loc = 3,ncol = 5, frameon = False)
      for legobj in leg2.legendHandles:
        legobj.set_linewidth(2.0)
    if indx_coeff != 3: 
      for ax in [ax1,ax2,ax4,ax5]:
        ax.set_ylim(0.39,1.15) 
    if indx_coeff == 3: 
        ax1.set_ylim(0.75,1.05) 
        ax4.set_ylim(0.75,1.05) 
        ax2.set_ylim(0.39,1.15) 
        ax5.set_ylim(0.39,1.15) 
        draw()
    for ax in [ax3, ax6]: 
      #ax.get_xaxis().get_major_formatter().set_scientific(False) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
      ax.set_ylim(-1.15,1.15) 
    for ax in [ax1, ax4]: 
      #ax.get_xaxis().get_major_formatter().set_scientific(False) 
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    
    # attach lines to plots
    axlist1 = [ax1,ax2,ax3]
    axlist2 = [ax4,ax5,ax6]
    line_kwargs = {"color": "k", "alpha": 0.25}
    #cent_g1 = 15770
    #cent_g2 = 16810
    for each in axlist1:
      each.set_xlim(cent_g1-14,cent_g1+14) # logg1
      draw() 
    axlist2 = [ax4,ax5,ax6]
    for each in axlist2:
      each.set_xlim(cent_g2-14, cent_g2+14) 
      draw() 
    #for each in [ax2,ax5]:
    for each in [ax2,ax5,ax1,ax4]:
      each.axhline(1., **line_kwargs)
    for each in [ax3,ax6]:
      each.axhline(0., **line_kwargs)
    #for each in [ax2, ax4]: 
    #  each.axhline(dataall.shape[1], **line_kwargs)
    ax2.set_ylabel(r"${\theta_0}$", fontsize = 20) 
    #ax3.set_ylabel( r"${\theta_1}$,"+r"${\theta_2}$, "+r"${\theta_3}$, "+r"${\theta_4}$, "+r"${\theta_5}$" ,fontsize = 20) 
    ax3.set_ylabel( r"${\theta_l}$" ,fontsize = 20) 
    ax1.set_ylabel( "Normalized Flux" ,fontsize = 20) 
    ax1.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = 20,labelpad = 10) 
    ax4.set_xlabel("wavelength $\lambda$" + r" (\mbox{\AA})", fontsize = 20,labelpad = 10) 
    #r"$\dot{\Theta}$[deg/s]"
    ax2.yaxis.set_major_locator(ymajorLocator2)
    ax3.yaxis.set_major_locator(ymajorLocator3)
    #ax1.semilogy()
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax6.set_xticklabels([])
    ax5.set_xticklabels([])
    ax5.set_yticklabels([])
    ax6.set_yticklabels([])
    ax4.set_yticklabels([])
    ax3.set_yticks([-1,-0.5,0,0.5,1]) 

    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0.05)
    savefig_mkn(fig, prefix, transparent=False, bbox_inches='tight', pad_inches=0.5)
    return 

def savefig_mkn(fig, prefix, **kwargs):
    suffix1 =  ".png"
    #suffix2 =  ".pdf"
    print "writing %s" % (prefix + suffix1)
    fig.savefig(prefix + suffix1, **kwargs)
    #print "writing %s" % (prefix + suffix2)
    #fig.savefig(prefix + suffix2, **kwargs)

if __name__ == "__main__": #args in command line 
    wl1,wl2,wl3,wl4,wl5,wl6 = 15392, 15697, 15958.8, 16208.6, 16120.4, 16169.5 
    wl2 = 15700
    # g 
    #cent_g1 = 16030 # log g max - there is a spike here - need to raise this issue  - only arises when training on log mass 
    cent_g1 = 15770 # log g max
    cent_g2 = 16810 # log g 2nd 
    #cent_g2 = 16030 # bad terrible spike
    plotdata('coeffs_2nd_order_5.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_g_3", cent_g1, cent_g2,1) 
    # t 
    cent_g1 = 15339 # teff 2nd max 
    cent_g2 = 15720 # teff max 
    plotdata('coeffs_2nd_order_5.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_t_3", cent_g1, cent_g2,0) 
    # feh,alpha 
    cent_g1 = 15221.5 # highest feh
    cent_g2 = 16369.1 # highest alpha 
    plotdata('coeffs_2nd_order_5.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_af_3", cent_g1, cent_g2,2) 
    # mass 
    #cent_g1 = 16904 # highest mass
    cent_g1 = 15241 # highest mass for _5 and _5 HWR
    cent_g2 = 15332 # second highest mass for _5
    ##cent_g2 = 15432.5 # highest mass for _5HWR
    plotdata('coeffs_2nd_order_5.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_m_3", cent_g1, cent_g2,3) 
    #plotdata('coeffs_2nd_order_5HWR.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_m_3", cent_g1, cent_g2,3) 
    #cent_g1 = 15929
    #cent_g2 = 16410.7
    #plotdata('coeffs_2nd_order_5.pickle', wl3,100, "/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/coeffs_V_3", cent_g1, cent_g2,3) 


def plotlines2(ax):
    ax.plot(15770,1,'yo', markersize = 10) 
    draw() 
    return 
def plotlines(wavelength_centre,ax,coeffs,xgrid,s1,showvector, cent_g1, cent_g2):
    bw = 10.
    pick = logical_and(s1 < 50, logical_and(xgrid > wavelength_centre-bw+2, xgrid < wavelength_centre+bw-2) ) # find low-element-number lines in the plot range 
    for jj in arange(len(pick))[pick]:
        if showvector[jj]:
            strpick = s1[jj]
            strval = mydict[strpick]
            xpick = xgrid[jj]
            ypick = 0.65+0.03*mod(strpick,5.33) 
            ypick2 = 0.7#+0.03*mod(strpick,5.33) 
            ax.text(xpick, ypick,strval, rotation='vertical',color = element_colors[strval], horizontalalignment= 'center',alpha=1.0,fontsize = 10) 
            draw() 
    for a1,b1 in  zip(xgrid[pick],s1[pick]):
      strval = mydict[b1]
      ax.vlines(a1, 1.1,1.1,linewidth=1.0, linestyle = '-', color = element_colors[strval] ,alpha=0.5)
      ax.vlines(a1, 0.8,0.9,linewidth=1.0, linestyle = '-', color = element_colors[strval] ,alpha=0.5)
      draw()
    return

 
