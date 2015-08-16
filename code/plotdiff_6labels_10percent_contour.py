#!/usr/bin/python 
import scipy 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import pickle 
import matplotlib.cm as cm 
import pyfits 
import glob 
from glob import glob 
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
rcParams["xtick.labelsize"] = 13
rcParams["ytick.labelsize"] = 13
rcParams['figure.figsize'] = 24.0, 3.5
#f = plt.figure()
#ax = f.add_subplot(111,frameon = 0 )
#ax.set_xlabel("Input Labels",labelpad = 40, fontsize = 20 )
#ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#ax1 = f.add_subplot(111,sharey = None,sharex = None)
#ax2 = f.add_subplot(121,sharey = None,sharex = None)
#ax3 = f.add_subplot(131,sharey = None,sharex = None)
#ax3 = f.add_subplot(132,sharey = None)
#ax4 = f.add_subplot(143,sharey = None)
#ax2 = f.add_subplot(112)
#ax3 = f.add_subplot(113)
#ax4 = f.add_subplot(114)
#ax5 = f.add_subplot(115)
#ax6 = f.add_subplot(116)
#    ax1 = fig.add_subplot(311)#,sharey = None)
#    ax2 = fig.add_subplot(312)
#    ax3 = fig.add_subplot(313)

f, (ax1, ax2, ax3,ax4,ax5,ax6) = plt.subplots(ncols=6)
ax = f.add_subplot(111,frameon = 0 )
ax.set_xlabel("Input Labels",labelpad = 5, fontsize = 20 )
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#f, (ax1, ax2, ax3,ax4,ax5,ax6) = plt.subplots(ncols=6)
f2, (ax7,ax8) = plt.subplots(ncols=2)


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
    variance = np.var(yd)
    residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(xd,yd)])
    rms = (np.sum([ (xx - yy)**2  for xx,yy in zip(xd,yd)])/len(xd))**0.5
    bias = (np.mean([ (xx - yy)  for xx,yy in zip(xd,yd)]))
    Rsqr = np.round(1-residuals/variance, decimals=2)
    yerr = [abs(slope*xx + intercept - yy)  for xx,yy in zip(xd,yd)]
    par = np.polyfit(xd, yerr, 2, full=True)
    yerrUpper = [(xx*slope+intercept)+(par[0][0]*xx**2 + par[0][1]*xx + par[0][2]) for xx,yy in zip(xd,yd)]
    yerrLower = [(xx*slope+intercept)-(par[0][0]*xx**2 + par[0][1]*xx + par[0][2]) for xx,yy in zip(xd,yd)]
    return bias, rms#**0.5

# read in files 
dir1 = '/Users/ness/new_laptop/Apogee_ages/crossvalidation/' 
starsout_files = glob(dir1+"random_not*log*")
training_files =  glob(dir1+"training*5_*_log.pickle*")
fn2 = dir1+'apokasc_ages_actual.txt'
fn = dir1+'training_apokasc_ages.list'
tin,gin,fehin,alphain, T_A, g_A, feh_A,rc_est = np.loadtxt(fn, usecols = (1,2,3,4,1,2,3,5), unpack =1)
massin = loadtxt(fn2, usecols = (-2,), unpack =1) 
snr = loadtxt(dir1+"training_apokasc_ages_SNR.list", usecols = (-1,), unpack =1) 

numlist = arange(0,10,1) 
t,g,feh,tout,gout,fehout = [],[],[],[],[],[]
alpha,alphaout,mass,ln_massout,snrout = [],[],[],[], []
for num in numlist: 
    a = open(training_files[num], 'r') 
    list_a = loadtxt(starsout_files[num], usecols = (0,), unpack =1) 
    list_a = list(list_a) 
    b = pickle.load(a) 
    a.close()
    tout1,ln_massout1,gout1,fehout1,alphaout1 = b[0][:,0], b[0][:,4], b[0][:,1], b[0][:,2], b[0][:,3]
    tout.append(tout1[list_a]) 
    gout.append(gout1[list_a]) 
    fehout.append(fehout1[list_a]) 
    alphaout.append(alphaout1[list_a]) 
    ln_massout.append(ln_massout1[list_a]) 
    mass.append( massin[list_a]) 
    t.append(tin[list_a]) 
    g.append( gin[list_a]) 
    feh.append( fehin[list_a]) 
    alpha.append( alphain[list_a]) 
    snrout.append( snr[list_a]) 

t,tout,g,gout = hstack((t)), hstack((tout)), hstack((g)), hstack((gout)) 
feh,fehout,alpha,alphaout, mass, ln_massout =  hstack((feh)), hstack((fehout)), hstack((alpha)), hstack((alphaout)), hstack((mass)), hstack((ln_massout) ) 
snrout = hstack((snrout))
vmin1,vmax1  = 1, 3.4
var = gout 
vmin1,vmax1  = -0.03, 0.3

cval = cm.cool
msize = 20.

# convert log mass to log age
plotage = False
plotage = True
if plotage:
    import masstoage_interp
    from masstoage_interp import *
    age = get_age_from_mass(feh,mass)
    ageout = get_age_from_mass(fehout,np.exp(ln_massout))

if plotage: 
  varin  =  [ t,g,feh,alpha,log10(mass), log10(age)]
  varout  =  [ tout,gout,fehout,alphaout,log10(e**ln_massout), log10(ageout)] 
else:
  varin  =  [ t,g,feh,alpha,log10(mass), log10(mass)]
  varout  =  [ tout,gout,fehout,alphaout,log10(e**ln_massout),log10(e**ln_massout)  ]
axall = [ax1,ax2,ax3,ax4,ax5,ax6]
for a,b,ax in zip(varin,varout,axall):
    counts,ybins,xbins,image = ax.hist2d(a,b,bins=30,norm=LogNorm(),cmin=1,cmap = cm.gray_r )# ,normed=True)
    test1,test2,test3,test4 = ax.hist2d(varin[0],varin[1],bins=30,norm = LogNorm(), cmin=1,cmap = cm.gray_r) #normed=True)
    extent = [ybins[0], ybins[-1], xbins[0], xbins[-1]]

fs = 18 
cbar_ax = f.add_axes([0.92, 0.20, 0.01, 0.70])
im = ax7.imshow(counts.T, cmap=plt.cm.gray_r, extent=extent, norm=LogNorm(vmin=0.01, vmax=1), interpolation = 'Nearest',origin = 'lower')
#im = ax7.imshow(test1.T, cmap=plt.cm.gray, extent=extent, norm = None, interpolation = 'Nearest',origin = 'lower')
test = f.colorbar(im, cax=cbar_ax)
test.set_label("density", fontsize = fs) 
biasage,rmsage = returnscatter(np.log10(age), np.log10(ageout))
biasmass,rmsmass = returnscatter(np.log10(mass), ln_massout/np.log(10.0))
plt.show()

#cbar_ax = f.add_axes([0.92, 0.20, 0.01, 0.65])
#test = colorbar(s1, cax=cbar_ax) 
#test.set_label("SNR", fontsize = 20) 

def mkn_set_ticks(ax, xticks):
    ax.set_xticks(xticks)
    ax.set_yticks(xticks)

ax1.set_xlim(3900,5200)
mkn_set_ticks(ax1, [4000,4500,5000])
ax2.set_xlim(1.2,3.5)
ax3.set_xlim(-1,0.5)
mkn_set_ticks(ax3, [-1.0,-0.5,-0,0.5])
ax4.set_xlim(-0.08,0.35)
mkn_set_ticks(ax4, [-0.1,0,0.1,0.2,0.3])
if plotage:
    ax5.set_xlim(-0.3,0.7)
    ax5.set_title(r"log 10 mass (M$_{\odot}$)", fontsize = fs)
    ax6.set_xlim(-0.5,1.5)
    ax6.set_title(r"log 10 age (Gyr) ", fontsize = fs)

for ax in [ax1, ax2, ax3, ax4, ax5,ax6]:
    ax.set_ylim(ax.get_xlim())
    ax.plot(ax.get_xlim(), ax.get_xlim(), "k-")


biast,rmst = returnscatter(t, tout)
biasg,rmsg = returnscatter(g, gout)
biasfeh,rmsfeh = returnscatter(feh, fehout)
biasalpha,rmsalpha = returnscatter(alpha, alphaout)

f1 = 11
ax1.text(0.45, 0.95,"bias, rms= "+str(round(biast,1))+", "+str(round(rmst,2)), ha='center', va='center', 
              transform=ax1.transAxes, fontsize = f1)
ax2.text(0.45, 0.95,"bias, rms= "+str(round(biasg,2))+", "+str(round(rmsg,2)), ha='center', va='center', 
              transform=ax2.transAxes, fontsize = f1)
ax3.text(0.45, 0.95,"bias, rms= "+str(round(biasfeh,3))+", "+str(round(rmsfeh,2)), ha='center', va='center', 
              transform=ax3.transAxes, fontsize = f1)
ax4.text(0.45, 0.95,"bias, rms= "+str(round(biasalpha,3))+", "+str(round(rmsalpha,2)), ha='center', va='center', 
              transform=ax4.transAxes, fontsize = f1)
ax5.text(0.45, 0.95,"bias, rms= "+str(round(biasmass,3))+", "+str(round(rmsmass,2)), ha='center', va='center', 
              transform=ax5.transAxes, fontsize = f1)
ax6.text(0.45, 0.95,"bias, rms= "+str(round(biasage,3))+", "+str(round(rmsage,2)), ha='center', va='center', 
              transform=ax6.transAxes, fontsize = f1)

fsize = 17
lpad = 7
lpad2 = 7 # -3
ax1.set_title("Teff (K) ", fontsize = 20)
ax2.set_title("log g (dex) ", fontsize = 20)
ax3.set_title("[Fe/H] (dex) ", fontsize = 20)
ax4.set_title(r"[$\alpha$/Fe] (dex) ", fontsize = 20)
#ax6.set_title("delta_nu", fontsize = 20)

ax1.set_ylabel("The Cannon Labels", fontsize = 20, labelpad = 10 ) 
##ax3.set_xlabel("Input Labels", fontsize = 20, labelpad = 10 ) 
f.subplots_adjust(hspace=0.42)
f.subplots_adjust(bottom=0.2)
close(f2) 
show()
draw()
#savefig('6labels_together.png', bbox = 'tight', fmt = "png") 
savefig("/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/validation_1639_6.pdf", bbox = 'tight', fmt = "pdf")
