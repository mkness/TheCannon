#!/usr/bin/python 
import scipy 
import matplotlib 
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
rc('text', usetex=True)
rc('font', family='serif')
rcParams["xtick.labelsize"] = 13
rcParams["ytick.labelsize"] = 13
rcParams['figure.figsize'] = 24.0, 8.0

f, ((ax1, ax2, ax3,ax4,ax5,ax6),(ax1a,ax2a,ax3a,ax4a,ax5a,ax6a)) = plt.subplots(ncols=6,nrows=2,sharex=False,sharey=False)
#ax = f.add_subplot(111,frameon = 0 )
#ax.set_xlabel("Input Labels",labelpad = 5, fontsize = 20 )
#ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
f2, (ax7,ax8) = plt.subplots(ncols=2)

# set labels 
fs = 20 
ax1.set_xlabel("Teff$_{\mbox{ASPCAP}}$", fontsize = fs) 
ax2.set_xlabel("logg$_{\mbox{KEPLER}}$", fontsize = fs) 
ax3.set_xlabel("[Fe/H]$_{\mbox{ASPCAP}}$", fontsize = fs) 
ax4.set_xlabel(r"[$\alpha$/Fe]$_{\mbox{ASPCAP}}$", fontsize = fs) 
ax5.set_xlabel("log$_{10}$ mass$_{\mbox{KEPLER}}$", fontsize = fs) 
ax6.set_xlabel("log$_{10}$ age$_{\mbox{KEPLER}}$", fontsize = fs) 
ax1a.set_xlabel(r"$\Delta_{\mbox{Teff}}$", fontsize = fs) 
ax2a.set_xlabel(r"$\Delta_{\mbox{logg}}$", fontsize = fs) 
ax3a.set_xlabel(r"$\Delta_{\mbox{[Fe/H]}}$", fontsize = fs) 
ax4a.set_xlabel(r"$\Delta_{\mbox{$\alpha$/Fe}}$", fontsize = fs) 
ax5a.set_xlabel(r"$\Delta_{\mbox{log$_{10}$ mass}}$", fontsize = fs) 
ax6a.set_xlabel(r"$\Delta_{\mbox{log$_{10}$ age}}$", fontsize = fs) 
ax4a.set_xticks([-0.1,0,0.1]) 
ax2a.set_xticks([-0.5,0,0.5]) 


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
dir2 = '/Users/ness/new_laptop/Apogee_ages/crossvalidation_quad/' 
starsout_files = glob(dir2+"random_not*log*")
training_files =  glob(dir2+"training*5_*_log.pickle*")
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

axall_bottom = [ax1a,ax2a,ax3a,ax4a,ax5a,ax6a]
for a,b,ax in zip(varin,varout,axall_bottom):
    ax.hist(a-b,bins=30,histtype = 'step', color = 'k',normed=0)
    lims = max(abs(a-b))*1.05
    ax.set_xlim(-1*lims,lims) 

#cbar_ax = f.add_axes([0.94, 0.20, 0.01, 0.70])
cbar_ax = f.add_axes([0.94, 0.595, 0.01, 0.30])
im = ax7.imshow(counts.T, cmap=plt.cm.gray_r, extent=extent, norm=LogNorm(vmin=0.01, vmax=1), interpolation = 'Nearest',origin = 'lower')
#im = ax7.imshow(test1.T, cmap=plt.cm.gray, extent=extent, norm = None, interpolation = 'Nearest',origin = 'lower')
test = f.colorbar(im, cax=cbar_ax)
test.set_label("Number of Stars", fontsize = fs) 
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
    ax5.set_title("log$_{10}$ mass (M$_{\odot}$)", fontsize = fsize)
    ax6.set_xlim(-0.5,1.5)
    ax6.set_title("Derived log$_{10}$ age (Gyr) ", fontsize = fsize)


for ax in [ax1, ax2, ax3, ax4, ax5,ax6]:
    ax.set_ylim(ax.get_xlim())
    ax.plot(ax.get_xlim(), ax.get_xlim(), "k-")


biast,rmst = returnscatter(t, tout)
biasg,rmsg = returnscatter(g, gout)
biasfeh,rmsfeh = returnscatter(feh, fehout)
biasalpha,rmsalpha = returnscatter(alpha, alphaout)

f1 = 12
#axall_bottom = [ax1a,ax2a,ax3a,ax4a,ax5a,ax6a]
biases = [str(round(biast,1)),str(round(biasg,2)),str(round(biasfeh,3)),str(round(biasalpha,3)),str(round(biasmass,2)), str(round(biasage,2))]
rmses = [str(round(rmst,1)),str(round(rmsg,2)),str(round(rmsfeh,2)),str(round(rmsalpha,2)),str(round(rmsmass,2)), str(round(rmsage,2))]
for ax,biasval,rmsval in zip(axall_bottom, biases, rmses):
    biasval = 'bias = '+biasval
    rmsval = 'rms = '+rmsval
    ax.annotate(biasval, xy=(0, 0), xycoords='axes fraction', fontsize=f1,
                xytext=(1, 175), textcoords='offset points',
                ha='left', va='top')
    ax.annotate(rmsval, xy=(0, 0), xycoords='axes fraction', fontsize=f1,
                xytext=(1, 164), textcoords='offset points',
                ha='left', va='top')
  
#ax1.text(0.45, 0.95,"bias, rms= "+str(round(biast,1))+", "+str(round(rmst,2)), ha='center', va='center', 
#              transform=ax1.transAxes, fontsize = f1)
#ax2.text(0.45, 0.95,"bias, rms= "+str(round(biasg,2))+", "+str(round(rmsg,2)), ha='center', va='center', 
#              transform=ax2.transAxes, fontsize = f1)
#ax3.text(0.45, 0.95,"bias, rms= "+str(round(biasfeh,3))+", "+str(round(rmsfeh,2)), ha='center', va='center', 
#              transform=ax3.transAxes, fontsize = f1)
#ax4.text(0.45, 0.95,"bias, rms= "+str(round(biasalpha,3))+", "+str(round(rmsalpha,2)), ha='center', va='center', 
#              transform=ax4.transAxes, fontsize = f1)
#ax5.text(0.45, 0.95,"bias, rms= "+str(round(biasmass,3))+", "+str(round(rmsmass,2)), ha='center', va='center', 
#              transform=ax5.transAxes, fontsize = f1)
#ax6.text(0.45, 0.95,"bias, rms= "+str(round(biasage,3))+", "+str(round(rmsage,2)), ha='center', va='center', 
#              transform=ax6.transAxes, fontsize = f1)

fsize = 20
lpad = 7
lpad2 = 7 # -3
ax1.set_title("Teff (K) ", fontsize = fsize)
ax2.set_title("log g (dex) ", fontsize = fsize)
ax3.set_title("[Fe/H] (dex) ", fontsize = fsize)
ax4.set_title(r"[$\alpha$/Fe] (dex) ", fontsize = fsize)
#ax6.set_title("delta_nu", fontsize = 20)

ax1.set_ylabel("The Cannon Labels", fontsize = fsize, labelpad = 10 ) 
ax1a.set_ylabel("Number of Stars", fontsize = fsize, labelpad = 10 ) 
##ax3.set_xlabel("Input Labels", fontsize = 20, labelpad = 10 ) 
f.subplots_adjust(hspace=0.3)
f.subplots_adjust(bottom=0.2)
print ax6.get_position() 
addval = 0.025
test = ax6.get_position()
testa = ax6a.get_position()
test = array(test) 
testa = array(testa) 
print test, testa
xa,xb,ya,yb = test[0][0],test[0][1], test[1][0], test[1][1]
xaa,xba,yaa,yba = testa[0][0],testa[0][1], testa[1][0], testa[1][1]
xa,xb,ya,yb = xa+addval,xb, ya+addval, yb 
xaa,xba,yaa,yba = xaa+addval,xba, yaa+addval, yba 
ax6.set_position([xa,xb,ya-xa,yb-xb]) 
ax6a.set_position([xaa,xba,yaa-xa,yba-xba]) 
#xa,xb,ya,yb = 0.7892851+addval, 0.2, 0.9+0.005, 0.9 
draw() 
close(f2) 
show()
draw()
#savefig('6labels_together.png', bbox = 'tight', fmt = "png") 
savefig("/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/validation_1639_6.pdf", bbox = 'tight', fmt = "pdf")
