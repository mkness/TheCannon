#!/usr/bin/python 
import scipy 
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
rcParams['figure.figsize'] = 23.0, 4.0
f, (ax1, ax2, ax3,ax4,ax5) = plt.subplots(ncols=5)


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
var = alphaout 
var = snr
vmin1,vmax1 = 80,300
#vmin1,vmax1  = 4050, 5100
#var =tout
cval = cm.cool
msize = 20.
c1 = ax1.scatter(t,tout, c = var,linewidth =0, s = msize ,vmin = vmin1, vmax = vmax1 , cmap = cval)
ax2.scatter(g,gout, c = var,linewidth =0, s = msize ,vmin = vmin1, vmax = vmax1 , cmap = cval)
ax3.scatter(feh,fehout, c = var,linewidth =0, s = msize ,vmin = vmin1, vmax = vmax1 , cmap = cval)
ax4.scatter(alpha,alphaout, c = var,linewidth =0, s = msize ,vmin = vmin1, vmax = vmax1 , cmap = cval)
#s1 = ax5.scatter(mass,np.exp(ln_massout), c = var, s = msize ,linewidth = 0 ,vmin = vmin1, vmax = vmax1, cmap = cval)
s1 = ax5.scatter(np.log10(mass),ln_massout/np.log(10.), c = var, s = msize ,linewidth = 0 ,vmin = vmin1, vmax = vmax1, cmap = cval)
cbar_ax = f.add_axes([0.92, 0.20, 0.01, 0.65])
test = colorbar(s1, cax=cbar_ax) 
#test.set_label("[alpha/Fe] ", fontsize = 20) 
test.set_label("SNR", fontsize = 20) 

ax1.set_xlim(3900,5200)
ax1.set_xticks([4000,4500,5000])
ax1.plot([3200,7000], [3200,7000], 'k') 
ax2.set_xlim(1.2,3.5)
ax2.plot([0,5],[0,5])
ax3.set_xlim(-1,0.5)
ax3.plot([-2.1,1],[-2.1,1])
ax3.set_xticks([-1.0,-0.5,-0,0.5])
ax4.set_xlim(-0.08,0.35)
ax4.plot([-0.1,1], [-0.1,1])
ax4.set_xticks([-0.1,0,0.1,0.2,0.3])
ax5.set_xlim(-0.3,0.7)
ax5.plot([-2,4], [-2,4], 'k') 
for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.set_ylim(ax.get_xlim())
    ax.set_yticks(ax.get_xticks())

biast,rmst = returnscatter(t, tout)
biasg,rmsg = returnscatter(g, gout)
biasfeh,rmsfeh = returnscatter(feh, fehout)
biasalpha,rmsalpha = returnscatter(alpha, alphaout)
biasnu1,rmsnu1 = returnscatter(np.log10(mass), ln_massout/np.log(10.0))

f1 = 12
ax1.text(0.45, 0.95,"bias, rms= "+str(round(biast,1))+", "+str(round(rmst,2)), ha='center', va='center', 
              transform=ax1.transAxes, fontsize = f1)
ax2.text(0.45, 0.95,"bias, rms= "+str(round(biasg,2))+", "+str(round(rmsg,2)), ha='center', va='center', 
              transform=ax2.transAxes, fontsize = f1)
ax3.text(0.45, 0.95,"bias, rms= "+str(round(biasfeh,3))+", "+str(round(rmsfeh,2)), ha='center', va='center', 
              transform=ax3.transAxes, fontsize = f1)
ax4.text(0.45, 0.95,"bias, rms= "+str(round(biasalpha,3))+", "+str(round(rmsalpha,2)), ha='center', va='center', 
              transform=ax4.transAxes, fontsize = f1)
ax5.text(0.45, 0.95,"bias, rms= "+str(round(biasnu1,3))+", "+str(round(rmsnu1,2)), ha='center', va='center', 
              transform=ax5.transAxes, fontsize = f1)

fsize = 17
lpad = 7
lpad2 = 7 # -3
ax1.set_title("Teff (K) ", fontsize = 20)
ax2.set_title("log g (dex) ", fontsize = 20)
ax3.set_title("[Fe/H] (dex) ", fontsize = 20)
ax4.set_title(r"[$\alpha$/Fe] (dex) ", fontsize = 20)
ax5.set_title(r"log_10 mass", fontsize = 20)
#ax6.set_title("delta_nu", fontsize = 20)

ax1.set_ylabel("THE CANNON LABELS", fontsize = 20, labelpad = 10 ) 
ax3.set_xlabel("INPUT LABELS", fontsize = 20, labelpad = 10 ) 
f.subplots_adjust(hspace=0.42)
f.subplots_adjust(bottom=0.2)
savefig('6labels.png', bbox = 'tight', fmt = "png") 
show()
