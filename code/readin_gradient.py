import glob
from glob import glob 
from scipy import interpolate 
import pickle
import numpy as np
import matplotlib.pyplot as plt

#To read ".npz"
#file_test = 'gradient_spectra_wrt_logg_t05000g2.50.npz'
file_cannon = 'coeffs_2nd_order_4.pickle'

#dir1 = '/Users/ness/new_laptop/Apogee_ages/YuanSen/Gradient_Spectra_elem/'
#file_test1 = dir1+'at12_feh+0.00_ele14+0.20_t05000g2.50.spec.gz.npz'
#file_test2 = dir1+'at12_feh+0.00_ele12+0.20_t05000g2.50.spec.gz.npz'
#file_test3 = dir1+'at12_feh+0.00_ele20+0.20_t05000g2.50.spec.gz.npz'

a = open(file_cannon,'r')
b = pickle.load(a)
a.close() 
coeffs = b[4]
alpha_coeffs = b[4][:,4]
wl_apogee = b[0][:,0,0]

def makeplot(num):
    dir1 = '/Users/ness/new_laptop/Apogee_ages/YuanSen/Gradient_Spectra_elem/'
    file_test = dir1+'at12_feh+0.00_ele'+str(num)+'+0.20_t05000g2.50.spec.gz.npz'
    data = np.load(file_test) 
    wavelength = data["wavelength"]
    gradient_value = data["gradient_value"]
   # data1 = np.load(file_test1) 
   # data2 = np.load(file_test2) 
   # data3 = np.load(file_test3) 
   # wavelength = data["wavelength"]
   # gradient_value1 = data1["gradient_value"]
   # gradient_value2 = data2["gradient_value"]
   # gradient_value3 = data3["gradient_value"]
   # gradient_value = mean([gradient_value1, gradient_value2, gradient_value3],axis=0) 

    #plt.plot(wavelength, gradient_value) #should look like plots in the overleaf.
    plot(wavelength, gradient_value) 
    title(str(num), fontsize = 20) 
    #f = interpolate.interp1d(wavelength, gradient_value) 
    #gradient_value_interpolated = f(wl_apogee[0:-1])
    #wl_apogee_cut = wl_apogee[0:-1]
    #alpha_coeffs_cut = alpha_coeffs[0:-1]
    #plot(wl_apogee_cut, alpha_coeffs_cut, 'k')
    #plot(wl_apogee_cut, gradient_value_interpolated, 'b')
    return

#makeplot(19) 

def makeplot2(num, wl_apogee, coeffs, save=None):
  dir_in = '/Users/ness/new_laptop/Apogee_ages/YuanSen/'
  fig, axes = plt.subplots(ncols=1, nrows=2,sharex = False,sharey = False)
  ax1 = axes[0] 
  ax2 = axes[1] 
  if num == 0: 
    file_in =dir_in+'gradient_spectra_wrt_teff_t04500g2.50.npz' 
    coeff_plot = coeffs[:,1]
    # Yuan-Sen units hardcoded (range +250K) 
    factor = 100.
  else:
    file_in =dir_in+'gradient_spectra_wrt_logg_t05000g2.50.npz'
    coeff_plot = coeffs[:,2]
    # Yuan-Sen units hardcoded  (range +0.5) 
    factor = 0.3
  data1 = np.load(file_in) 
  wavelength = data1["wavelength"]
  gradient_value1 = data1["gradient_value"]
  ax1.plot(wavelength, gradient_value1 / factor, 'b-',alpha =0.5)  
  ax1.plot(wl_apogee, coeff_plot, 'k-',alpha =0.5) 
  ax2.plot(gradient_value1 / factor, coeff_plot, 'k.', alpha =0.5)
  ax2.set_xlabel("Y-S gradient")
  ax2.set_ylabel("The Cannon gradient")
  ax2.plot([-1,1], [-1,1], 'k-') 
  val1,val2 = 16600, 16900 
  ax1.set_xlim(val1, val2)
  if num == 0:
      ax1.set_title("Teff")
      ax1.set_ylim(-.0002,.0002)
  else:
      ax1.set_title("logg")
      ax1.set_ylim(-0.08, 0.08)
  ax2.set_xlim(ax1.get_ylim())
  ax2.set_ylim(ax1.get_ylim())
  if save is not None:
      plt.savefig(save)

makeplot2(1, wl_apogee, coeffs, save="logg.png") 
makeplot2(0, wl_apogee, coeffs, save="teff.png")

def CR_bounds():
    """
    Cramer-Rao bounds for signal to noise of 100 spectra 
    """
    filenames = glob("./gradient*npz") + glob("./Gradient_Spectra_elem/*npz") 
    for fn in filenames: 
        data = np.load(fn) 
        gradient_value = data["gradient_value"] 
        crb = 1.0/np.sqrt(np.sum(gradient_value**2 * 1.e4) ) 
        print crb, fn, crb 

#CR_bounds() 
