from scipy import interpolate
import numpy as np

feh_list,mass_list = loadtxt('/Users/ness/new_laptop/Apogee_ages/apokasc_ages_actual.txt', usecols = (4,12), unpack =1) 

def makeage(fehval, massval):
    interpflag = [] 
    feh_list=np.arange(-0.9,0.6, 0.2)
    pickmodel1_ind = argsort(abs(fehval - feh_list))[0] 
    pickmodel2_ind = argsort(abs(fehval - feh_list))[1]
    pickmodel1_feh = round(feh_list[pickmodel1_ind],1)
    pickmodel2_feh = round(feh_list[pickmodel2_ind],1)
    dist1 = sort(abs(fehval-feh_list))[0]
    dist2 = sort(abs(fehval-feh_list))[1]
    mass_1,age_1=np.loadtxt('age_mass_PARSEC_eta02_'+str(pickmodel1_feh)+'.txt',unpack=True,usecols=(1,4))
    mass_2,age_2=np.loadtxt('age_mass_PARSEC_eta02_'+str(pickmodel2_feh)+'.txt',unpack=True,usecols=(1,4))
    fage1=interpolate.interp1d(mass_1,age_1,kind='linear')
    fage2=interpolate.interp1d(mass_2,age_2,kind='linear')
    if massval < min(mass_1):
      interpflag.append(1)
    else: 
      interpflag.append(0)
    fage_1 = fage1(clip(massval,min(mass_1),4)) 
    fage_2 = fage1(clip(massval,min(mass_1),4)) 
    fage = (fage_1*dist2 + fage_2*dist1) / (dist1+dist2) 
    print fehval, massval 
    return fage , interpflag

ages = []
flags = []
for a,b in zip(feh_list,mass_list):
    age, flag = makeage(a,b) 
    ages.append(age)
    flags.append(flag)
#feh_list,mass_list = loadtxt('/Users/ness/new_laptop/Apogee_ages/apokasc_ages_actual.txt', usecols = (4,12), unpack =1) 
#def makeage(fehval, massval):
