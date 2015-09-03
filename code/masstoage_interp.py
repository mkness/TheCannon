from scipy import interpolate
import numpy as np

feh_list,mass_list = np.loadtxt('/Users/ness/new_laptop/Apogee_ages/apokasc_ages_actual.txt', usecols = (4,12), unpack =1) 
dir1 = '/Users/ness/new_laptop/Apogee_ages/Marie/'

def makeage(fehval, massval):
    interpflag = [] 
    feh_list=np.arange(-0.9,0.6, 0.2)
    pickmodel1_ind, pickmodel2_ind = np.argsort(np.abs(fehval - feh_list))[0:2] 
    pickmodel1_feh = np.round(feh_list[pickmodel1_ind],1)
    pickmodel2_feh = np.round(feh_list[pickmodel2_ind],1)
    dist1 = np.abs(fehval-feh_list)[pickmodel1_ind]
    dist2 = np.abs(fehval-feh_list)[pickmodel2_ind]
    mass_1,age_1=np.loadtxt(dir1+'age_mass_PARSEC_eta02_'+str(pickmodel1_feh)+'.txt',unpack=True,usecols=(1,4))
    mass_2,age_2=np.loadtxt(dir1+'age_mass_PARSEC_eta02_'+str(pickmodel2_feh)+'.txt',unpack=True,usecols=(1,4))
    if massval < np.min(mass_1) or massval < np.min(mass_2):
        interpflag.append(1)
        fage_1 = np.max(age_1)
        fage_2 = np.max(age_2)
    elif massval > np.max(mass_1) or massval > np.max(mass_2):
        interpflag.append(1)
        fage_1 = np.min(age_1)
        fage_2 = np.min(age_2)
    else: 
        fage1=interpolate.interp1d(mass_1,age_1,kind='linear')
        fage2=interpolate.interp1d(mass_2,age_2,kind='linear')
        interpflag.append(0)
        fage_1 = fage1(np.clip(massval,np.min(mass_1),4)) 
        fage_2 = fage2(np.clip(massval,np.min(mass_1),4)) 
    if (fehval < pickmodel1_feh) and (fehval < pickmodel2_feh):
        print "makeage():", fehval, fage_1
        fage = fage_1
    elif (fehval > pickmodel1_feh) and (fehval > pickmodel2_feh):
        print "makeage():", fehval, fage_1
        fage = fage_1
    else:
        fage = (fage_1*dist2 + fage_2*dist1) / (dist1+dist2) 
    print fehval, massval , fage
    return fage , interpflag


def get_age_from_mass(feh_list,mass_list):
    ages = np.zeros_like(mass_list)
    for jj,(feh,mass) in enumerate(zip(feh_list,mass_list)):
        age, flag = makeage(feh,mass) 
        ages[jj] = age
    return ages

#feh_list,mass_list = loadtxt('/Users/ness/new_laptop/Apogee_ages/apokasc_ages_actual.txt', usecols = (4,12), unpack =1) 
#def makeage(fehval, massval):
