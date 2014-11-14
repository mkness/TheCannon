#!/usr/bin/python 
from glob import glob
import pyfits 
import matplotlib.pyplot as plt
import matplotlib as mpl 
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
rcParams["xtick.labelsize"] = 16
rcParams["ytick.labelsize"] = 16
#font = {'family':'serif','size':16, 'weight':'normal'}
#plt.rc('font',**font)
#plt.rc('legend',**{'fontsize':14})

def makeiso(number):

    filein = glob("/Users/ness/old_laptop/workdir/Apogee/*allStar*304*fits*")[0]
    hdulist = pyfits.open(filein)
    datain = hdulist[1].data
    apstarid= datain['APSTAR_ID'] 
    J= datain['J'] 
    K= datain['K'] 
    fields = datain['FIELD'] 
    at2 = datain['APOGEE_TARGET2'] 
    loc = datain['LOCATION_ID']
    pickfield = fields == '4102'
    Fehall = datain['METALS']
    Fehall_err = datain['METALS_ERR']
    loggall = datain['LOGG']
    loggall_err = datain['LOGG_ERR']
    teffall = datain['TEFF']
    teffall_err = datain['TEFF_ERR']
    vscatter = datain['VSCATTER']
    velall = datain['VHELIO_AVG']
    velall_err = datain['VERR']
    lval = datain['GLON']
    bval = datain['GLAT']
    Gal = 220*cos(bval*pi/180)*sin(lval*pi/180)+16.5*(sin(bval*pi/180)*sin(25*pi/180)+cos(bval*pi/180)*cos(25*pi/180)*cos((lval-53)*pi/180))
    velall = array(velall)
    vgal = velall + Gal #+ 10
    alphaall = datain['ALPHAFE']
    FLAG_1all = datain['PARAMFLAG']
    loc = array(loc)
    if logical_or(logical_or(filein == 4260 ,filein == 4161), logical_or(filein == 4202, filein == 4201)): 
      pickit = logical_and(loc == filein, Fehall < 2.) 
    else:
      pickit = logical_and(loc == filein, Fehall < 200.) 
    appick = apstarid[pickit]
    Fehpick = Fehall[pickit]
    vscatterpick = vscatter[pickit]
    Fehpick_err = Fehall_err[pickit]
    gpick = loggall[pickit]
    gpick_err = loggall_err[pickit]
    Tpick = teffall[pickit]
    Tpick_err = teffall_err[pickit]
    velpick = velall[pickit]
    fehsort = sort(Fehpick)
    argis = argsort(Fehpick) 
    lval = array(lval)
    bval = array(bval)
    lvalpick = lval[pickit]
    bvalpick = bval[pickit]
    at2pickit= array(at2)[pickit]
    JmK = array(J) - array(K) 
    JmK_pickit = JmK[pickit]
    JmK_pickit = array(JmK_pickit) 
    t,t_err,g,g_err,feh,feh_err = Tpick,Tpick_err,gpick,gpick_err,Fehpick, Fehpick_err
    t2_data,g2_data,feh2_data,vel2,l2,b2,starloc2,targ2,J2,K2,SNR2,STARFLAG2,tA,gA,fehA,vscatter,chi2 = loadtxt('play_v19.txt', usecols = (0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17), unpack =1)
    rotation_warn = loadtxt('starflags_v19.txt', usecols = (0,), unpack =1) 
    t2_data,g2_data,feh2_data,vel2,l2,b2,starloc2,targ2,J2,K2,SNR2,STARFLAG2,tA,gA,fehA,vscatter,chi2 = loadtxt('play_v20mkn.txt', usecols = (0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17), unpack =1)
    rotation_warn = loadtxt('starflags_v20mkn.txt', usecols = (0,), unpack =1) 
    fehA = array(fehA) 
     

    rcParams['figure.figsize'] = 12.0, 12.0
    fig, axs = plt.subplots(2,2, figsize=(12, 10), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .001, wspace=.001)
    axs = axs.ravel()
    if number ==1 :# all 
      pick1 =  array([ np.int(a) & 2**1 for a in targ2]) == 0 
      pick2 =  array([ np.int(a) & 2**2 for a in targ2]) == 0 
      pick3 =  array([ np.int(a) & 2**9 for a in targ2]) == 0 
      pick4 =  logical_and(vscatter < 10.2 , rotation_warn == 0) 
      pick4a =  vscatter > 10.0 
      pick =  array([ np.int(a) & 2**9 for a in targ2]) == 0 
      pick = logical_and(logical_and( logical_and(pick1, pick2) , logical_and(pick4, pick3) ) , chi2 < 3000) 
      pick_other = logical_and( logical_and(pick1, pick2) , logical_and(pick4a, pick3) ) 
    if number ==2 : #aspcap cross over
      picka = logical_and(fehA < 1.4, fehA > -3) 
      pickb = array([ np.int(a) * 2**9 for a in targ2]) == 0 
      pick = logical_and(picka, pickb) 
    if number ==3 :# not aspcap 
      picka = logical_or(fehA > 1.4, fehA < -3) 
      pickb = array([ np.int(a) & 2**9 for a in targ2]) == 0 
      pick = logical_and(picka, pickb) 
    if number ==4 :# in no-man's land  
      picka = logical_and(fehA < 1.4, fehA > -3) 
      m = 460
      c = 3200 
      pick_gt =logical_and( g2_data < 4.1,  t2_data < m*g2_data + c ) 
      pickb = array([ np.int(a) * 2**9 for a in targ2]) == 0 
      pick = logical_and(pick_gt, logical_and(picka, pickb) ) 
    cval = ['k', 'b', 'r'] 

    ##
    file_pt3 = '10Gyr_0.030.txt'   
    file_pt2 = '10Gyr_019.txt'   
    file_0 = '10Gyr_0.0155.txt'
    file_mpt9 = '10Gyr_0.0019.txt' 
    file_mpt5 = '10Gyr_0.0048.txt' 
    file_m1 = '10Gyr_0.0015.txt' 
    file_pt1 = '10Gyr_019.txt' 
    file_m2 = '10Gyr_0.00017.txt' 
    file_m2pt2 = '10Gyr_0_0.0001.txt' 
    file_m1pt5 =  '10Gyr_0.00045.txt'
    dirin = '/Users/ness/old_laptop/workdir/Apogee/isochrones/' 
    t0,g0 = loadtxt(dirin+file_pt3, usecols = (5,6), unpack =1)#,label = '+0.2') 
    t1,g1 = loadtxt(dirin+file_pt2, usecols = (5,6), unpack =1)#,label = '+0.2') 
    t2,g2 = loadtxt(dirin+file_0,   usecols = (5,6), unpack =1)#, label = '+0.0') 
    t3,g3 = loadtxt(dirin+file_mpt5,usecols = (5,6), unpack =1)#, label = '-0.5') 
    t4,g4 = loadtxt(dirin+file_m1,  usecols = (5,6), unpack =1)#, label = '-1.0') 
    t5,g5 = loadtxt(dirin+file_m1pt5,  usecols = (5,6), unpack =1)#, label = '-2.0') 
    ts = [t2,t3,t4,t5]
    gs = [g2,g3,g4,g5]
    cl = [ '#FF0000', '#FFCC00', '#7FFFD4' , 'blue' ]
    label1 = [ '[Fe/H] = +0.0', '[Fe/H] = -0.5', '[Fe/H] = -1.0', '[Fe/H] =  -1.5' ] 
    label2 = [ '[Fe/H] $>$ -0.25 ', '-0.75 $<$ [Fe/H] $<$ -0.25 ', '-1.25 $<$ [Fe/H] $<$ -0.75 ', '[Fe/H] $<$ -1.25 ' ] 
    cl = [0,-0.5,-1,-2] 
    mycmap = cm.jet 
    fehinds = [-3.75, -1.25, -0.75, -0.25, 2.25]    
    fehinds = [2.25, -0.25, -0.75, -1.25, -3.75]
    label21 = ['[Fe/H] < -1.25', '-1.25 < [Fe/H] < -0.75','-0.75 < [Fe/H] < -0.25', '[Fe/H] > -0.25']
    for i in range(0,4):
        pickfeh = logical_and(feh2_data > fehinds[i+1], feh2_data < fehinds[i])
        pickboth = logical_and(pickfeh, pick) 
        print len(pickboth[pickboth]) 
        s1 = axs[i].plot(10**ts[i],gs[i], color = 'grey', label = label1[i] ,alpha  = 1.0,linewidth =1) 
        axs[i].legend(loc = 2,fontsize = 14)#,frameon = False) 
        #axs[i].plot(10**ts[i],gs[i], color = 'k', alpha  = 1.0,linewidth =0.5,linestyle= 'dashed') 
        pickboth = logical_and(pick, pickfeh) 
        #s2 = axs[i].plot(t2_data[pickboth], g2_data[pickboth], 'ko', ms = 2,alpha = 0.05)
        image =axs[i].hist2d(t2_data[pickboth], g2_data[pickboth],bins = 40,cmin =0, cmap = cm.gray_r)#, 'ko', ms = 2,alpha = 0.05)
        #bad = isnan(image[0])
        #image[0][bad] = 0
        #data1 = [] 
        #data2 = [] 
        #for i in range(0,len(image[0])):
        #  data1.append( mean(image[1][i]+image[1][i+1])*0.5)
        #  data2.append( mean(image[2][i]+image[2][i+1])*0.5)
        #test1 = ndimage.gaussian_filter(image[0]/image[0], order = 0, sigma = 1 ) 
        #data1 = array(data1) 
        #data2 = array(data2) 
        #test = image[0] 
        #test2 = image[0] == 0.
        #cm.jet.set_bad(color='white', alpha=None)
        #masked_data = ma.masked_array(test,test2)
        ##axs[i].pcolormesh(data2,data1,masked_data, cmap = cm.gray_r , alpha =0.6) 
        #axs[i].contour(data2, data1, masked_data, 40,linewidths = 1.5, cmap = cm.gray_r,alpha = 1.0,linestyle = 'dashed') 
        #axs[0].pcolormesh(data2,data1,masked_data, cmap = cm.gray , alpha =0.6) 
        #axs[i].contour(data2, data1, masked_data, 40,linewidths = 1.5, cmap = cm.jet,alpha = 1.0,linestyle = 'dashed') 
        #leg = round(mean(feh2_data[pickboth]) , 1) 
        axs[i].text(0.7, 0.2,str(label2[i]), ha='center', va='center',transform=axs[i].transAxes,fontsize = 18)
        axs[i].set_xlim(6000, 3100) 
        axs[i].set_ylim(5.5, -1) 
    
    print len(t2_data)
    print len(t2_data[pick]) 
   
    axs[0].set_xticklabels([])
    axs[1].set_xticklabels([])
    axs[3].set_yticklabels([])
    axs[1].set_yticklabels([])
    #axs[3].set_yticklabels([])
    #axs[i].set_xlabel(r'\textsl{The Cannon}  \mbox{$T_{eff}$}, (K)', fontsize = 20,labelpad = 10) 
    #axs[2].set_xlabel('\mbox{$T_{eff}$}, (K)', fontsize = 20,labelpad = 10) 
    #axs[3].set_xlabel('\mbox{$T_{eff}$}, (K)', fontsize = 20,labelpad = 10) 
    #axs[0].set_ylabel(r'$\log g$, (dex)', fontsize = 20,labelpad = 5) 
    #axs[2].set_ylabel(r'$\log g$, (dex)', fontsize = 20,labelpad = 5) 
    axs[2].set_xlabel('\mbox{$T_{eff}$}, (K)', fontsize = 24,labelpad = 10) 
    axs[3].set_xlabel('\mbox{$T_{eff}$}, (K)', fontsize = 24,labelpad = 10) 
    axs[0].set_ylabel('log g, (dex)', fontsize = 24,labelpad = 5) 
    axs[2].set_ylabel('log g, (dex)', fontsize = 24,labelpad = 5) 
            
    fig.subplots_adjust(hspace=0.0)
    if number ==1: 
      prefix = str("DR10_isochrone_v18")
    if number ==2: 
      prefix = str("xover_DR10_isochrone")
    if number ==3: 
      prefix = str("no_ASPCAP_params_DR10_isochrone")
    if number ==4: 
      prefix = str("Cannon_triangle_cut_DR10_isochrone")
    savefig3(fig, prefix, transparent=False, bbox_inches='tight', pad_inches=0.5)
    savefig3(fig, prefix, transparent=False, bbox_inches='tight', pad_inches=0.5)
    return 

def savefig3(fig, prefix, **kwargs):
 #   for suffix in (".png"):
    suffix = ".png"
    print "writing %s" % (prefix + suffix)
    fig.savefig(prefix + suffix)#, **kwargs)

makeiso(1) 

