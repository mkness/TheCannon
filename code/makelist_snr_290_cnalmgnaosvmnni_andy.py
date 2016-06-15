import pyfits 
filein = '/Users/ness/new_laptop/Apogee_DR12_play/allStar-v603.fits'
hdulist = pyfits.open(filein)
datain = hdulist[1].data
apstarid= datain['APSTAR_ID'] 
targetid= datain['TARGET_ID'] 
vscatter= datain['VSCATTER'] 
aspcap_id= datain['ASPCAP_ID'] 
fields = datain['FIELD'] 
loc = datain['LOCATION_ID']
params = datain['PARAM']
t = params[:,0]
g = params[:,1]
#M = params[:,3]
feh = datain['FE_H']
SNR = datain["SNR"]
c = datain['C_H']
n = datain['N_H']
cfe  = c - feh
nfe =  n - feh
cm = datain['ELEM'][:,0]
nm = datain['ELEM'][:,1]
om = datain['ELEM'][:,2]
nah = datain['ELEM'][:,3]
mgm = datain['ELEM'][:,4]
alh= datain['ELEM'][:,5]
sm = datain['ELEM'][:,7]
vh = datain['ELEM'][:,11]
mnh = datain['ELEM'][:,12]
nih = datain['ELEM'][:,14]
pickit1 = logical_and(abs(cm) < 2, abs(nm) < 2)
pickit2 = logical_and(abs(om) < 2, abs(nah) < 2)
pickit3 = logical_and(abs(mgm) < 2, abs(alh) < 2)
pickit4 = logical_and(abs(sm) < 2, abs(vh) < 2)
pickit5 = logical_and(abs(nih) < 2, abs(mnh) < 2)
pickitall1 = logical_and(logical_and(pickit1,pickit2), pickit3) 
takeelem = logical_and(logical_and(pickit4,pickit5), pickitall1) 
aspcapflag = datain["ASPCAPFLAG"]
starflag = datain["STARFLAGS"]
alpha = params[:,-1]
alphatest = alpha
fehtest = feh 
rot = []
bad = []

for each in aspcapflag:
  rot.append(each & 2**10) 
  bad.append(each & 2**23) 
apid = []
for each in apstarid:
  apid.append( each.split('.')[-1][2:]) 
bad = array(bad) 
apid = array(apid) 
take1 = logical_and(alphatest > -0.1,logical_and(abs(fehtest) < 3,logical_and(SNR > 285, SNR < 305)))
take1 = logical_and(alphatest > -0.4,logical_and(abs(fehtest) < 3,logical_and(SNR > 285, SNR < 305)))
take2 = logical_and(vscatter <1,logical_and(aspcapflag ==  0, vscatter  < 1))
take3 = logical_and(bad == 0, loc != 1) 
ind = logical_and(logical_and(takeelem, take3), logical_and(take1,take2)) 

loctake = loc[ind]
idtake = apid[ind]  
t_take = t[ind]
g_take = g[ind]
feh_take = feh[ind]
alpha_take = alpha[ind]
c_take = cm[ind]
n_take = nm[ind]
o_take = om[ind]
na_take = nah[ind]
mg_take = mgm[ind]
al_take = alh[ind]
s_take = sm[ind]
v_take = vh[ind]
mn_take = mnh[ind]
ni_take = nih[ind]
snr_take = SNR[ind]
vscatter_take = vscatter[ind]
rot = array(rot) 
rot_take = rot[ind]
aspcap_take = aspcapflag[ind]
starflag_take = starflag[ind]

text2 ='/home/ness/new_laptop/Apogee_DR12/data.sdss3.org/sas/dr12/apogee/spectro/redux/r5/stars/l25_6d/v603/'
text3 = '/aspcapStar-r5-v603-2M' 
text4 = '.fits'
listin = [] 
for a,b in zip(loctake, idtake):
  listin.append(text2+str(a)+text3+str(b)+text4) 
data = zip(listin, t_take, g_take, feh_take, alpha_take,c_take, n_take, al_take, mg_take, na_take, o_take, s_take, v_take, mn_take, ni_take, snr_take,vscatter_take,rot_take,aspcap_take,starflag_take) 
savetxt('training_SNR290_snr_cnalmgnaosvmni_andy.list', data, fmt = "%s" ) 
print len(data) 
