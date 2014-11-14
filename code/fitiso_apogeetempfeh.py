a = open('test14_edited_own.txt', 'r')
ta,ga,feha =loadtxt('test14_edited_own.txt', usecols = (3,5,7), unpack =1) 
tc,gc,fehc =loadtxt('test14_edited_own.txt', usecols = (4,6,8), unpack =1) 
#a = open('test14.txt', 'r')
al = a.readlines() 
name = []
bl = []
for each in al:
    name.append(each.split()[1])
    bl.append(each.split()[0])

bl =array(bl) 
unames = unique(name)
name = array(name)
take = [] 
for each in unames: 
  take.append(name == each ) 

J,K = loadtxt("all_2mass_download_edit_cut.txt", usecols = (11, 19), unpack =1) 
ebmv = loadtxt("all_extinction_edit_cut.txt", usecols = (3,), unpack =1) 
J = array(J)
K = array(K)
JmK = J-K 
ebmv = array(ebmv) 
ejmk = 0.535*ebmv 
JmKo = J-K-ejmk

dir1= '/Users/ness/old_laptop/workdir/Apogee/isochrones/'
#f_M107 ,f_M13, f_M15, f_M2, f_M3, f_M5, f_M53, f_M71, f_M92, f_N2158, f_N2420, f_N4147, f_N5446, f_M45, f_M67, f_N188, f_N6791, f_N6819, f_N7789 = \
#        [
#names1 = ["f_M92","f_M15", "f_M53", "f_N5466", "f_N4147", "f_M2", "f_M13", "f_M3", "f_M5", "f_M107", "f_M71", "f_N2158", "f_M35", "f_N2420", "f_N188", "f_M67", "f_N7789", "f_M45","f_N6819",
#"f_N6791" ] 
#metals1 =  [ -2.35, -2.33,-2.06, -1.98, -1.78, -1.66, -1.58, -1.50, -1.33, -1.03, -0.82, -0.28, -0.21, -0.35, -0.03, -0.01, 0.02, 0.03, 0.09, 0.47] 
names1 = ["f_M92","f_M15", "f_M53", "f_N5466", "f_N4147", "f_M2", "f_M13", "f_M3", "f_M5", "f_M107", "f_M71", "f_N2158", "f_M35", "f_N2420", "f_N188", "f_M67", "f_N7789", "f_N6819",
"f_N6791" ] 
metals1 =  [ -2.35, -2.33,-2.06, -1.98, -1.78, -1.66, -1.58, -1.50, -1.33, -1.03, -0.82, -0.28, -0.21, -0.35, -0.03, -0.01, 0.02, 0.09, 0.47] 

feh_dict = {}
for i in range(len(names1)):
    feh_dict[names1[i]] = metals1[i]

M107 = 'M107_edit.txt'
M13 = 'M13.txt'
M15 = 'M15.txt'
M2 = 'M2.txt'
M3 = 'M3.txt'
M5 = 'M5.txt'
M53 = 'M53.txt'
M71 = 'M71.txt'
M92 = 'M92.txt'
N2158 = 'N2158.txt'
N2420 = 'N2420.txt'
N4147 = 'N4147.txt'
N5466 = 'N5466.txt'
#Pleiades = 'Pleiades.txt'
M67 = 'M67.txt'
N118 = 'N118.txt'
N6791 = 'N6791.txt'
N6819 = 'N6819.txt'
N7789 = 'N7789.txt'
K = JmKo
def theta(b0,b1,b2,b3,b4,b5, feh, JmKo):
    result = b0+b1*JmKo + b2*JmKo*JmKo + b3*JmKo*feh + b4*feh + b5*feh*feh
    return 5040.0/result
b0,b1,b2,b3,b4,b5,Na,Nb = 0.6517,0.6312,0.0168, -0.0381, 0.0256,0.0013,145, 94
nameall = [] 
tall = [] 
fehall = [] 
gall = [] 
gaall =[]
taall = []
fehaall = [] 
blall = [] 
for each in unames:
    fehinput = "f_"+each
    isoinput = each+".txt"
    fehval = feh_dict[(fehinput)]
    takeit = array(name) == each 
    JmKotake = JmKo[takeit]
    bltake = bl[takeit]
    fehatake = feha[takeit]
    fehval = fehc[takeit]
    tatake = ta[takeit]
    gatake = ga[takeit]
    temperatures = tc[takeit]
    isoreadin = dir1+isoinput
    logt, g_iso  = genfromtxt(isoreadin, usecols = (5,6), unpack =1) 
    t_iso = 10**logt
    take1 = list(g_iso).index(min(g_iso))
    #take1 = diff(g_iso) > 0
    #nums = arange(0,len(g_iso))[take1][0]
    t_iso = t_iso[0:take1]
    g_iso = g_iso[0:take1]
    #temperatures = theta(b0,b1,b2,b3,b4,b5,fehval , JmKotake) 
    #temp_take = temperatures > min(t_iso) 
    #temperatures2 = temperatures[temp_take] 
    y_new_all = [] 
    for t_each in temperatures: 
        sel_t = logical_and(t_iso > t_each -400, t_iso < t_each + 400 ) 
        sel = logical_and(g_iso < 4, sel_t) 
        t_pick = t_iso[sel]
        g_pick = g_iso[sel]
        t_new = arange(min(t_pick), max(t_pick), 1) 
        g_new = arange(min(g_pick), max(g_pick), 0.01) 
        fa = interpolate.interp1d(sort(t_pick), g_pick[argsort(t_pick)])
        if t_each > min(t_pick):
          new_ydata = fa(t_each) 
        if t_each <= min(t_pick):
          new_ydata = 9999.9
        y_new_all.append(new_ydata) 
        #fehval2 = [fehval]*len(y_new_all)
        fehval2 = fehval
        each2 = [each]*len(y_new_all) 
    nameall.append(each2)
    tall.append(temperatures)
    gall.append(y_new_all)
    fehall.append(fehval2)
    taall.append(tatake)
    gaall.append(gatake)
    fehaall.append(fehatake)
    blall.append(bltake)

blall = hstack((blall))
tall = hstack((tall))
gall = hstack((gall))
fehall = hstack((fehall))
taall = hstack((taall))
gaall = hstack((gaall))
fehaall = hstack((fehaall))
nameall = hstack((nameall))
g_new = [round(a, 2) for a in gall] 
t_new = [round(a, 2) for a in tall] 
g_new =array(g_new) 
t_new =array(t_new) 
arangeit = argsort(blall)
blall =array(blall) 
data = zip(blall[arangeit], nameall[arangeit], taall[arangeit], t_new[arangeit], gaall[arangeit], g_new[arangeit], fehaall[arangeit], fehall[arangeit]) 
savetxt("mkn_labels_Atempfeh.txt", data, fmt = "%s" ) 
al = array(al) 
