params_labels = [params[:,0], params[:,1], params[:,2] , covs_params[:,0,0]**0.5, covs_params[:,1,1]**0.5, covs_params[:,2,2]**0.5 ]  
testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
file2 = '4332_data_all_more.txt'
file2in = testdir+file2
t,t_err,g,g_err,feh,feh_err = loadtxt(file2in, usecols = (1,2,3,4,5,6), unpack =1) 
rcParams['figure.figsize'] = 10.0, 8.0

pick = logical_and(input_ASPCAP[1] > 0, logical_and(input_ASPCAP[3] < 300, input_ASPCAP[2] > -4.0) ) 
x = input_ASPCAP[1][pick]
y = params_labels[1][pick]
yerr = params_labels[2+2][pick]
xerr = input_ASPCAP[2+2][pick]
c_val = input_ASPCAP[0][pick]

starsin = open('starsin_new_all_ordered.txt','r') 
starsin = open('test4_selfg.txt','r') 
al = starsin.readlines()
names = []
for each in al:
  names.append(each.split()[1]) 
unames = unique(names) 
starind = arange(0,len(names), 1) 
name_ind = [] 
name_ind2 = [] 
names = array(names) 
#for each in unames:
#  takeit = each == names 
#  name_ind.append(starind[takeit][-1]+1. ) 
#  name_ind2.append(starind[takeit][0] ) 
#names = array(names) 
for each in unames:
  takeit = each == names 
  name_ind.append(starind[takeit][-1]+1. ) 
  name_ind2.append(starind[takeit][0]+1. ) 
name_ind_sort = sort(name_ind)
name_ind2_sort = sort(name_ind2)

#zip(sort(name_ind), unames[argsort(name_ind)]) 

#T_est,g_est,feh_est = np.loadtxt("starsin_new_all_ordered.txt", usecols = (4,6,8), unpack =1) 
T_est,g_est,feh_est = np.loadtxt("test4_selfg.txt", usecols = (4,6,8), unpack =1) 
feh_clusters = [] 
g_clusters = [] 
t_clusters = [] 
name_ind2_sort1 = list(name_ind2_sort) + list([len(al)] ) 
name_ind2_sort1 = array(name_ind2_sort1) 
name_ind2_sort1 = [a - 1 for a in name_ind2_sort1] 
for i in range(len(name_ind2_sort1)-1):
  ind1 = name_ind2_sort1[i]
  ind2 = name_ind2_sort1[i+1]
  feh_clusters.append(feh_est[ind1:ind2]) 
  g_clusters.append(g_est[ind1:ind2]) 
  t_clusters.append(T_est[ind1:ind2]) 

#gc = [] 
#for jj, each in enumerate(g_clusters):
#  gc[jj,:] = each

n_bins = 10
n_bins2 = 10
x = np.random.randn(1000, 3)
#for each in g_clusters:
#  axHistx.hist(each, bins=bins,alpha  = 0.5)
#hist(g_clusters, n_bins, normed=1, histtype='bar', stacked=True)
#hist(t_clusters, n_bins2, normed=1, histtype='bar', stacked=True)
hist(feh_clusters, n_bins2, normed=1, histtype='bar', stacked=True)
#xlabel("log g", fontsize = 20) 
#xlabel("Teff ", fontsize = 20) 
xlabel("[Fe/H] ", fontsize = 20) 
ylabel("Number of Stars", fontsize = 20) 
