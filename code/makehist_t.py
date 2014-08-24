params_labels = [params[:,0], params[:,1], params[:,2] , covs_params[:,0,0]**0.5, covs_params[:,1,1]**0.5, covs_params[:,2,2]**0.5 ]  
testdir = "/Users/ness/Downloads/Apogee_raw/calibration_fields/4332/apogee/spectro/redux/r3/s3/a3/v304/4332/"
file2 = '4332_data_all_more.txt'
file2in = testdir+file2
t,t_err,g,g_err,feh,feh_err = loadtxt(file2in, usecols = (1,2,3,4,5,6), unpack =1) 

pick = logical_and(g > 0, logical_and(t_err < 300, feh > -4.0) ) 
cval = ['k', 'b', 'r'] 
input_ASPCAP = [t, g, feh, t_err, g_err, feh_err] 
input_ASPCAP = [t, g, feh, t_err, g_err, feh_err] 
pick = logical_and(input_ASPCAP[1] > 0, logical_and(input_ASPCAP[3] < 300, input_ASPCAP[2] > -4.0) ) 
x = input_ASPCAP[0][pick]
y = params_labels[0][pick]
yerr = params_labels[1+2][pick]
xerr = input_ASPCAP[1+2][pick]
c_val = input_ASPCAP[1][pick]

starsin = open('starsin_new_all_ordered.txt','r') 
al = starsin.readlines()
names = []
for each in al:
  names.append(each.split()[1]) 
unames = unique(names) 
starind = arange(0,len(names), 1) 
name_ind = [] 
name_ind2 = [] 
names = array(names) 
for each in unames:
  takeit = each == names 
  name_ind.append(starind[takeit][-1]+1. ) 
  name_ind2.append(starind[takeit][0] ) 
names = array(names) 
for each in unames:
  takeit = each == names 
  name_ind.append(starind[takeit][-1]+1. ) 
  name_ind2.append(starind[takeit][0]+1. ) 
name_ind_sort = sort(name_ind)
name_ind2_sort = sort(name_ind2)

#zip(sort(name_ind), unames[argsort(name_ind)]) 

T_est,g_est,feh_est = np.loadtxt("starsin_new_all_ordered.txt", usecols = (4,6,8), unpack =1) 
feh_clusters = [] 
g_clusters = [] 
t_clusters = [] 
name_ind2_sort1 = list(name_ind2_sort) + list([len(al)] ) 
name_ind2_sort1 = array(name_ind2_sort1) 
for i in range(len(name_ind2_sort1)-1):
  ind1 = name_ind2_sort1[i]
  ind2 = name_ind2_sort1[i+1]
  feh_clusters.append(feh_est[ind1:ind2]) 
  g_clusters.append(g_est[ind1:ind2]) 
  t_clusters.append(T_est[ind1:ind2]) 

nullfmt   = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(1, figsize=(11,11))

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx,alpha = 0.4)
axHisty = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
s1 = axScatter.scatter(x, y, c = c_val, s = 30)
cind = c_val
cax = axes([0.11, 0.55, 0.012, 0.18])
a,b,c1 = axScatter.errorbar(x, y,yerr= yerr,marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
a,b,c2 = axScatter.errorbar(x, y,xerr= xerr,marker='',ls='',zorder=0, fmt = None,elinewidth = 1,capsize = 0)
c_T = colorbar(s1,cax=cax) 
c_T.set_label('APSCAP logg', fontsize = 16,labelpad = 6) 
g_color = c_T.to_rgba(cind) 
c1[0].set_color(g_color)
c2[0].set_color(g_color)


# now determine nice limits by hand:
binwidth = 50
xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
lim = ( int(xymax/binwidth) + 1) * binwidth

axScatter.plot([0,6000], [0,6000],'k',linewidth = 2) 

bins = np.arange(-lim, lim + binwidth, binwidth)
bins = np.arange(3300,6000,50) 
for each in t_clusters:
  axHistx.hist(each, bins=50,alpha  = 0.5,range = (3300,6000) )

axHisty.hist(y, bins=50, orientation='horizontal',range = (0,6000))

axScatter.set_xlim(3300,5500) 
axScatter.set_ylim(0,5500) 
axHisty.set_ylim(0,5500) 
axHistx.set_xlim(3300,5500) 
axScatter.set_xlabel("ASPCAP Teff", fontsize = 20)
axScatter.set_ylabel("NHR+ Teff", fontsize = 20) 
axHistx.set_title("calibrators: Teff", fontsize = 20 ) 
#axHisty.set_ylabel("[Fe/H] hisogram of NHR+ [Fe/H]", fontsize = 20) 
axHisty.yaxis.set_label_position("right")
draw()
show()
savefig('/Users/ness/Downloads/Apogee_Raw/calibration_apogeecontinuum/documents/plots/cal_teff.eps', transparent=True, bbox_inches='tight', pad_inches=0)
