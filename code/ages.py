a = open("starsin_new_all_ordered.txt",'r')
a = open("test4_selfg.txt",'r')
a = open("test14.txt",'r')
al = a.readlines()
bl = [] 
name = [] 
for each in al:
  name.append(each.split()[1]) 
name_ind = [] 
name_ind_first = [] 
name_u = unique(name) 
name = array(name) 
for each in name_u:
  pick = name == each
  name_ind.append(arange(0,len(name))[pick] ) 
  name_ind_first.append(arange(0,len(name))[pick] [0]) 

name_ind_sort = sort(name_ind_first) 
name_order = name_u[argsort(name_ind_first)] 
#name_ind_sort = [np.int(a)-1 for a in name_ind_sort] 
#cluster_ages = names[name_ind_sort[0:-1]]
name = array(name)
    
a_M107 = 14.0
a_M13 = 11.65
a_M15 = 12.0
a_M2 = 13
a_M3 = 11.4
a_M5 = 13.0
a_M53 = 12.67
a_M67 = 3.2
a_M71 = 10.0
a_M92 = 14.
a_N188 = 5.
a_N2158 = 1.05
a_N2420 = 2.
a_N4147 = 14.
a_N5466 = 12.5
a_N6791 = 5
a_N6819 = 2.5
a_N7789 = 1.6
a_Pleiades = 0.15
d = {}
# Fill in the entries one by one
# Static lists for purpose of illustration
names = name_u
ages = [a_M107, a_M13, a_M15, a_M2, a_M3, a_M5, a_M53, a_M67, a_M71, a_M92, a_N188, a_N2158, a_N2420, a_N4147, a_N5466, a_N6791, a_N6819, a_N7789, a_Pleiades] 
ages_dict = {}
for i in range(len(names)):
    ages_dict[names[i]] = ages[i]

age_vals = ones(len(name) ) 
for one,two in zip(name_u, name_ind):
 # one = "a_"+one
  age_vals[two] = ages_dict[one]

savetxt("ages.txt", age_vals, fmt = '%s' ) 


n_each = [0] + list(name_ind_sort) 
num_each = diff(n_each) 
#array(['a_M15', 'a_M53', 'a_N5466', 'a_N4147', 'a_M13', 'a_M2', 'a_M3', 'a_M5', 'a_M107',
#       'a_M71', 'a_N2158', 'a_N2420', 'a_Pleiades', 'a_N7789', 'a_M67', 'a_N6819',
#       'a_N188', 'a_N6791']
