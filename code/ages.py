name_ind_sort
name_ind
names = array(names)
names[name_ind_sort] 
name_ind_sort = [np.int(a)-1 for a in name_ind_sort)] 
cluster_ages = names[name_ind_sort[0:-1]]
       
a_M107 = 14.0
a_M13 = 11.65
a_M15 = 12.0
a_M2 = 13
a_M3 = 11.4
a_M5 = 10.6/13.
a_M53 = 12.67
a_M67 = 4.
a_M71 = 9.5
a_M92 = 14.
a_N188 = 7.
a_N2158 = 1.05
a_N2420 = 2.
a_N4147 = 14.
a_N5466 = 12.5
a_N6791 = 5
a_N6819 = 2.5
a_N7789 = 1.6
a_Pleaides = 0.15 
n_each = [0] + list(name_ind_sort) 
num_each = diff(n_each) 
array(['a_M15', 'a_M53', 'a_N5466', 'a_N4147', 'a_M13', 'a_M2', 'a_M3', 'a_M5', 'a_M107',
       'a_M71', 'a_N2158', 'a_N2420', 'a_Pleiades', 'a_N7789', 'a_M67', 'a_N6819',
       'a_N188', 'a_N6791']
