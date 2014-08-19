a = open("starsin_new_all_ordered.txt", 'r')
al = a.readlines() 
bl = []
for each in al:
  bl.append(each.split()[0]) 
each = bl[0]
#each = each.strip('\n')
a = pyfits.open(each) 
b = pyfits.getheader(each) 
print median(a[1].data), median(a[2].data) 
print np.max(a[1].data), np.max(a[2].data) 
print np.min(a[1].data), np.min(a[2].data) 
