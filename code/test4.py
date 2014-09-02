pick1 = logical_and(t_f > 4000, t_f < 4300)
pick2 = logical_and(g_f > 1.4, g_f < 2.0)
pick3 = logical_and(feh_f > -0.4, feh_f < -0.25) 
pick1 = logical_and(t > 4500, t < 4700)
pick2 = logical_and(g > 1.3, g < 1.5)
pick2 = logical_and(g > 0.8, g < 1.4)
pick3 = logical_and(feh > -2.2, feh < -1.9) 
#pick1 = logical_and(t > 4500, t < 5000)
#pick2 = logical_and(g > 2.4, g < 2.8)
#pick3 = logical_and(feh > -0.4, feh < 0.2) 
pick = logical_and( logical_and(pick1, pick2), pick3) 
names[pick] 
print names[pick] 
#['M92' 'M92' 'M92' 'M15' 'M15' 'M53' 'M53' 'M53' 'M53' 'M53' 'N5466'
# 'N5466' 'N5466' 'N5466']

num_stars = len(pick[pick]) 

#testdata_pick = testdata[:,pick,:]
mean_spectra = coeffs[:,0]
testdata_pick = dataall[:,pick,:]
good =  array([(max(testdata_pick[ii,:,2])) for ii in range(8575)]) < 0.5
test = testdata_pick[:,:,1] - mean_spectra[:,None]
testnew = test[good]
colorind = ['b','b', 'b', 'r', 'r','k','k','k','k','k','m','m','m','m' ,'m','m','c']#,'k','k']
colorind = ['M92' 'M92' 'M15' 'M15' 'M53' 'M53' 'M53' 'M53' 'M53' 'N5466']
colorind = ['b' ,'b', 'r', 'r', 'k', 'k', 'k', 'k', 'k', 'm']
colorind = ['b', 'b', 'r', 'r', 'k', 'k', 'k', 'k', 'k', 'k', 'c', 'c' ] 
#colorind = ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b',
# 'b', 'b', 'b', 'b', 'r', 'c', 'c', 'c', 'c', 'c',
# 'c' ,'c', 'c', 'c', 'k', 'k', 'k', 'k', 'k', 'k',
# 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k',
# 'k', 'k', 'm', 'y']

meanspectra =  [(mean(testnew[ii,:])) for ii in range(len(testnew))]
for jj,col in zip(range(shape(testnew)[1]),colorind):
#for jj in range(shape(testnew)[1]):
  testnew = array(testnew) 
  pick2 = testnew[:,jj] > -0.3
  meanspectra = array(meanspectra) 
  z = np.polyfit(meanspectra[pick2], testnew[:,jj][pick2], 2)
  f = np.poly1d(z)
  #plot(meanspectra[pick2],f(meanspectra[pick2]))#, color = col) 
  x = meanspectra[pick2]
  y = f(meanspectra[pick2])
  x1 = x[argsort(y)]
  y1 = sort(y) 
  plot(x1,y1, color = col, linewidth = 2 ) 

print num_stars
print shape(testnew) 
