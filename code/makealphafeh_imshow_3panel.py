def makeim(xval,yval,wval,rad, binnum):
  fs = 20
  hist1,x2,y2,temp = hist2d(xval, yval, weights = wval, bins= binnum,cmin = 5) 
  hist1_norm,x3,y3,temp = hist2d(xval, yval, bins = binnum,cmin = 5) 
  #close("all") 
  #plt.figure(figsize = (10,6)) 
  rcParams['figure.figsize'] = 16.0, 6.0
  fig, (ax1, ax2, ax3,ax4) = plt.subplots(ncols=4, sharey = False)
  hist1a = hist1/hist1_norm
  vmin1 = 0
  vmax1 = 13
  image = hist1a 
  image_density = hist1_norm 
  bad = isnan(image)
  image[bad] = None
  s1 = ax1.imshow(image.T, interpolation="nearest" ,aspect = 'auto',origin = 'lower', cmap=plt.cm.YlGnBu,vmin =vmin1,vmax =vmax1, extent = (x3.min(), x3.max(), y3.min(), y3.max() ))
  ax1.set_xlabel('[Fe/H]',fontsize = fs) 
  ax1.set_ylabel(r'[$\alpha$/Fe]',fontsize = fs) 
  test1 = colorbar(s1,ax=ax1, orientation = 'horizontal') 
  test1.set_label("Age (Gyr) ", fontsize = fs) 
  s2 = ax2.imshow(image_density.T, interpolation="nearest" ,aspect = 'auto',origin = 'lower', cmap=plt.cm.YlGnBu, extent = (x3.min(), x3.max(), y3.min(), y3.max() ))
  ax2.set_xlabel('[Fe/H]',fontsize = fs) 
  test2= colorbar(s2,ax=ax2, orientation = 'horizontal') 
  test2.set_label("Density", fontsize = fs) 
  s2 = ax2.imshow(image_density.T, interpolation="nearest" ,aspect = 'auto',origin = 'lower', cmap=plt.cm.YlGnBu, extent = (x3.min(), x3.max(), y3.min(), y3.max() ))
  for ax in [ax1,ax2]:
    ax.set_xlim(-1.1,0.6)
    ax.set_ylim(-0.1,0.4) 
  ycut = -0.2*xval + 0.12
  pick_low = logical_and(wval < 1, yval < ycut ) 
  pick_high = logical_and(logical_and(wval > 5, wval < 14) , yval < ycut ) 
  x_values = arange(-2,1,0.1)
  y_values = -0.15*x_values + 0.08
  ax2.plot(x_values, y_values, 'm',linestyle = 'dashed')
  ax2.vlines(-0.2,0 ,0.11, 'm',linestyle = 'dashed')
  ax2.vlines(0.0, 0,0.08, 'm',linestyle = 'dashed')
  ax2.hlines(0.0, -0.2,0, 'm' , linestyle = 'dashed') 
 
  figure() 
  hist1_rad_low_norm,x4,y4,temp1 = hist2d(rad[pick_low], xval[pick_low], bins = binnum,cmin = 4) 
  hist1_rad_high_norm,x5,y5,temp2 = hist2d(rad[pick_high], xval[pick_high], bins = binnum,cmin = 4) 
  #ax.scatter(galr[pick], rc_feh[pick], color = 'k', alpha = 0.8) 
  image_rad_young_density = hist1_rad_low_norm 
  image_rad_old_density = hist1_rad_high_norm 
  s3 = ax3.imshow(image_rad_young_density.T, interpolation="nearest" ,aspect = 'auto',origin = 'lower', cmap=plt.cm.YlGnBu, extent = (x4.min(), x4.max(), y4.min(), y4.max() ))
  s4 = ax4.imshow(image_rad_old_density.T, interpolation="nearest" ,aspect = 'auto',origin = 'lower', cmap=plt.cm.YlGnBu, extent = (x5.min(), x5.max(), y5.min(), y5.max() ))
  #ax4.scatter(rad[pick_high], xval[pick_high],marker = 'x')
  test3= colorbar(s3,ax=ax3, orientation = 'horizontal') 
  test4= colorbar(s4,ax=ax4, orientation = 'horizontal') 
  test3.set_label("Density", fontsize = fs) 
  test4.set_label("Density", fontsize = fs) 
  ax4.yaxis.set_label_position("right")
  ax4.set_ylabel("[Fe/H]", fontsize = fs) 
  for ax in [ax3,ax4]:
    ax.set_xlabel('R$_{GAL}$ (kpc)' , fontsize = fs,labelpad = 10) 
  ax4.yaxis.tick_right() 
  ax3.text(5.5,0.5, "Age < 1 Gyr", fontsize = fs) 
  ax4.text(5.5,0.5, "Age > 5 Gyr", fontsize = fs) 
  ax3.set_ylim(-0.6,0.6)
  ax4.set_ylim(-0.6,0.6)
  ax4.set_xlim(4.5,16.4) 
  ax3.set_xlim(4.5,16.4) 
  ax3.set_yticklabels([])
  ax2.set_yticklabels([])
  fig.subplots_adjust(left=None, bottom=0.15)
  fig.tight_layout()
  #test = ax.scatter(galr[pick], rc_feh[pick], c = rc_alpha[pick] , s = 30, linewidth  =  0 ) 
  #fig.subplots_adjust(hspace=0)
  fig.subplots_adjust(wspace=0.01)
  savefig('/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/redclump_4panel.png', fmt = "png",bbox = 'Tight')
#  savefig('/Users/ness/new_laptop/TheCannon/TheCannon/documents/mass_and_age/plots/redclump_4panel.pdf', fmt = "pdf",bbox = 'Tight')
  return s1 

pick = galr < 8 
pick= galr < 100.
makeim(rc_feh[pick],rc_alpha[pick],rc_ages[pick], galr[pick], 30 ) 
#fig = plt.figure(figsize = (10,6)) 
#ax = fig.add_subplot(111)
#ycut = -0.2*rc_feh + 0.12
#pick = rc_ages < 2. 
##pick = logical_and(rc_ages > 12, rc_alpha < ycut ) 
#pick = logical_and(rc_ages < 1, rc_alpha < ycut ) 
#pick = logical_and(rc_ages < 1, rc_alpha < ycut ) 
#pick = logical_and(rc_ages < 1, rc_alpha < ycut ) 
#pick = logical_and(rc_ages > 12, rc_alpha < ycut ) 
#pick = logical_and(rc_ages < 1, rc_alpha < ycut ) 
#pick = logical_and(rc_ages > 4, logical_and(rc_ages < 6, rc_alpha < ycut ) ) 
##ax.scatter(galr[pick], rc_feh[pick], color = 'k', alpha = 0.8) 
##test = ax.scatter(galr[pick], rc_feh[pick], c = rc_alpha[pick] , s = 30, linewidth  =  0 ) 
#test = ax.scatter(galr[pick], rc_feh[pick], c = rc_alpha[pick] , s = 30, linewidth  =  0 ) 
#colorbar(test, ax=ax) 
#ax.set_xlim(4,18)
#ax.set_ylim(-1.0,0.55)
#ax.set_xlabel("R$_{GC}$ (kpc) ", fontsize = 28,labelpad = 8 ) 
#ax.set_ylabel("[Fe/H]", fontsize = 28) 
#fig.subplots_adjust(left=None, bottom=0.18)
draw() 
