from pylab import *

#data
time = array(g[pick]) 
signal = time**2
error = ones(len(time))*1000

figure(1)
sc = scatter(time,signal,s=20,c=time)

#create colorbar according to the scatter plot
clb = colorbar(sc)

#create errorbar plot and return the outputs to a,b,c
a,b,c = errorbar(t[pick],tme[pick],yerr=t_err[pick],marker='',ls='',zorder=0)

#convert time to a color tuple using the colormap used for scatter
time_color = clb.to_rgba(time)

#adjust the color of c[0], which is a LineCollection, to the colormap
c[0].set_color(time_color)

fig = gcf()
fig.show()
