import numpy as np
import matplotlib.pyplot as pl
import scipy.stats as st


#data = np.random.multivariate_normal((0,0), [[0.8, 0.05], [0.05, 0.7]], 1000)
#datab = np.random.multivariate_normal((-.1,-.1), [[1.0, 0.0], [0.0, 1.0]], 1000)

## read in beta-coal points in file rnorm0
data = []
with open("rnorm0") as h:
    for l in h:
        data.append([float(z) for z in l.split()])

h.close()
del h
data = np.array(data)

x=data[:, 0]
y=data[:,1]

xmin, xmax= -3.0, 3.0
ymin, ymax= -3.0, 3.0

xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([ xx.ravel(), yy.ravel()])

values = np.vstack([x,y])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
f = f/np.max(f)

del data, x,y, values, kernel

## read in exp growth data points in file rnorm1 
data=[]
with open("rnorm1") as h:
    for l in h:
       data.append([float(z) for z in l.split()])

h.close()
del h
data = np.array(data)

x = data[:, 0]
y = data[:, 1]
values = np.stack([x,y])
kernel = st.gaussian_kde(values)
g = np.reshape(kernel(positions).T, xx.shape)
g = g/np.max(g)


fig = pl.figure()
ax = fig.gca()
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
##cfset = ax.contourf(xx, yy, f )
##ax.imshow( np.rot90(f), cmap='Greys', extent=[xmin,xmax,ymin,ymax])

## contour plot
cset= ax.contour(xx, yy, f, colors=['.1','0.6'],  levels=[0.01, 0.11] )
cset= ax.contour(xx, yy, g, colors=['.1','0.6'],  levels=[0.01, 0.11], linestyles = 'dashed' )


ax.clabel(cset, inline=1, fontsize=10)


## example for adding actual data points to contour plot
pl.plot(0.5, 0.5, 'b+' )
pl.plot(-.5, 0.5, 'k<' )
pl.plot(-1.5, 0.5, 'ko' )


pl.show()


del data, x, y, f, g
