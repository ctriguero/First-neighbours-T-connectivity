#!/usr/bin/env python

from pylab import *
from numpy import *
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
     
rc('text', usetex=True)
rc('font', family='serif')
figure(1, figsize=(16,10))

z = loadtxt('Distances.dat')

n, bins, patches = hist(z, 400, normed=1,  linewidth=0.1, edgecolor='gray', color='red', alpha=0.75)#, alpha=0.75 ,edgecolor='black' ,, histtype='stepfilled'

of = open('PDF_D.dat', mode='w')

for i in range(n.shape[0]):
    s = '%12.6f %12.6f\n' % ( bins[i] , n[i] )
    of.write(s)
    print s,
of.close()
    
    
xlabel(r'First neighbours distance, $D$',fontsize=30)
ylabel(r'Probability den. dist., $P(D)$',fontsize=30)
title(r'First neighbours distance probability density distribution',fontsize=25, color='gray')

#Axis Range
#       xmin xmmax ymin ymax
plt.axis([0.0,50.0,0.0,1.3],fontsize=40)

#Text
#plt.text(0.05, 30, r'$\sigma_{\rm S}=|\beta-d|-\displaystyle{\frac{1}{2}}$',fontsize=40)
#axvline(x=0.131122448979592,linewidth=1.5, color='black' , linestyle='dashed') # Limit of well
#axvline(x=0.36887756,linewidth=1.5, color='black' , linestyle='dashed') # Limit of well

#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#img=mpimg.imread('x.png')
#imgplot = plt.imshow(img, (x[0.05],y[30]), xycoords='data', frameon=False)


#def get_bars(request)
#    fig = Figure(facecolor='#F0F0F0',figsize=(4.6,4))
#    ax1 = fig.add_subplot(111,ylabel="Valeur",xlabel="Code",autoscale_on=True)
#    ax1.bar(ind,values,width=width, color='#FFCC00',edgecolor='#B33600',linewidth=1)
#    canvas = FigureCanvas(fig)
#    response = HttpResponse(content_type='x.png')
#    canvas.print_png(response)
#    return response


savefig('PDF_D.pdf')
#savefig('PDF_Lambda.png')
#show()

