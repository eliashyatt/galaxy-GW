from ligo.skymap.postprocess import find_greedy_credible_levels
import healpy as hp
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
#import matplotlib.gridspec as gridspec
from ligo.skymap import plot, io
from ligo.skymap.io import fits
from ligo.skymap import postprocess
from astropy import units as u
import os
import pandas as pd

approx_name = "SEOBNRv4PHM"
contour_levels = [50,90]

Events = []
Real_50 =[]
Real_90 =[]
Inj_10s_before_50 =[]
Inj_10s_before_90 =[]

print(os.listdir("skymaps/"))

#for fin in os.listdir("skymaps/10s_before/"):
for fin in os.listdir("skymaps/"):
    if fin.endswith(".fits"):
        
        print (fin)
        f = "skymaps/"+fin
        event = fin[0:len(fin)-5] # before
#    event = fin[0:len(fin)-22] # before
        print (event)
        Events.append(event)
    
 ### Fake skymaps
        skymap, metadata = fits.read_sky_map(f, nest=None) #also check if skymap normalized and if sum is normalized
        cls = 100 * postprocess.find_greedy_credible_levels(skymap)
    
        nside = hp.npix2nside(len(skymap))
        deg2perpix = hp.nside2pixarea(nside, degrees=True)
        probperdeg2 = skymap / deg2perpix
#    vmax = probperdeg2.max()
    
        ax = plt.axes(projection=('astro hours mollweide'))
        ax.imshow_hpx((probperdeg2, 'ICRS'), nested=True,
#                   vmax=vmax
                       cmap="cylon")
        text = []
        pp = np.round(contour_levels).astype(int)
        ii = np.round(
                np.searchsorted(np.sort(cls), contour_levels) * deg2perpix).astype(int)
        Inj_10s_before_50.append(ii[0])
        Inj_10s_before_90.append(ii[1])
   
        for i, p in zip(ii, pp):
            text.append(u'{:d}% area: {:d} deg²'.format(p, i, grouping=True))
   
        
        ax.text(1, 1.05, '\n'.join(text), transform=ax.transAxes, ha='right',
                    fontsize=8)  
  
        ax.grid()
        ax.set_title(event,fontsize=8)
        plot.outline_text(ax)
        plt.savefig('./skymaps_plot/'+event+".png")
        plt.close()

