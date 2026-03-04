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

#for fin in os.listdir("skymaps/10s_before/"):
for fin in os.listdir("skymaps/"):
    print (fin)
    f = "skymaps/"+fin
    event = fin[0:len(fin)-23] # before
#    event = fin[0:len(fin)-22] # before
    print (event)
    Events.append(event)
    
 ### Fake skymaps
    skymap, metadata = fits.read_sky_map(f, nest=None)
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
    ax.set_title(event+"_10s_before_from_lalapp",fontsize=8)
    plot.outline_text(ax)
    plt.savefig('./skymaps_plot/'+event+"_10s_before_lalapp.png")
    plt.close()
    
##### Real skymaps

    f_real = '/home/yanyan.zheng/others/skymap/O3a_pe_samples/all_skymaps/' + event +'_C01:'+ approx_name +'.fits'
    sm_real, md_real = fits.read_sky_map(f_real, nest=None)
    cls_real = 100 * postprocess.find_greedy_credible_levels(sm_real)
    
    nside_real = hp.npix2nside(len(sm_real))
    deg2perpix_real = hp.nside2pixarea(nside_real, degrees=True)
    probperdeg2_real = sm_real / deg2perpix_real

    ax2 = plt.axes(projection=('astro hours mollweide'))
    ax2.imshow_hpx((probperdeg2_real, 'ICRS'), nested=True, 
                   cmap="cylon")
    text_real= []
#    pp = np.round(contour_levels).astype(int)
    ii_real = np.round(
            np.searchsorted(np.sort(cls_real), contour_levels) * deg2perpix_real).astype(int)
    Real_50.append(ii_real[0])
    Real_90.append(ii_real[1])
    
    for i, p in zip(ii_real, pp):
        text_real.append(u'{:d}% area: {:d} deg²'.format(p, i, grouping=True))

        
    ax2.text(1, 1.05, '\n'.join(text_real), transform=ax2.transAxes, ha='right',
                fontsize=8)
    ax2.grid()
    ax2.set_title(event+"_real_skymap",fontsize=8)
    plot.outline_text(ax2)
    plt.savefig('./skymaps_plot/'+event+".png")
    plt.close()

# Comparison skymaps

#    fig = plt.figure(figsize=(8, 8))

    ax1 = plt.axes(projection='astro hours mollweide')
    ax1.contour_hpx((cls_real, 'ICRS'), nested=md_real['nest'],
            colors='red', linewidths=0.5,levels=contour_levels)
    ax1.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors='blue', linewidths=0.5, levels=contour_levels)
    ax1.set_title(event+"_real(red)_vs_lalinjection(blue)",fontsize=8)
    text = []
    for i_real,i_inj, p in zip(ii_real, ii, pp):
        text.append(u'Diff in {:d}% area: {:d} deg²'.format(p, i_real-i_inj, grouping=True))
    
    
    ax1.text(1.05, 1.05, '\n'.join(text), transform=ax1.transAxes, ha='right',
                fontsize=8)
    plot.outline_text(ax1)
    plt.savefig('./skymaps_plot/'+event+"_10s_before_lalapp_comparison.png",fontsize=8)
    plt.close()
    
    
    # compare with the previous one
    f_fake2 = '/home/yanyan.zheng/others/skymap/O3a_HLV_skymaps/skymaps/10s_before/'+event+'_10s_before_lalapp.fits'
    sm_fake2, md_fake2 = fits.read_sky_map(f_fake2, nest=None)
    cls_fake2 = 100 * postprocess.find_greedy_credible_levels(sm_fake2)

    ax2 = plt.axes(projection='astro hours mollweide')
    ax2.contour_hpx((cls_real, 'ICRS'), nested=md_real['nest'],
            colors='red', linewidths=0.5,levels=contour_levels)
    ax2.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors='blue', linewidths=0.5, levels=contour_levels)
    ax2.contour_hpx(
            (cls_fake2, 'ICRS'), nested=md_fake2['nest'],
            colors='green', linewidths=0.5, levels=contour_levels)
    ax2.set_title(event+"_real(red)_vs_actual_psd(blue)_vs design(green)",fontsize=8)
    plot.outline_text(ax1)
    plt.savefig('./skymaps_plot/'+event+"_10s_before_lalapp_comparison_psd.png",fontsize=8)
    plt.close()

df = pd.DataFrame(columns = Events,index=['Real_50%','Inj_10s_before_50%','Diff_10s_before_50%','Real_90%','Inj_10s_before_90%','Diff_10s_before_90%'])
df.loc['Real_50%',Events] = Real_50
df.loc['Real_90%',Events] = Real_90
df.loc['Inj_10s_before_50%',Events] = Inj_10s_before_50
df.loc['Inj_10s_before_90%',Events] = Inj_10s_before_90
df.loc['Diff_10s_before_50%',Events] = df.loc['Real_50%',Events] - df.loc['Inj_10s_before_50%',Events]
df.loc['Diff_10s_before_90%',Events] = df.loc['Real_90%',Events] - df.loc['Inj_10s_before_90%',Events]
df.to_csv('./skymaps_plot/O3a_HLV_checking_symaps_10s.csv')
