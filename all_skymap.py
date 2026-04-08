import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import os
import sys
import argparse
import pandas as pd
from scipy.special import eval_legendre
from  ligo.skymap.io.fits import read_sky_map
import pickle
#===============================================================
#----define variables ----------------------------------------
#===============================================================

parser = argparse.ArgumentParser()
parser.add_argument('--directory', action='store',default = '../skymap_database/synthetic_skymaps/all_skymaps_sets_34/', help='Directory of fits files')
parser.add_argument('--outdir', action='store',default = './', help='savepath')
parser.add_argument('--savename', action='store',default = 'synthetic', help='savepath')

args = parser.parse_args()
folder_path = args.directory + '/'

def main():
  
    file_paths = (os.path.join(root, file) for root, _, files in os.walk(folder_path) for file in files)
    event_filenames = list((file_paths))
    
    print (len(event_filenames), "will be normalized and stacked")
    all_skymap_norm = stack_all_skymaps_norm(event_filenames)  

    output = args.outdir
    check_mkdir_floder(output)
    savepath = output +"/all_map_"+args.savename 
    plot_skymaps(all_skymap_norm, savepath) 


    return


def stack_all_skymaps_norm(event_filenames):

    set_map_resolution = 256
   # res_cutoff = 1.0/30.0  #rad  - corresponding to freq ~ 200 for LV
    number_pixel = hp.nside2npix(set_map_resolution)
    
    skymaps = []
    for event in event_filenames:
        print (event)
        skymap,header = read_sky_map(event, nest=False, distances=False, moc=False)
        skymap_resized = hp.pixelfunc.ud_grade(skymap, set_map_resolution, pess=False, order_in='RING', order_out='RING', power=-2, dtype=None)
        norm_map_resized = np.sum(skymap_resized)
        skymaps.append(skymap_resized)
    all_skymap = [sum(i) for i in zip(*skymaps)]
    all_skymap_norm = np.divide(all_skymap, len(event_filenames))
    
    return (all_skymap_norm)


def plot_skymaps(skymap,save_path):

    hp.mollview(
        skymap,
        cbar='',
        title='',
        norm="hist",
        min=0,
        max=0.1)

    plt.savefig(save_path +'.png')
    hp.write_map(save_path+'.fits', skymap, overwrite=True)
    
    return
   


def check_mkdir_floder(f):
    check_folder = os.path.isdir(f)    
    # If folder doesn't exist, then create it.
    if not check_folder:
        os.makedirs(f)
        print("creating folder : ", f)
    else:
        pass
        
    return

if __name__ == "__main__":    
    main()                                    


