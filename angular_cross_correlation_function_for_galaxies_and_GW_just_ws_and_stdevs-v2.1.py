#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
os.environ["CORRFUNC_USE_AVX512"] = "0"
os.environ["CORRFUNC_USE_AVX2"] = "0"
os.environ["CORRFUNC_USE_AVX"] = "0"
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
import astropy.units as u
import healpy as hp
from joblib import Parallel, delayed
import argparse
import configparser
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter
from scipy import integrate
from scipy.interpolate import UnivariateSpline
from scipy.stats import gaussian_kde
import scipy.stats as stats
import astropy.constants as const
import time


def parse_config():    
    parser = argparse.ArgumentParser()    
    parser.add_argument('--config', required=True)    
    args = parser.parse_args(["--config", "config.ini"])
    config = configparser.ConfigParser() 
    config.read(args.config) 
    return config

def make_clean():    
    global galaxy_csv_name, GW_csv_name, deg_or_rad, GW_data_size, galaxy_data_size, redshift_min, redshift_bin_size, luminosity_bin_min, luminosity_bin_size, theta_max, nside, z_bin_num, d_L_bin_num, r2dr, nthreads, savename
    galaxy_csv_name = config['Config'].get('galaxy_csv_name')
    GW_csv_name = config['Config'].get('GW_csv_name')
    deg_or_rad = config['Config'].get('deg_or_rad')
    GW_data_size = config['Config'].getint('GW_data_size')    
    galaxy_data_size = config['Config'].getint('galaxy_data_size')    
    redshift_min = config['Config'].getfloat('redshift_min')    
    redshift_bin_size = config['Config'].getfloat('redshift_bin_size')
    luminosity_bin_min = config['Config'].getfloat('luminosity_bin_min')
    luminosity_bin_size = config['Config'].getfloat('luminosity_bin_size')
    theta_max = config['Config'].getfloat('theta_max')
    nside = config['Config'].getint('nside')
    z_bin_num = config['Config'].getint('z_bin_num')
    d_L_bin_num = config['Config'].getint('d_L_bin_num')
    r2dr = config['Config'].getint('random_2_data_ratio')
    nthreads = config['Config'].getint('nthreads')
    savename = config['Config'].get('savename')


#get variables from config file
config = parse_config()
make_clean()

#csv data
galaxy_df_full = pd.read_csv(galaxy_csv_name)
GW_df_full = pd.read_csv(GW_csv_name)

#distance label
if GW_csv_name[-3:] == 'txt':
    GW_df_full = pd.read_csv(GW_csv_name, comment="#", names=["dec", "ra", "distance", "weight"])
    GW_df_full['luminosity_distance'] = GW_df_full['distance']

possible_names = [
    'luminosity_distance',
    'distance',
    'luminosity distance'
]

distance_col_galaxy = next(
    (col for col in possible_names if col in galaxy_df_full.columns),
    None
)
if distance_col_galaxy is None:
    raise ValueError("No distance column found")

galaxy_df_full['redshift'] = galaxy_df_full[distance_col_galaxy]

distance_col_GW = next(
    (col for col in possible_names if col in GW_df_full.columns),
    None
)
if distance_col_GW is None:
    raise ValueError("No distance column found")
    
GW_df_full['luminosity_distance'] = GW_df_full[distance_col_GW]



print("Checking parameters...")
#testing parameters
if GW_data_size > len(GW_df_full):
    raise ValueError("GW_data_size too big. Input must be <=", len(GW_df_full))
if int(GW_data_size) != GW_data_size or GW_data_size < 1:
    raise ValueError("GW_data_size must be a positive integer")
if galaxy_data_size > len(galaxy_df_full):
    raise ValueError("galaxy_data_size too big. Input must be <=", len(galaxy_df_full))
if int(galaxy_data_size) != galaxy_data_size or galaxy_data_size < 1:
    raise ValueError("galaxy_data_size must be a positive integer")
if redshift_min < 0:
    raise ValueError("redshift_min must be positive")
if redshift_bin_size < 0:
    raise ValueError("redshift_bin_size must be positive")
if luminosity_bin_min < 0:
    raise ValueError("luminosity_bin_min must be positive")
if luminosity_bin_size < 0:
    raise ValueError("luminosity_bin_size must be positive")
if theta_max < 0:
    raise ValueError("theta_max must be positive")
if nside <= 0 or (nside & (nside - 1)) != 0:
    raise ValueError("nside must be a power of 2")
if int(z_bin_num) != z_bin_num or z_bin_num < 1:
    raise ValueError("z_bin_num must be a positive integer")
if int(d_L_bin_num) != d_L_bin_num or d_L_bin_num < 1:
    raise ValueError("d_L_bin_num must be a positive integer")
if r2dr < 1:
    raise ValueError("r2dr must be larger than 1")
if deg_or_rad != 'rad' and deg_or_rad != 'deg':
    raise ValueError("deg_or_rad must be either 'deg' or 'rad'")

#converting ra and dec to deg
if deg_or_rad == 'rad':
    GW_df_full['RA_deg'] = GW_df_full['ra'] * 180 / np.pi
    
    if np.max(GW_df_full['dec']) > (np.pi / 2):
        GW_df_full['Dec_deg'] = 90.0 - (GW_df_full['dec'] * 180 / np.pi)
    else:
        GW_df_full['Dec_deg'] = GW_df_full['dec']
    
    galaxy_df_full['RA_deg'] = galaxy_df_full['ra'] * 180 / np.pi
    if np.max(galaxy_df_full['dec']) > (np.pi / 2):
        galaxy_df_full['Dec_deg'] = 90.0 - (galaxy_df_full['dec'] * 180 / np.pi)
    else:
        galaxy_df_full['Dec_deg'] = galaxy_df_full['dec']


#accounts for no weight column
if 'weight' not in GW_df_full.columns:
    GW_df_full['weight'] = 1.0

#defines the redshift bin for an iteration
def redshift_bin(i):
    z_min = redshift_min + i * redshift_bin_size
    z_max = z_min + redshift_bin_size
    z_bin = galaxy_df[galaxy_df['redshift'] > z_min]
    z_bin = z_bin[z_bin['redshift'] < z_max]
    return z_bin


#defines the luminosity distance bin for an iteration
def luminosity_distance_bin(i):
    d_L_min = luminosity_bin_min + i * luminosity_bin_size
    d_L_max = d_L_min + luminosity_bin_size
    filtered_GW = GW_df[GW_df['luminosity_distance'] < d_L_max]
    filtered_GW = filtered_GW[filtered_GW['luminosity_distance'] > d_L_min]
    d_L_bin = filtered_GW
    return d_L_bin


#creates an array of ra and dec values the length of num
def random_func(num):
    ra_random = np.random.uniform(0,360,num)
    u = np.random.uniform(-1,1,num)
    dec_random = np.degrees(np.arcsin(u))
    random_coords = np.array([np.array(ra_random),np.array(dec_random)])
    return random_coords


#creates bins using the random dataset that is proportional to the size of the "real" dataset
def random_bins(random_coords, bin_size):
    random_indices = np.random.choice(len(random_coords[0]), bin_size, replace=False)

    random_ra = random_coords[0][random_indices]
    random_dec = random_coords[1][random_indices]
    random_coords_filtered = np.array([random_ra, random_dec])
    return random_coords_filtered


#creates random weights for the random GW dataset using KDE
def random_weights(GW_weights, d_L_bin):
    kde = gaussian_kde(GW_weights)
    GW_random_weights = kde.resample(len(d_L_bin)*r2dr).flatten()
    return GW_random_weights

#cross correlation function
def dd_rr_func(theta_bin, galaxy, GW, GW_weights):
    galaxy_ra, galaxy_dec = galaxy
    GW_ra, GW_dec = GW
    dd_rr = DDtheta_mocks(
        autocorr=0,
        nthreads=nthreads, 
        binfile=theta_bin,
        RA1=galaxy_ra,
        DEC1=galaxy_dec,
        RA2=GW_ra,
        DEC2=GW_dec,
        weights2=GW_weights,
        weight_type='pair_product'
    )
    return dd_rr


#normalizes dd and rr based on the bin sizes
def norm(dd, rr, gdr, wdr, z_bin, d_L_bin, GW_weights, GW_random_weights):
    gd_w = np.sum(np.ones(len(z_bin)))
    gr_w = np.sum(np.ones(len(z_bin)*r2dr))
    gwd_w = np.sum(GW_weights)
    gwr_w = np.sum(GW_random_weights)

    dd_norm = (dd['npairs'] * dd['weightavg']) / (gd_w * gwd_w)
    rr_norm = (rr['npairs'] * rr['weightavg']) / (gr_w * gwr_w)
    gdr_norm = (gdr['npairs'] * gdr['weightavg']) / (gd_w * gwr_w)
    wdr_norm = (wdr['npairs'] * wdr['weightavg']) / (gr_w * gwd_w)
    return dd_norm, rr_norm, gdr_norm, wdr_norm

#Landy-Szalay estimator
def w_func(dd_norm, rr_norm, gdr_norm, wdr_norm):
    #w = (dd_norm / rr_norm) - 1
    w = (dd_norm - gdr_norm - wdr_norm + rr_norm) / rr_norm
    return w


#w calculation without jackknife estimate
def w_no_jk(z_bin, d_L_bin, all_data):
    galaxy, GW, galaxy_random, GW_random, GW_weights, GW_random_weights = all_data
    dd = dd_rr_func(theta_bin, galaxy, GW, GW_weights)
    rr = dd_rr_func(theta_bin, galaxy_random, GW_random, GW_random_weights)
    gdr = dd_rr_func(theta_bin, galaxy, GW_random, GW_random_weights)
    wdr = dd_rr_func(theta_bin, galaxy_random, GW, GW_weights)
    dd_norm, rr_norm, gdr_norm, wdr_norm = norm(dd, rr, gdr, wdr, z_bin, d_L_bin, GW_weights, GW_random_weights)
    w = w_func(dd_norm, rr_norm, gdr_norm, wdr_norm)
    return w


#assigns healpix regions
def assign_healpix_region(ra_deg, dec_deg):
    theta = np.radians(90.0 - dec_deg)  # co-latitude
    phi = np.radians(ra_deg)            # longitude
    return hp.ang2pix(nside, theta, phi)


#assigns region numbers to each data point based on the nside resolution
def regions(coords):
    ra, dec = coords
    assigned_regions = assign_healpix_region(ra, dec)
    return assigned_regions


#masks the data by removing one region
def masking(coords, region, i):
    ra, dec = coords

    mask = region != i

    sample_ra = ra[mask]
    sample_dec = dec[mask]
    sample = np.array([sample_ra, sample_dec])
    return sample

#weight mask
def mask_weights(weights, regions, i):
    return weights[regions != i]


#performs the w calculation after removing one region at a time, returns the w value without jackknife if no points are found in the region
def jk_loop(all_data, all_regions, i, w0):
    galaxy, GW, galaxy_random, GW_random, GW_weights, GW_random_weights = all_data
    galaxy_regions, GW_regions, galaxy_random_regions, GW_random_regions = all_regions

    if i in galaxy_regions or i in GW_regions:
        galaxy_sample = masking(galaxy, galaxy_regions, i)
        GW_sample = masking(GW, GW_regions, i)
        galaxy_random_sample = masking(galaxy_random, galaxy_random_regions, i)
        GW_random_sample = masking(GW_random, GW_random_regions, i)
        GW_weights_sample = mask_weights(GW_weights, GW_regions, i)
        GW_random_weights_sample = mask_weights(GW_random_weights, GW_random_regions, i)

        assert len(GW_sample[0]) == len(GW_weights_sample)

        galaxy_random_ra, galaxy_random_dec = galaxy_random_sample

        dd = dd_rr_func(theta_bin, galaxy_sample, GW_sample, GW_weights_sample)
        rr = dd_rr_func(theta_bin, galaxy_random_sample, GW_random_sample, GW_random_weights_sample)
        gdr = dd_rr_func(theta_bin, galaxy_sample, GW_random_sample, GW_random_weights_sample)
        wdr = dd_rr_func(theta_bin, galaxy_random_sample, GW_sample, GW_weights_sample)

        dd_norm, rr_norm, gdr_norm, wdr_norm = norm(dd, rr, gdr, wdr, galaxy_sample[0], GW_sample[0], GW_weights_sample, GW_random_weights_sample)
        w = w_func(dd_norm, rr_norm, gdr_norm, wdr_norm)
    else:
        w = w0
    return w


#runs the jk_loop function and calculates the mean and variance from the results
def jackknife(all_data, all_regions, w0):
    galaxy, GW, galaxy_random, GW_random, GW_weights, GW_random_weights = all_data
    galaxy_regions, GW_regions, galaxy_random_regions, GW_random_regions = all_regions

    n_regions = len(galaxy_regions)

    jackknife_ws = Parallel(n_jobs=1)(
        delayed(jk_loop)(all_data, all_regions, i, w0) for i in range((nside**2)*12)
    )

    jackknife_mean = np.nanmean(jackknife_ws)
    jackknife_variance = ((n_regions - 1) / n_regions) * np.sum((jackknife_ws - jackknife_mean) ** 2)
    return jackknife_variance


#returns the w value and jackknife variance for pair of a redshift bin and a luminosity distance bin
def run_code(z_iter, d_L_iter):
    z_bin = redshift_bin(z_iter)
    d_L_bin = luminosity_distance_bin(d_L_iter)
    #print(len(z_bin))
    #print(len(d_L_bin))

    galaxy_ra, galaxy_dec = np.array(z_bin["RA_deg"]), np.array(z_bin["Dec_deg"])
    galaxy = np.array([galaxy_ra, galaxy_dec])
    GW_ra, GW_dec, GW_weights = np.array(d_L_bin["RA_deg"]), np.array(d_L_bin["Dec_deg"]), np.array(d_L_bin['weight'])
    GW = np.array([GW_ra, GW_dec])

    galaxy_random = random_bins(all_galaxy_random, len(z_bin)*r2dr)
    GW_random = random_bins(all_GW_random, len(d_L_bin)*r2dr)
    if GW_weights.all() == 1.0:
        GW_random_weights = np.ones(len(d_L_bin)*r2dr)
    else:
        GW_random_weights = random_weights(GW_weights, d_L_bin)

    inputs = [galaxy, GW, galaxy_random, GW_random]
    results = Parallel(n_jobs=4)(delayed(regions)(data) for data in inputs)
    galaxy_regions, GW_regions, galaxy_random_regions, GW_random_regions = results

    all_data = [galaxy, GW, galaxy_random, GW_random, GW_weights, GW_random_weights]
    all_regions = [galaxy_regions, GW_regions, galaxy_random_regions, GW_random_regions]

    w = w_no_jk(z_bin, d_L_bin, all_data)
    var = jackknife(all_data, all_regions, w)
    return w, var
    

print("Starting")
start_time = time.time()

print("Sampling data sets")

gamma = stats.gamma
#galaxy_a, galaxy_loc, galaxy_scale = 2.4211, 0, 0.1753
galaxy_a, galaxy_loc, galaxy_scale = 1.9094, 0, 1319.36
galaxy_probs = gamma.pdf(galaxy_df_full['redshift'], galaxy_a, galaxy_loc, galaxy_scale)
galaxy_probs = galaxy_probs / np.sum(galaxy_probs)

#random sample of GW and galaxy data sets
if galaxy_data_size == len(galaxy_df_full):
    galaxy_df = galaxy_df_full
    print("Galaxy sample size same as population")
else:
    z_least_bin_count = 0
    while z_least_bin_count < 2:
        galaxy_df_indices = np.random.choice(len(galaxy_df_full), galaxy_data_size, replace=False, p=galaxy_probs)
        galaxy_df = galaxy_df_full.iloc[galaxy_df_indices]
        z_least_bin_count = (galaxy_df['redshift'] < redshift_bin_size).sum()

if GW_data_size == len(GW_df_full):
    GW_df = GW_df_full
    print("GW sample size same as population")
else:
    #pdf for GW events from GWTC-4
    GW_a, GW_loc, GW_scale = 1.9094, 0, 1319.36
    GW_probs = gamma.pdf(GW_df_full['luminosity_distance'], GW_a, GW_loc, GW_scale)
    GW_probs = GW_probs / np.sum(GW_probs)
    #GW sample
    GW_df_indices = np.random.choice(len(GW_df_full), GW_data_size, replace=False, p=GW_probs)
    GW_df = GW_df_full.iloc[GW_df_indices]


print("GW sources:", len(GW_df))
print("Galaxies:", len(galaxy_df))
print("Number of luminosity distance bins:", d_L_bin_num)
print("Number of redshift bins:", z_bin_num)


ws = []
variances = []

#establishes the theta bin for the cross correlation function
theta_max_deg = np.degrees(theta_max)
theta_bin = np.array([0, theta_max_deg])

galaxy_size = len(galaxy_df)
GW_size = len(GW_df)

#establishes the random datasets
galaxy_ra_random, galaxy_dec_random = random_func(galaxy_size*r2dr)
GW_ra_random, GW_dec_random = random_func(GW_size*r2dr)
all_galaxy_random = np.array([galaxy_ra_random, galaxy_dec_random])
all_GW_random = np.array([GW_ra_random, GW_dec_random])

num_done = 0
total_bin_pairs = z_bin_num * d_L_bin_num

#loops the run_code function over all pairs of redshift bins and luminosity distance bins
for d_L_iter in range(d_L_bin_num):
    for z_iter in range(z_bin_num):
        w, var = run_code(z_iter, d_L_iter)
        ws.append(w)
        variances.append(var)
        num_done += 1
        print(num_done, "/", total_bin_pairs)
    #print(d_L_iter + 1, "/", d_L_bin_num)

#print(ws)
#print(len(ws))
#print(variances)
#print(len(variances))


stdevs = np.sqrt(variances)

ws_arr = np.concatenate(ws)
stdevs_arr = np.array(stdevs)


#print(ws_arr)
#print(stdevs_arr)


output = np.column_stack((ws_arr, stdevs_arr))

from io import StringIO
buffer = StringIO()
config.write(buffer)
config_string = buffer.getvalue()

header_text = f"{config_string}\nws,stdevs"

np.savetxt(savename, output, delimiter=",", header=header_text, comments="")

end_time = time.time()
execution_time = end_time - start_time
print("The execution time was:", execution_time, "seconds")

print("Done!")

