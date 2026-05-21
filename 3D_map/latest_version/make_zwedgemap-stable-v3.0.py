#!/usr/bin/env python3
# coding: utf-8

# Usage example: ./make_zwedgemap.py --map cumulative_reduced3Dmap-32.npz --distance_step 1 --verbose --phi 60

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from numba import njit, prange

import numpy as np
import healpy as hp

import argparse
import time
import sys

def plot3Dmap(filename):

	global number_pixel, distance_bin_number, distance_bins, distance_step

	fig = plt.figure(figsize=(12, 12))
	ax3D = fig.add_subplot(projection='3d')		
	cm = plt.get_cmap(args.cmp)

	filemap = np.load(filename)['cumulative_reduced3Dmap']
	number_pixel = len(filemap)	
	
	distance_bin_number = len(filemap[1,:])-1
	distance_step = args.distance_step
	distance_bins = np.arange(0,distance_bin_number,distance_step)
	
	nside = hp.npix2nside(number_pixel)
		
	if args.verbose:
		print('Plotting the wedge cumulative map. It may take a while.')

	theta_values = np.linspace(0,np.pi,args.angular_resolution,endpoint=True)
	selected_pix = hp.pixelfunc.ang2pix(nside, theta_values, args.phi*np.pi/180, nest=True, lonlat=False)

	fig = plt.figure(figsize=(24, 24))
	ax = fig.add_subplot(projection='polar')
	
	selected_col = np.empty((distance_bin_number,args.angular_resolution))
	for n in range(distance_bin_number):
		selected_col[n,:] = filemap[:,n+1][selected_pix]
		c = ax.scatter(theta_values, np.ones(args.angular_resolution)*distance_step*n, c=selected_col[n,:], marker='.', cmap='YlOrBr', alpha=args.alpha)

	ax.set_thetamin(0)
	ax.set_thetamax(180)
	ax.tick_params(axis='both', which='major', labelsize=24)
	#ax.set_xlabel('Mpc', fontsize=36)
	plt.title('Slice at longitude '+ str(args.phi) + '${}^\\circ$ degrees', fontsize=36)
	plt.savefig('zwedge_3Dlowres-'+str(nside)+'-'+str(args.phi).replace('.','d')+'deg.png')
	plt.close()

	if args.verbose:
		print('Done!')
	sys.exit()
				
		
if __name__ == "__main__":

	global args

	parser = argparse.ArgumentParser(description='',prog='bayestar_3D-current.py')
	parser.add_argument('--verbose', help='Verbose mode. Default: NO', action='count')
	parser.add_argument('--map', help='3D map file', required=True)
	parser.add_argument('--distance_step', help='Delta distance of the map (in Mpc). Default = 1.', required=False, default=1,type=int)
	parser.add_argument('--cmp', help='Color scheme. Default: YlOrBr', required=False, default='YlOrBr')
	parser.add_argument('--alpha', help='Alpha value. Default: 0.1', required=False, default=0.1)
	parser.add_argument('--angular_resolution', help='Angular resolution of the plot.', required=False, default=361, type=int)
	parser.add_argument('--phi', help='Value of the longitude angle for the slice map in degrees. Default: prime meridian = 0 degrees', required=False, default=0, type=float)
	
	args = parser.parse_args()
	
	plot3Dmap(args.map)

sys.exit()
