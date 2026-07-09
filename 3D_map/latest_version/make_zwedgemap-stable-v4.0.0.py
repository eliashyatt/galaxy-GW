#!/usr/bin/env python3
# coding: utf-8

#usage: make_zwedgemap-stable-v3.1.py [-h] [--verbose] --map MAP [--cmp CMP] [--alpha ALPHA] [--res RES] [--phi PHI]
#
#options:
#  -h, --help           show this help message and exit.
#  --verbose            Verbose mode. Default: NO.
#  --map MAP            3D map in npz format created with make_cumulative3D-stable-v3.1.py.
#  --cmp CMP            Color scheme. Default: YlOrBr.
#  --alpha ALPHA        Alpha value. Default: 0.1.
#  --res RES		Angular resolution of the plot in degrees. Default: 1 degree.
#  --phi PHI            Value of the longitude angle for the slice map in degrees. Default: prime meridian = 0 degrees.

# Usage example: ./make_zwedgemap-stable-v3.1.py --map cumulative_reduced3Dmap-32.npz --verbose --phi 60

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import healpy as hp

import argparse
import time
import sys

def read_data():
	global number_pixel, distance_bin_number, distance_bins, distance_step, distance_max, dim
	global nside

	data = np.load(args.map)['cumulative_reduced3Dmap']
		
	parameters = args.map.split('-')
	distance_max = int(parameters[2][:-3])
	distance_step = int(parameters[3][:-8])
	distance_bin_number = distance_max/distance_step
	
	if not distance_bin_number.is_integer():
		print('The number of bins does not seem to be an integer. Something really wrong in your map. I am quitting now!')
		sys.exit()
	else:
		distance_bin_number = int(distance_bin_number)
		
	nside = int(parameters[1][:-5])

	number_pixel = len(data)

	if distance_bin_number != len(data[1,:])-2:
		print('The number of bins I read from the map does not match the namefile. Something really wrong in your map. I am quitting now!')
		sys.exit()
	
	if nside != hp.npix2nside(number_pixel):
		print('The nside I read from the map does not match the namefile. Something really wrong in your map. I am quitting now!')
		sys.exit()	
	
	distance_bins = np.linspace(0,distance_max,distance_bin_number+1)

	dim = number_pixel * distance_bin_number

	if args.verbose:
		print('The 2D map has',str(number_pixel),'pixels and',str(distance_bin_number),'distance bins with size',str(distance_step),'Mpc from', \
		str(distance_bins[0]), 'Mpc to',str(distance_bins[-1]), 'Mpc for a total of',str(dim),'cells.')	

	return data

def plotwedgemap(data):

	angular_resolution = int(np.round(360/args.res,0))+1
	theta_values = np.linspace(0, np.pi, angular_resolution, endpoint=True)

	selected_pix = hp.pixelfunc.ang2pix(nside,theta_values,args.phi * np.pi / 180,nest=True,lonlat=False)

	r_values = np.arange(distance_bin_number) * distance_step

	# data[:, 1:] is assumed (skip first column)
	data_slice = data[:, 1:]  # shape: (npix, nbins)

	selected_col = data_slice[selected_pix, :].T

	fig = plt.figure(figsize=(24, 24))
	ax = fig.add_subplot(projection='polar')

	cmap = plt.get_cmap(args.cmp)

	for n in range(distance_bin_number):
		ax.scatter(theta_values,np.full_like(theta_values, r_values[n]),c=selected_col[n],marker='.',cmap=cmap,alpha=args.alpha)

	ax.set_thetamin(0)
	ax.set_thetamax(180)
	ax.tick_params(axis='both', which='major', labelsize=24)

	plt.title(f'Slice at longitude {args.phi}$^\\circ$. Resolution = {args.res}$^\\circ$.',fontsize=36)

	plt.savefig(
		f'cumulative_wedgefullmap-{nside}-{distance_max}-{distance_step}-{str(args.phi).replace(".","d")}deg.png')
	plt.close()

	
				
		
if __name__ == "__main__":

	global args

	parser = argparse.ArgumentParser(description='Program to make a 2D slice of a 3Dmap at a given longitude.',prog='make_zwedgemap-stable-v3.1.py')
	parser.add_argument('--verbose', help='Verbose mode. Default: NO.', action='count')
	parser.add_argument('--map', help='3D map in npz format created with make_cumulative3D-stable-v3.1.py.', required=True)
	parser.add_argument('--cmp', help='Color scheme. Default: YlOrBr.', required=False, default='YlOrBr')
	parser.add_argument('--alpha', help='Alpha value. Default: 0.1.', required=False, default=0.1)
	parser.add_argument('--res', help='Angular resolution of the plot in degrees. Default: 1 degree.', required=False, default=1, type=float)
	parser.add_argument('--phi', help='Value of the longitude angle for the slice map in degrees. Default: prime meridian = 0 degrees.', required=False, default=0, type=float)
	
	args = parser.parse_args()

	if args.verbose:
		print('Reading the cumulative map...')
		
	data = read_data()
	
	if args.verbose:
		print('Plotting the full wedge map...')

	plotwedgemap(data)

	if args.verbose:
		print('done!')

	sys.exit()
