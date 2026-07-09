#!/usr/bin/env python3
# coding: utf-8

#usage: make_cumulative3D-stable-v3.1.py [-h] [--verbose] --mapdir MAPDIR [--nside_out NSIDE_OUT] [--max_distance MAX_DISTANCE] [--num_distance_bins NUM_DISTANCE_BINS] [--plot3D]
#                                        [--plot2D] [--plot3Dsum] [--plot2Dsum] [--norm_check_initial] [--norm_check_reduced] [--time] [--save]
#
#Program to create the cumulative 3D map of a set of events
#
#options:
#  -h, --help            show this help message and exit
#  --verbose             Verbose mode. Default: NO.
#  --mapdir MAPDIR       Directory of the flattened event maps in fits format.
#  --nside_out NSIDE_OUT
#                        Nside value for the (reduced) cumulative 3D map. Default: 64.
#  --max_distance MAX_DISTANCE
#                        Max distance for the 3D map (in Mpc). Default: 2000 Mpc.
#  --num_distance_bins NUM_DISTANCE_BINS
#                        Number of distance bins (linspace style). Default: 51.
#  --plot3D              Makes plots of the 3D map of all events and a slice at a random chosen distance. Default: NO.
#  --plot2D              Makes plots of 2D maps of all events at both initial and final resolution. Default: NO.
#  --plot3Dsum           Makes the plots of the cumulative 3D map. Default: NO.
#  --plot2Dsum           makes the plots of the cumulative 2D map. Default: NO.
#  --norm_check_initial  Checks the initial map normalization (Beware of large memory usage!). Default: NO.
#  --norm_check_reduced  Checks the reduced map normalization. Default: NO.
#  --time                Measures the execution time. Default: NO.
#  --save                Saves the 3D maps of all events. Default: NO.


# Usage example: ./make_cumulative3D-stable-v3.1.py --verbose --mapdir test --nside_out 8 --max_distance 8000 --num_distance_bins 2001 \
#					--norm_check_reduced --norm_check_initial --time --plot2D --plot3D --plot3Dsum --plot2Dsum 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from numba import njit, prange

import numpy as np
import healpy as hp
import pandas as pd
from scipy.stats import norm

import argparse
import time
import sys

from glob import glob


from astropy.io import fits
from astropy.table import Table as fits_Table

@njit(fastmath=True,error_model='numpy',parallel=True) 
def spheric2cartesian(angles):
	x = np.sin(angles[0])*np.cos(angles[1])
	y = np.sin(angles[0])*np.sin(angles[1])
	z = np.cos(angles[0])
	return x, y, z



def makes_everything(filename):
	root_filename = filename.split('/')[-1].split('.')[0]
	
	if verbose:
		print('Event:', root_filename)
		
	with fits.open(filename, memmap=True) as hdul:
		data = hdul[1].data
		header = hdul[1].header

		nside = header["NSIDE"]
		number_pixel = 12 * nside**2

		bayestar_map = np.empty((number_pixel, 4), dtype=np.float64)
	
		bayestar_map[:, 0] = data["PROB"]
		bayestar_map[:, 1] = data["DISTMU"]
		bayestar_map[:, 2] = data["DISTSIGMA"]
		bayestar_map[:, 3] = data["DISTNORM"]	
		
	#np.nan_to_num(bayestar_map, copy=False, nan=0.0,posinf=0.0)

	order=hp.nside2order(nside)
	pixel_area = hp.nside2pixarea(nside, 'degrees')

	if verbose:
		print('Initial nside:', nside)
		print('Initial number of pixels:', number_pixel)
		print('Pixel area:', pixel_area, 'square degrees')

	order_diff = order - order_low
	number_pixel_ratio = 4**order_diff


#prob2D = hp.pixelfunc.ud_grade(prob2D, nside, pess=False, order_in='RING', order_out='RING', power=-2, dtype=None)
#skymap_resized = np.arange(number_pixel)
#hp.get_map_size(prob2D)

	if args.plot2D:

		if verbose:
			print('Plotting the 2D map at initial resolution ('+root_filename+'_2Dhires.png)')

		hp.mollview(bayestar_map[:,0], nest=True, cbar='', title='2D map at initial resolution ('+str(number_pixel)+' pixels)', norm='hist')
	 
		plt.savefig(root_filename+'_2Dhires.png')
		plt.close()


	if args.norm_check_initial:

		if verbose:
			print('Checking the normalization of the initial 3D map')

		@njit(fastmath=True,error_model='numpy')	
		def dpdV_func_hires(i):
			r = distance_bins[i]
			dpdV = bayestar_map[:,0] * \
			   bayestar_map[:,3]/np.sqrt(2*np.pi) / bayestar_map[:,2] * \
			   np.exp(-(r-bayestar_map[:,1])**2/2/bayestar_map[:,2]**2)
			return dpdV

		@njit(fastmath=True,error_model='numpy',parallel=True) 
		def make_3Dmap_hires(distance_bin_number,prob2D):
			bayestar_3Dmap_hires = np.empty(shape=(number_pixel, distance_bin_number+1))
			bayestar_3Dmap_hires[:,0] = prob2D
			for i in prange(distance_bin_number):		
				bayestar_3Dmap_hires[:,i+1] = dpdV_func_hires(i)
			return bayestar_3Dmap_hires

		bayestar_3Dmap_hires = make_3Dmap_hires(distance_bin_number,bayestar_map[:,0])


		pixel_distance_normalization_hires = np.empty(shape=(number_pixel))
		for npix in range(number_pixel):
			pixel_distance_normalization_hires[npix] = distance_step * np.sum(bayestar_3Dmap_hires[npix,1:]*distance_bins**2)

		#The total 3D norm can also be calculated this way: norm3D = distance_step * np.sum(bayestar_3Dmap_hires[:,1:]*distance_bins**2)
	
		if verbose:
			print('Total 2D normalization of the initial map:',np.sum(bayestar_3Dmap_hires[:,0]))
			print('Total 3D normalization of the initial map:',np.sum(pixel_distance_normalization_hires))


		del bayestar_3Dmap_hires

		if verbose:
			print('Making the histogram of the pixel normalizations')

		fig = plt.figure(figsize=(12, 12))
		plt.xlabel('Normalization in distance')
		plt.ylabel('Number of pixels')
		plt.yscale('log')
		plt.title('Normalization in distance for initial map')
		plt.hist(pixel_distance_normalization_hires,bins=np.linspace(0,1,101))
		plt.savefig(root_filename+'_pixel_normalization-hist-hires.png')
		plt.close()


	if verbose:
		print('Computing the 3D map. Slicing in distance')

	bayestar_map_lowres = bayestar_map.reshape(-1, number_pixel_ratio,4)
	prob2D = np.sum(bayestar_map_lowres[:,:,0],axis=1)


	@njit(fastmath=True,error_model='numpy')
	def dpdV_func(i):
		r = distance_bins[i]
		dpdV = bayestar_map_lowres[:,:,0]* \
		   bayestar_map_lowres[:,:,3]/np.sqrt(2*np.pi) / bayestar_map_lowres[:,:,2] * \
		   np.exp(-(r-bayestar_map_lowres[:,:,1])**2/2/bayestar_map_lowres[:,:,2]**2)
		return dpdV.sum(axis=1)

	if args.plot2D:

		if verbose:
			print('Plotting the 2D map at lower resolution ('+root_filename+'_2Dlowres.png)')

		hp.mollview(prob2D, nest=True, cbar='', title='2D map at lower resolution ('+str(number_pixel_low)+' pixels)', norm='hist')

		plt.savefig(root_filename+'_2Dlowres.png')
		plt.close()


	if verbose:
		print('Reducing the resolution')

	@njit(fastmath=True,error_model='numpy',parallel=True) 
	def make_3Dmap(distance_bin_number,prob2D):
		bayestar_3Dmap = np.empty(shape=(number_pixel_low, distance_bin_number+1))
		bayestar_3Dmap[:,0] = prob2D
		for i in prange(distance_bin_number):		
			bayestar_3Dmap[:,i+1] = dpdV_func(i)
		return bayestar_3Dmap

	bayestar_3Dmap = make_3Dmap(distance_bin_number,prob2D)

	if verbose:
		print('Reduced 3D map is done')

	if args.save:
		np.savez_compressed(root_filename+'_reduced3Dmap-'+str(nside_low)+'.npz', bayestar_3Dmap=bayestar_3Dmap)
		if verbose:
			print('Saving the reduced 3D map')


	if args.plot2D:

		if verbose:
			print('Plotting a random slice of the 3D map as a visual check ('+root_filename+'_3Dlowres_slice.png)')

		random_distance = np.random.randint(1,len(distance_bins))
		hp.mollview(bayestar_3Dmap[:,random_distance], nest=True, cbar='', title='Slice of lowres 3D map at '+str(distance_bins[random_distance])+' Mpc', norm="hist")
		plt.savefig(root_filename+'_3Dlowres_slice.png')
		plt.close()



#@njit(fastmath=True,error_model='numpy',parallel=True) 
#def make_3Dmap(distance_bin_number,prob2D):
#	bayestar_3Dmap = np.empty(shape=(number_pixel_low, distance_bin_number+1))
#	bayestar_3Dmap[:,0] = prob2D
#	for i in prange(distance_bin_number):		
#		bayestar_3Dmap[:,i+1] = dpdr_func(i)
#	return bayestar_3Dmap


	if args.plot3D:

		if verbose:
			print('Passing to Cartesian coordinates to plot the 3D map')

		angles = hp.pix2ang(nside_low, np.arange(number_pixel_low),nest=True)
 
		i, j, k = spheric2cartesian(angles)
	
		if verbose:
			print('Plotting the 3D map at reduced resolution ('+root_filename+'_3Dlowres.png). It may take a while.')

		fig = plt.figure(figsize=(12, 12))
		ax = fig.add_subplot(projection='3d')
		cm = plt.get_cmap('viridis')#'YlOrBr')
		#rgba = cm(1/distance_bin_number)
		#ax.w_xaxis.set_pane_color(rgba)
		#ax.w_yaxis.set_pane_color(rgba)
		#ax.w_zaxis.set_pane_color(rgba)
		norm_alpha = 10*np.max(bayestar_3Dmap[:,1:])
 
		for d,dist in enumerate(distance_bins[1:]):
			col = bayestar_3Dmap[:,d+1]
			alpha = col/norm_alpha
			#print(np.max(alpha))
			scatter = ax.scatter3D(i*dist,j*dist,k*dist,c=col,marker='o',alpha=alpha,cmap=cm,zorder=-d*10**5,depthshade=False)#,s=distance_bins[i]**2)#*np.sin(angles[0]));
			#ax.set_box_aspect(None, zoom=2)

		fig.colorbar(scatter)
		plt.savefig(root_filename+'_3Dlowres.png')
		plt.close()


	if args.norm_check_reduced:

		if verbose:
			print('Checking the normalization of the reduced map')

		pixel_distance_normalization = np.empty(shape=(number_pixel_low))
		for npix in range(number_pixel_low):
			pixel_distance_normalization[npix] = distance_step * np.sum(bayestar_3Dmap[npix][1:]*distance_bins**2)

		#The total 3D norm can also be calculated this way: norm3D = np.sum(bayestar_3Dmap[:,1:]*distance_bins**2)

		if verbose:
			print('Total 3D normalization of the reduced map:',np.sum(pixel_distance_normalization))


		if verbose:
			print('Making the histogram of the pixel normalizations ('+root_filename+'_pixel_normalization-'+str(nside_low)+'-hist.png')


		fig = plt.figure(figsize=(12, 12))
		plt.xlabel('Normalization in distance')
		plt.ylabel('Number of pixels')
		plt.yscale('log')
		plt.title('Normalization in distance for reduced map')
		plt.hist(pixel_distance_normalization,bins=np.linspace(0,1,101))
		plt.savefig(root_filename+'_pixel_normalization-'+str(nside_low)+'-hist.png')
		plt.close()
	
	
		if verbose:
			print('Plotting the normalization of a random pixel in distance as a check ('+root_filename+'_pixel_normalization-'+str(nside_low)+'-random.png)')

		npix = np.random.randint(number_pixel_low)
		#summed = np.sum(bayestar_3Dmap,axis=0)
		fig = plt.figure(figsize=(12, 12))
		plt.xlabel('D (Mpc)')
		plt.ylabel('dP/dD ('+str(distance_step)+r' Mpc)${}^{-1}$')
		plt.title('Normalization in distance for pixel number '+str(npix))
		plt.plot(distance_bins,bayestar_3Dmap[npix][1:])
		#plt.plot(distance_bins,summed[1:])
		plt.savefig(root_filename+'_pixel_normalization-'+str(nside_low)+'-random.png')
		plt.close()

	return bayestar_3Dmap	

if __name__ == "__main__":

	global distance_step,distance_bins
	global number_pixel,number_pixel_low
	global order, order_low
	global nside_low
	global bayestar_map_lowres
	global args

	parser = argparse.ArgumentParser(description='Program to create the cumulative 3D map of a set of events',prog='make_cumulative3D-stable-v3.1.py')
	parser.add_argument('--verbose', help='Verbose mode. Default: NO.', action='count')
	parser.add_argument('--mapdir', help='Directory of the flattened event maps in fits format.', required=True)
	parser.add_argument('--nside_out', help='Nside value for the (reduced) cumulative 3D map. Default: 64.', required=False, type=int,default=64)
	parser.add_argument('--max_distance', help='Max distance for the 3D map (in Mpc). Default: 2000 Mpc.', required=False, type=int, default=2000)
	parser.add_argument('--num_distance_bins', help='Number of distance bins (linspace style). Default: 11.', required=False, type=int, default=51)
	parser.add_argument('--plot3D', help='Makes plots of the 3D map of all events and a slice at a random chosen distance. Default: NO.', action='count')
	parser.add_argument('--plot2D', help='Makes plots of 2D maps of all events at both initial and final resolution. Default: NO.', action='count')
	parser.add_argument('--plot3Dsum', help='Makes the plots of the cumulative 3D map. Default: NO.', action='count')
	parser.add_argument('--plot2Dsum', help='makes the plots of the cumulative 2D map. Default: NO.', action='count')
	parser.add_argument('--norm_check_initial', help='Checks the initial map normalization (Beware of large memory usage!). Default: NO.', action='count')
	parser.add_argument('--norm_check_reduced', help='Checks the reduced map normalization. Default: NO.', action='count')
	parser.add_argument('--time', help='Measures the execution time. Default: NO.', action='count')
	parser.add_argument('--save', help='Saves the 3D maps of all events. Default: NO.', action='count')
	

	args = parser.parse_args()

	verbose = args.verbose
	nside_low = args.nside_out
	number_pixel_low = 12*nside_low**2   #or hp.nside2npix(nside_low)
	pixel_area_low = hp.nside2pixarea(nside_low, 'degrees')
	order_low=hp.nside2order(nside_low)


	if verbose:
		print('Starting!')

	if args.time:
		start_time = time.time()

	mapdir = args.mapdir

	if verbose:
		print('Working with list of events in',mapdir)
		print('Reduced nside:',nside_low)
		print('Reduced number of pixels:', number_pixel_low)
		print('Pixel area for reduced map:', pixel_area_low, 'square degrees')

	max_distance_Mpc = args.max_distance 
	distance_bin_number = args.num_distance_bins

	distance_bins = np.linspace(0,max_distance_Mpc, distance_bin_number)
	distance_step = distance_bins[1]-distance_bins[0]
	
	if not distance_step.is_integer():
		print('The width of the distance bins is not an integer. You must choose the argument of --num_distance_bins accordingly. Try again. I am quitting now!')
		sys.exit()
	else:
		distance_step = int(distance_step)

	if verbose:
		print('Max distance:', max_distance_Mpc, 'Mpc')
		print('Distance step:', distance_step, 'Mpc')

	

	sum_reduced3D = np.zeros(shape=(number_pixel_low, distance_bin_number+1))
	skipped_maps = 0 

	map_list = glob(mapdir+'/*.fits')  + glob(mapdir+'/*.fits.gz')

	for mapfile in map_list:
		try:
			reduced3D = makes_everything(mapfile)
			sum_reduced3D+=reduced3D
		except:
			skipped_maps+=1
			print('WARNING: The file', mapfile, 'is unreadable. Skipping!')
			pass

	sum_reduced3D/=(len(map_list)-skipped_maps)
	
	if verbose:
		print('Saving the reduced 3D map. I skipped', skipped_maps, 'map(s) that were unreadable.')

	np.savez_compressed('cumulative_reduced3Dmap-'+str(nside_low)+'nside-'+str(max_distance_Mpc)+'max-'+str(distance_step)+'step.npz', cumulative_reduced3Dmap=sum_reduced3D)

	if args.plot2Dsum:

		if verbose:
			print('Plotting the 2D cumulative map (cumulative_2Dlowres.png')

		hp.mollview(sum_reduced3D[:,0], nest=True, cbar='', title='Cumulative 2D map at reduced resolution ('+str(number_pixel_low)+' pixels)', norm='hist')
	 
		plt.savefig('cumulative_2Dlowres.png')
		plt.close()

	if args.plot3Dsum:

		if verbose:
			print('Plotting the 3D cumulative map (cumulative_3Dlowres.png). It may take a while.')

		angles = hp.pix2ang(nside_low, np.arange(number_pixel_low),nest=True)
 
		i, j, k = spheric2cartesian(angles)
	
		fig = plt.figure(figsize=(12, 12))
		ax = fig.add_subplot(projection='3d')
		cm = plt.get_cmap('viridis')#YlOrBr')
		#rgba = cm(1/distance_bin_number)
		#ax.w_xaxis.set_pane_color(rgba)
		#ax.w_yaxis.set_pane_color(rgba)
		#ax.w_zaxis.set_pane_color(rgba)
		ax.set_xlim((-args.max_distance,args.max_distance))
		ax.set_ylim((-args.max_distance,args.max_distance))
		ax.set_zlim((-args.max_distance,args.max_distance))

		norm_alpha = 10*np.max(sum_reduced3D[:,1:])
 
		for d,dist in enumerate(distance_bins[1:]):
			col = sum_reduced3D[:,d+1]
			alpha = col/norm_alpha
			#print(np.max(alpha))
			scatter = ax.scatter3D(i*dist,j*dist,k*dist,c=col,marker='.',alpha=alpha,cmap=cm,zorder=-d*10**5,depthshade=True)#,s=distance_bins[i]**2)#*np.sin(angles[0]));
			#ax.set_box_aspect(None, zoom=2)

		plt.colorbar(scatter)
		plt.savefig('cumulative_3Dlowres.png')
		plt.close()


			

	if args.time:
		end_time = time.time()
		execution_time = end_time - start_time
		print('The execution time was:', execution_time, 'seconds')

	if verbose:
	   print('The End!')

	sys.exit()
