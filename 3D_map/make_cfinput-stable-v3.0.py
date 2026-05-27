#!/usr/bin/env python3
# coding: utf-8

# Usage example: python make_cfinput-stable-v3.0.py --verbose --map cumulative_reduced3Dmap-256nside-2000max-5.0step.npz --distance_step 5 --npoints 100000 --sampling dist --plotspherical --plotwedge 90

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from numba import njit, prange

import numpy as np
import healpy as hp

import argparse
import time
import sys

import gc


def read_data():
	global number_pixel, distance_bin_number, distance_bins, distance_step, distance_max, dim
	global nside

	data = np.load(args.map)['cumulative_reduced3Dmap']
	
	#data = data[:48,:12]					#FOR TESTING 	
	#np.savetxt('test.csv',data,delimiter=',')		#FOR TESTING

	number_pixel = len(data)
	distance_bin_number = len(data[1,:])-1
	distance_step = args.distance_step
	distance_bins = np.arange(0,(distance_bin_number)*(distance_step),distance_step)
	distance_max = distance_bins[-1]
		
	nside = hp.npix2nside(number_pixel)

	dim = number_pixel * distance_bin_number
	#print(distance_bin_number);print(distance_bins)	#FOR TESTING

	if args.verbose:
		print('The 2D map has',str(number_pixel),'pixels and',str(distance_bin_number),'distance bins with size',str(distance_step),'Mpc from', \
		str(distance_bins[0]), 'Mpc to',str(distance_max), 'Mpc for a total of',str(dim),'cells.')	

	return data


@njit(fastmath=True,error_model='numpy',parallel=True)
def select_data(data, npoints, sampling):
	flattened_data = data[:,1:].T.flatten()
	
	#distances = np.empty(dim)
	#for n in prange(distance_bin_number):
	#	distances[n*number_pixel:(n+1)*number_pixel] = n*distance_step
	distances = np.repeat(np.arange(distance_bin_number) * distance_step, number_pixel)

	total = np.sum(flattened_data)
	probabilities = flattened_data / total

	if sampling > 0:
	    cdf = np.cumsum(probabilities)
	    cdf /= cdf[-1]
	    selected_points = np.searchsorted(cdf,np.random.random(npoints))
	else:
	    npoints = min(npoints, len(flattened_data))
	    selected_points = np.argpartition(-flattened_data,npoints - 1)[:npoints]
		
	selected_weights = flattened_data[selected_points]
	selected_distances = distances[selected_points]
	selected_points = selected_points % number_pixel

	return selected_points, selected_distances, selected_weights
	
def dist2z(dist):
	H0 = 72  #km/s/MPc
	c = 3*10**5 #km/s
	return dist*H0/c

#@njit(fastmath=True,error_model='numpy',parallel=True) 
#def make_distance_bins(i,j,k):
#    dim = number_pixel * distance_bin_number
#    dist_i = np.empty(dim)
#    dist_j = np.empty(dim)
#    dist_k = np.empty(dim)
#    for n in prange(distance_bin_number):
#        dist_i[n*number_pixel:(n+1)*number_pixel] = i*n*distance_step
#        dist_j[n*number_pixel:(n+1)*number_pixel] = j*n*distance_step
#       dist_k[n*number_pixel:(n+1)*number_pixel] = k*n*distance_step
#   return dist_i,dist_j,dist_k	


def plot3Dmap(ras, decs, zs, weights):

	fig = plt.figure(figsize=(12, 12))
	ax3D = fig.add_subplot(projection='3d')		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced 3D map (cumulative_3Dreducedmap-nside-max-step-points.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=ax3D)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	alpha = 0.5*(weights-np.min(weights))/(np.max(weights)-np.min(weights))

	scatter = ax3D.scatter(ras, np.cos(decs), zs, c=weights, marker=',', cmap=cm, alpha=alpha)#, depthshade=False)
	#scatter = ax3D.tricontourf(ras, np.cos(decs), zs, marker='.', cmap=cm)

	rgba = [0,0,0,0]
#	plot_axis_max = args.max_distance 
	ax3D.set_xlim((0,2*np.pi))
	ax3D.set_ylim(-1,1)#(0,np.pi))
#	ax3D.set_zlim((-plot_axis_max,plot_axis_max))
	ax3D.set_xlabel('R.A. (rad)', fontsize=16)
	ax3D.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_3Dreducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()

@njit(fastmath=True,error_model='numpy',parallel=True) 
def spheric2cartesian(theta,phi):
	x = np.sin(theta)*np.cos(phi)
	y = np.sin(theta)*np.sin(phi)
	z = np.cos(theta)
	return x, y, z

def plot2Dmap(ras, decs, zs, weights):

	fig = plt.figure(figsize=(12, 12))
	ax2D = fig.add_subplot()		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced 2D map (cumulative_2Dreducedmap-nside-max-step-npoints.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=ax2D)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	alpha = 0.05*(zs-np.min(zs))/(np.max(zs)-np.min(zs))
	
	if np.isnan(alpha).any():
		alpha = 1
		print('Warning! Using alpha = 1 for the 2D map (all points are on a shell at fixed z).')
	
	size = 10**(3*(zs-np.max(zs))/(np.min(zs)-np.max(zs)))

	if np.isnan(size).any():
		size = 1
		print('Warning! Using size = 1 for the 2D map  (all points are on a shell at fixed z).')

	color = (weights-np.min(weights))/(np.max(weights)-np.min(weights))

	scatter = ax2D.scatter(ras, np.cos(decs), c=color, s=size, marker=',', cmap=cm, alpha=alpha,zorder=-100*np.mean(zs)) #, depthshade=False)

	rgba = [0,0,0,0]
	ax2D.set_xlim((0,2*np.pi))
	ax2D.set_ylim(-1,1)#(0,np.pi))
	ax2D.set_xlabel('R.A. (rad)', fontsize=16)
	ax2D.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_2Dreducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()

def plotspheremap(ras, decs, zs, weights):

	fig = plt.figure(figsize=(12, 12))
	axsphere = fig.add_subplot(projection='3d')		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced map on a sphere (cumulative_spherereducedmap-nside-max-step-npoints.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=axsphere)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	alpha = 0.05*(zs-np.min(zs))/(np.max(zs)-np.min(zs))
	
	if np.isnan(alpha).any():
		alpha = 1
		print('Warning! Using alpha = 1 for the sphere map (all points are on a shell at fixed z).')
	
	size = 10**(3*(zs-np.max(zs))/(np.min(zs)-np.max(zs)))

	if np.isnan(size).any():
		size = 1
		print('Warning! Using size = 1 for the sphere map  (all points are on a shell at fixed z).')

	color = (weights-np.min(weights))/(np.max(weights)-np.min(weights))

	x,y,z= spheric2cartesian(decs,ras)
	scatter = axsphere.scatter(x,y,z, c=color, s=size, marker=',', cmap=cm, alpha=alpha,zorder=-100*np.mean(zs))
	
	rgba = [0,0,0,0]
	axsphere.set_xlim(-1,1)
	axsphere.set_ylim(-1,1)
	#axsphere.set_xlabel('R.A. (rad)', fontsize=16)
	#axsphere.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_spherereducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()

def plotsphericalmap(ras, decs, distances, weights):

	fig = plt.figure(figsize=(12, 12))
	axspherical = fig.add_subplot(projection='3d')		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced spherical map (cumulative_sphericalreducedmap-nside-max-step-npoints.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=axspherical)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	color = (weights-np.min(weights))/(np.max(weights)-np.min(weights))

	x,y,z = spheric2cartesian(decs,ras) * distances
	
	
	mask = ((x >= -distance_max/np.sqrt(2)) & (x <= distance_max/np.sqrt(2)) &(y >= -distance_max/np.sqrt(2)) & (y <= distance_max/np.sqrt(2)) \
		&(z >= -distance_max/np.sqrt(2)) & (z <= distance_max/np.sqrt(2)))
	x = x[mask]
	y = y[mask]
	z = z[mask]
	weights = weights[mask]
	color = color[mask]
	distances = distances[mask]

	alpha = 0.5*(distances-np.min(distances))/(np.max(distances)-np.min(distances))	
	if np.isnan(alpha).any():
		alpha = 1
		print('Warning! Using alpha = 1 for the spherical map.')
		
	size = 2
	scatter = axspherical.scatter(x,y,z, c=color, s=size, marker=',', cmap=cm, alpha=alpha,zorder=-np.mean(distances))

#	size = 6
#	scatter = axspherical.scatter(x,y,z,c=color,s=size,marker='o',cmap=cm,alpha=0.25,linewidths=0)

	
	rgba = [0,0,0,0]
	axspherical.set_xlim(-distance_max/np.sqrt(2),distance_max/np.sqrt(2))
	axspherical.set_ylim(-distance_max/np.sqrt(2),distance_max/np.sqrt(2))
	axspherical.set_zlim(-distance_max/np.sqrt(2),distance_max/np.sqrt(2))
	
	#axsphere.set_xlabel('R.A. (rad)', fontsize=16)
	#axsphere.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_sphericalreducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()


def plotwedgemap(data):

    if args.verbose:
        print('Plotting full wedge map...')

    angular_resolution = 361
    theta_values = np.linspace(0, np.pi, angular_resolution, endpoint=True)

    selected_pix = hp.pixelfunc.ang2pix(nside,theta_values,args.plotwedge * np.pi / 180,nest=True,lonlat=False)

    r_values = np.arange(distance_bin_number) * distance_step

    # data[:, 1:] is assumed (skip first column)
    data_slice = data[:, 1:]  # shape: (npix, nbins)

    # result shape: (nbins, angular_resolution)
    selected_col = data_slice[selected_pix, :].T


    fig = plt.figure(figsize=(24, 24))
    ax = fig.add_subplot(projection='polar')

    cmap = plt.get_cmap(args.cmp)

    for n in range(distance_bin_number):
        ax.scatter(theta_values,np.full_like(theta_values, r_values[n]),c=selected_col[n],marker='.',cmap=cmap,alpha=args.alpha)

    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.tick_params(axis='both', which='major', labelsize=24)

    plt.title(f'Slice at longitude {args.plotwedge}$^\\circ$ degrees',fontsize=36)

    plt.savefig(
        f'cumulative_wedgefullmap-{nside}-{distance_max}-{distance_step}-'
        f'{args.npoints}-{str(args.plotwedge).replace(".","d")}deg.png')
    plt.close()
		
		
if __name__ == "__main__":

	global args

	parser = argparse.ArgumentParser(description='',prog='make_cfinput.py')
	parser.add_argument('--verbose', help='Verbose mode. Default: NO', action='count')
	parser.add_argument('--map', help='3D map file', required=True)
	parser.add_argument('--distance_step', help='Delta distance of the map (in Mpc). Default = 1.', required=False, default=1,type=int)
	parser.add_argument('--cmp', help='Color scheme. Default: YlOrBr', required=False, default='YlOrBr')
	parser.add_argument('--alpha', help='Alpha value. Default: 0.1', required=False)
	parser.add_argument('--npoints', help='Number of points in the weigth map.', required=True, type=int)
	parser.add_argument('--plot3D', help='Plots the 3D map. Default: NO', action='count')
	parser.add_argument('--plot2D', help='Plots the 2D map. Default: NO', action='count')
	parser.add_argument('--plotsphere', help='Plots the map on a sphere. Default: NO', action='count')
	parser.add_argument('--plotspherical', help='Plots the 3D map in spherical coordinates. Default: NO', action='count')
	parser.add_argument('--save', help='Saves the input file.', action='count')
	parser.add_argument('--sampling', help='Point sampling scheme (dist or max). Default: max', required=False, default='max')
	parser.add_argument('--plotwedge',nargs='?',type=float, const=60, default=None, help='Plot wedge slice at given longitude (degrees). If used without a value, defaults to 60.')	
	args = parser.parse_args()
	
	if args.verbose:
		print('Reading the cumulative map.')
	data = read_data()
		
	if args.npoints >= dim:
		print('Warning! You are selecting either all cells or an unphysical number of cells. This will not work. Please choose a value of --npoints smaller than', \
			str(args.npoints),'and rerun.')
		sys.exit()
		
	if args.verbose:
		print('Selecting ' + str(args.npoints) + ' pixels according to the weight probability distribution ('+str(np.round(args.npoints/dim*100,1))+'% of the total).')

	if args.sampling == 'dist':
		sampling = 1 
	else:
		sampling = 0
	
	points, distances, weights = select_data(data,args.npoints,sampling)
	
	if args.verbose:
		print('Computing coordinates and redshifts.')
	zs = dist2z(distances)
	decs,ras = hp.pix2ang(nside, points, nest=True)
	
	if args.save:
		if args.verbose:
			print('Saving the input file.')
		header = 'Highest ' + str(args.npoints) + ' probability pixels from ' + args.map + '\n' + 'dec,ra,z,weight'
		matrix_output = np.c_[decs,ras,distances,weights]
		matrix_output = matrix_output[matrix_output[:,3].argsort()][::-1]   #to sort in descending order
		np.savetxt('cf_input_'+args.map.split('.')[0]+'-'+str(args.npoints)+'.txt', matrix_output, fmt='%.6f,%.6f,%.6f,%.8e', delimiter=',', newline='\n', header=header, footer='', comments='# ', encoding=None)

	#print(ras);print(decs);print(points);print(distances);print(zs);print(weights)	#FOR TESTING

	if args.plot3D:
		plot3Dmap(ras, decs, distances, weights)
	
	if args.plot2D:
		plot2Dmap(ras, decs, distances, weights)
	
	if args.plotsphere:
		plotspheremap(ras, decs, zs, weights)

	if args.plotspherical:
		plotsphericalmap(ras, decs, distances, weights)

	if args.plotwedge:
		plotwedgemap(data)


sys.exit()
