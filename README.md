# Procedure  

## Generating skymaps
1. Run `get_inj_params.py` with a command like ```python get_inj_params.py --N 1000``` where N = 10x the number of skymaps wanted since depending on the detector settings, only about 1/10 of the injections will generate skymaps. This will create a csv file with the name `inj_parameters_{N}.csv` using parameters from the file `BBH_15_parameters_O4_prior.json` and distance and sky location from `2700_GW_sources.csv` which has GW sources up to 2700 Mpc (lum_dist)
2. Run `generating_skymaps-v2.py` with a command like ```python generating_skymaps-v2.py --input_file inj_parameters_100.csv --n_skymaps 100 --foldername 100_maps --psdname test_high``` using the csv from the previous step. This will create n_skymaps number of skymaps in the folder named in the command which is inside the folder `/skymaps/`. Each skymap will be a fits file named based on the event number assigned in the previous step and the psd used by Bayestar.
3. Run `flatten_skymaps.py` with a command like ```python flatten_skymaps.py --mapdir ./100_maps --nside 64``` using the directory of skymaps from the previous step. This will create skymaps with a fixed pixel size that will be used for the sampling steps.
## Sampling method
1. COM:  
    a. Run `make_separate_npz.py` with a command like ```python make_separate_npz.py --mapdir /skymaps/100_maps/flattened_maps --nside_out 64 --max_distance 2000 --num_distance_bins 21 --save_folder 100_npz``` using the folder of skymaps created in the previous step. This will use the script `make_cumulative3D-stable-v4.0.0.py` to create a 3D probability distribution map of each skymap in the folder resulting in an npz file for each skymap.  
    b. Run `com_finder.py` with a command like ```python com_finder.py --mapdir ./100_npz``` using the folder of npz files created in the previous step. This will find the center of mass of the probability volume described in each npz file resulting in a distance and sky location for each npz file and therefore each skymap.  
2. Other sampling methods:  
    a. Run `make_cumulative3D-stable-v4.0.0.py` with a command like ```python make_cumulative3D-stable-v4.0.0.py --mapdir ./100_skymaps/flattened_maps --nside_out 64 --max_distance 2000 --num_distance_bins 21``` with the folder of skymaps created in the previous step. This will create a single npz file which describes the cumulative probability distribution of all skymaps. If more than ~250 skymaps are in the folder, consider making the npz files in batches and then combining them afterwards using the script `combine_npz_test.py`  
    b. Run `make_cfinput-stable-v5.1.py` with a command like ```python make_cfinput-stable-v5.1.py --map cumulative_reduced3Dmap-64nside-2000max-100step.npz --npoints 1000 --sampling dist --save``` using the npz file created in the previous step. This will create a txt file with a sky location, distance, and weight for npoints. There are 5 different sampling methods which are as follows:  
    - dist: This method flattens the 3D map and then samples based on the entire distribution from the npz file.
    - max: This method flattens the 3D map and then samples the highest probability points from the npz file.
    - distance_dist: This method first calculates the probability of each distance bin based on the values in the bin then assigns a number of points to each bin based on the probability. It then does the normal dist method inside of each bin.
    - distance_max: This method first calculates the probability of each distance bin based on the values in the bin then assigns a number of points to each bin based on the probability. It then does the normal max method inside of each bin.
    - div_dist: This method assigns an equal number of points to each distance bin then does the dist method inside of each bin.
    - div_max: This method assigns an equal number of points to each distance bin then does the max method inside of each bin.
## Correlation function  
1. Input parameters into `config.ini`. This config file will be read by the script that runs the correlation function and contains:  
    - The names of both files of data
    - Whether the GW file is in degrees or radians (the galaxy file is already in degrees)
    - The sample size of the GW data, the sample size of the galaxy data
    - The minimum redshift, the redshift bin size, and number of bins
    - The minimum luminosity distance, the luminosity distance bin size, and number of bins
    - The maximum value of theta used for the pair counts
    - The ratio of the amount of points in the random data to the actual data
    - The number of threads used for DDtheta_mocks
    - The save name of the csv that will contain the results.  
2. Run `angular_ccf_galaxy_GW-v3.0.py` which will calculate the value of the angular correlation function between different pairs of d_L bins and z bins using the jackknife method to find w and its error. Specifically, the function counts the pairs between datasets that are within theta_max of each other and uses the Landy-Szalay estimator to find w (w = (dd - dr - rd + rr) / rr) in which pair counts are normalized based on weights
3. Run `angular_ccf_analysis-v2.0.py` with a command like ```python angular_cross_correlation_analysis.py --output_file cf_result.csv --save``` which will plot the values of w using the stdevs as error bars. Each subplot represents a luminosity distance bin and shows the w values for each redshift bin where the y-axis is the value of w and the x-axis is the redshift with each point at the middle z value of each redshift bin.
