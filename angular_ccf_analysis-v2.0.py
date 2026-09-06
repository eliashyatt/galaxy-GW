import pandas as pd
import numpy as np
import argparse
import configparser
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter
from scipy import integrate
from scipy.interpolate import UnivariateSpline
import scipy.stats as stats
import astropy.constants as const


def get_config(output_file):
    config = {}

    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()

            if not line:
                break
        
            if line == '[Config]':
                continue
        
            if '=' in line:
                key, value = line.split('=', 1)
                config[key.strip()] = value.strip()
    return config

def get_ws_and_stdevs(output_file):
    df = pd.read_csv(
        output_file,
        skiprows=18   # skips lines before "ws,stdevs"
    )

    ws = df["ws"].to_numpy()
    stdevs = df["stdevs"].to_numpy()

    #splits w and stdevs into the different d_L bins
    split_w = np.array_split(np.array(ws),d_L_bin_num)
    split_stdevs = np.array_split(np.array(stdevs),d_L_bin_num)
    return split_w, split_stdevs

#finds factors closest to the square root of the input number (used to make the subplots look better)
def find_factors(num):
    sqrt_num = np.sqrt(num)
    if sqrt_num % 1 == 0:
        high_factor = int(sqrt_num)
        low_factor = int(sqrt_num)
    else:
        low_num = num // sqrt_num
        while num % low_num != 0:
            low_num = low_num - 1
        low_factor = int(low_num)
        high_factor = int(num / low_factor)
    return low_factor, high_factor

#defines gaussian
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


def find_z_means():
    GW_df = pd.read_csv('GW_inj_parameters/2700_GW_sources.csv')
    z_means = []
    for i in range(d_L_bin_num):
        d_L_min = luminosity_bin_min + i * luminosity_bin_size
        d_L_max = d_L_min + luminosity_bin_size
        filtered_GW = GW_df[GW_df['luminosity_distance'] < d_L_max]
        filtered_GW = filtered_GW[filtered_GW['luminosity_distance'] >= d_L_min]
        d_L_bin = filtered_GW
        z_means.append(d_L_bin['redshift'].mean())
    return z_means

def plot_w(split_w, split_stdevs, split_w2, split_stdevs2, split_w3, split_stdevs3, z_means):
    #each subplot is a d_L bin and graphs the w value for each redshift bin with a gaussian fit
    z_min = redshift_min + redshift_bin_size/2
    z_max = redshift_min + z_bin_num * redshift_bin_size - redshift_bin_size/2
    z_values = np.array(np.linspace(z_min,z_max,z_bin_num))
    print(z_values)
    sp_rows, sp_columns = find_factors(d_L_bin_num)
    fig,axs = plt.subplots(sp_rows,sp_columns,figsize=(10,8))
    z_means_est = []
    stdev_est = []

    for m in range(sp_rows):
        for n in range(sp_columns):
            i = m*sp_columns + n
            d_L_min = luminosity_bin_size*i+luminosity_bin_min
            d_L_max = d_L_min + luminosity_bin_size
            title = str(d_L_min) + "-" + str(d_L_max)
            x_data = z_values
            y_data = np.ravel(split_w[i])
            y_error = split_stdevs[i]
            z_mean = z_means[i]

            if split_w2:
                y_data2 = np.ravel(split_w2[i])
                y_error2 = split_stdevs2[i]
            if split_w3:
                y_data3 = np.ravel(split_w3[i])
                y_error3 = split_stdevs3[i]

            #y_smooth = savgol_filter(y_data, 9, 2)
            #start_cut = 0
            #peaks, _ = find_peaks(y_smooth[start_cut:], prominence=0.1)
            #peak_idx = start_cut + peaks[np.argmax(y_smooth[start_cut:][peaks])]

            inverse_error = np.where(y_error>0)
            #print(inverse_error)
            z_stdev = np.std((y_data[inverse_error] - np.average(y_data[inverse_error]))/y_error[inverse_error])
            #print(z_stdev)
        
            #peaks, _ = find_peaks(y_data, height=None, prominence=[0, z_stdev])
            #print(peaks)
            #peak_idx = peaks[np.argmax(y_data[peaks])]
            #mu0 = x_data[peak_idx]
            #print(x_data[np.argmax(y_data)],mu0,x_data[peaks])
            #A0 = y_data[peak_idx]

            #popt, pcov = curve_fit(gaussian, x_data, y_data, sigma=y_error, absolute_sigma=True, p0 = [A0, mu0, redshift_bin_size])
            #A_fit, mu_fit, sigma_fit = popt
            #z_means_est.append(mu_fit)
            #stdev_est.append(abs(sigma_fit))
            #g_x_fit = np.linspace(min(x_data), max(x_data), 500)
            #g_y_fit = gaussian(g_x_fit, *popt)

            if split_w2 == None and split_w3 == None:
                axs[m,n].errorbar(x_data, y_data, yerr=y_error, fmt='.', color='b')
            else:
                axs[m,n].errorbar(x_data, y_data, yerr=y_error, fmt='.', color='b', label=f'{output_file}', alpha=0.5)
            
            if split_w2:
                axs[m,n].errorbar(x_data, y_data2, yerr=y_error, fmt='.', color='r', label=f'{output_file2}', alpha=0.5)
            if split_w3:
                axs[m,n].errorbar(x_data, y_data3, yerr=y_error, fmt='.', color='g', label=f'{output_file3}', alpha=0.5)
            
            axs[m,n].axvline(z_mean,linestyle='--',color='k')
            #axs[m,n].plot(g_x_fit, g_y_fit, label='Gaussian fit', color='r')
            axs[m,n].set_title(title)
            axs[m,n].set_xlabel('z')
            axs[m,n].set_ylabel('w')
            axs[m,n].set_ylim(-1, 1)
            #axs[m,n].legend()

    handles, labels = plt.gca().get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    #plt.show()

    if args.save:
        #print("Saving plot as a PNG")
        png_name = args.output_file.split('.')[0]
        if split_w2 == None and split_w3 == None:
            png_savename = f'{png_name}.png'
            #plt.savefig(f'{png_name}.png')
        else:
            if split_w3:
                png_name2 = args.output_file2.split('.')[0]
                png_name3 = args.output_file3.split('.')[0]
                png_savename = f'{png_name}_and_{png_name2}_and_{png_name3}_comparison.png'
                #plt.savefig(f'{png_name}_and_{png_name2}_and_{png_name3}_comparison.png')
            else:
                png_name2 = args.output_file2.split('.')[0]
                png_savename = f'{png_name}_and_{png_name2}_comparison.png'
                #plt.savefig(f'{png_name}_and_{png_name2}_comparison.png')
        
        print(f"Saving plot as {png_savename}")
        plt.savefig(png_savename)

    plt.show()            


if __name__ == "__main__":
    global redshift_min, redshift_bin_size, z_bin_num, d_L_bin_num, luminosity_bin_min, luminosity_bin_size, args
    parser = argparse.ArgumentParser(description="Plots the angular cross-correlation function w and its standard deviation for different luminosity distance bins.")
    parser.add_argument("--output_file", help="CSV file containing the output from the angular cross-correlation function.", required=True)
    parser.add_argument("--output_file2", help="CSV file for a second dataset to compare.", required=False)
    parser.add_argument("--output_file3", help="CSV file for a third dataset to compare.", required=False)
    parser.add_argument("--save", help="Save the plot as a PNG file.", action="count")
    args = parser.parse_args()

    output_file = args.output_file
    config = get_config(output_file)

    output_file2 = args.output_file2
    output_file3 = args.output_file3

    print("Starting...")

    if output_file3:
        if output_file2 == None:
            raise ValueError('No second output file listed but there is one in the third slot. Consider moving it to the second slot instead.')

    redshift_min = float(config["redshift_min"])
    redshift_bin_size = float(config["redshift_bin_size"])
    z_bin_num = int(config["z_bin_num"])
    d_L_bin_num = int(config["d_l_bin_num"])
    luminosity_bin_min = float(config["luminosity_bin_min"])
    luminosity_bin_size = float(config["luminosity_bin_size"])

    print("Finding values from first file")

    split_w, split_stdevs = get_ws_and_stdevs(output_file)
    z_means = find_z_means()

    if output_file2:
        print("Finding values from second file")
        split_w2, split_stdevs2 = get_ws_and_stdevs(output_file2)
    else:
        split_w2, split_stdevs2 = None, None

    if output_file3:
        print("Finding values from third file")
        split_w3, split_stdevs3 = get_ws_and_stdevs(output_file3)
    else:
        split_w3, split_stdevs3 = None, None

    print("Plotting...")
    plot_w(split_w, split_stdevs, split_w2, split_stdevs2, split_w3, split_stdevs3, z_means)
    




