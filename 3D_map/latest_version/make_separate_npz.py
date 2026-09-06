import numpy as np
import pandas as pd
from glob import glob
import os
import shutil
import argparse
import subprocess
from pathlib import Path


def create_npz(filename, mapdir):
    workspace = os.path.join(mapdir, 'workspace')
    os.makedirs(workspace, exist_ok=True)

    basename = os.path.basename(filename)
    newname = os.path.join(workspace, basename)
    shutil.move(filename, newname)

    print(newname)

    event_num = basename[:5]

    #print(os.getcwd())

    make_cumulative = [
        'python3',
        'make_cumulative3D-stable-v4.0.0.py',
        '--verbose',
        '--mapdir', workspace,
        '--nside_out', str(nside_out),
        '--max_distance', str(max_distance),
        '--num_distance_bins', str(num_distance_bins)
    ]
    subprocess.run(make_cumulative, check=True)

    shutil.move(newname, filename)

    npzfilename = Path('cumulative_reduced3Dmap-'+str(nside_out)+'nside-'+str(max_distance)+'max-'+str(distance_step)+'step.npz')
    npzfilename_new = npzfilename.with_name(f'{npzfilename.stem}_{event_num}{npzfilename.suffix}')
    npzfilename.rename(npzfilename_new)

    new_npz_path = os.path.join(save_folder, npzfilename_new)
    shutil.move(npzfilename_new, new_npz_path)

if __name__ == "__main__":

    global mapdir
    global nside_out
    global max_distance
    global num_distance_bins
    global distance_step
    global save_folder

    parser = argparse.ArgumentParser(description='Creates a 3D prob distribution for each fits file in a directory in the form of an npz file')
    parser.add_argument('--mapdir', help='Directory of the flattened event maps in fits format.', required=True)
    parser.add_argument('--nside_out', help='Nside value for the (reduced) cumulative 3D map. Default: 64.', required=True, type=int)
    parser.add_argument('--max_distance', help='Max distance for the 3D map (in Mpc). Default: 2000 Mpc.', required=True, type=int)
    parser.add_argument('--num_distance_bins', help='Number of distance bins (linspace style). Default: 11.', required=True, type=int)
    parser.add_argument('--save_folder', help="Name of save folder for npz files.", required=True)

    args = parser.parse_args()

    mapdir = args.mapdir
    nside_out = args.nside_out
    max_distance = args.max_distance
    num_distance_bins = args.num_distance_bins
    distance_step = int(max_distance / (num_distance_bins - 1))
    save_folder = args.save_folder

    files = glob(f'{mapdir}/*.fits')

    os.makedirs(save_folder, exist_ok=True)

    for filename in files:
        create_npz(filename, mapdir)
