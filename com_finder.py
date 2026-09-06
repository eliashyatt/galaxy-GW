import numpy as np
import pandas as pd
import healpy as hp
import argparse
from glob import glob
import os


def find_com(filename, mapdir):

    npz_data = np.load(filename)["cumulative_reduced3Dmap"]
    #npz_data = np.load('./3D_map/latest_version/cumulative_reduced3Dmap-64nside-2000max-100step_0014.npz')["cumulative_reduced3Dmap"]

    npz_data = npz_data[:,1:]

    basename = os.path.basename(filename)
    splitname = basename.split('-')
    distance_step = int(splitname[-1][:-14])

    flattened_data = npz_data.T.flatten()

    npix = npz_data.shape[0]
    ndist = npz_data.shape[1]

    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, np.arange(npix), nest=True)

    distances = np.repeat(np.arange(ndist) * distance_step, npix)

    probs = (distances**2) * distance_step * flattened_data

    probs = probs / probs.sum()

    theta = np.tile(theta, ndist)
    phi = np.tile(phi, ndist)

    x = np.sum(distances * np.sin(theta) * np.cos(phi) * probs)
    y = np.sum(distances * np.sin(theta) * np.sin(phi) * probs)
    z = np.sum(distances * np.cos(theta) * probs)

    r = np.sqrt(x**2 + y**2 + z**2)

    ra = np.arctan2(y, x)
    if ra < 0:
        ra += 2*np.pi

    dec = np.arcsin(z / r)

    com_vec = {"ra": ra, "dec": dec, "distance": r}

    return com_vec

    #df = pd.DataFrame([com_vec])
    #df.to_csv("com_position_0014.csv", index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Finds the center of mass points of 3D prob dist of skymaps in the form of npz files')
    parser.add_argument('--mapdir', help='Directory of npz files', required=True)

    args = parser.parse_args()

    mapdir = args.mapdir
    folder_name = mapdir.split('/')[-1]
    print(folder_name)

    #/home/vboxuser/Desktop/Code_for_Thesis/galaxy-GW/GW_bayestar/v2/skymaps/1000_fits

    files = glob(f'{mapdir}/*.npz')

    print("Finding COM for", len(files), "files")

    coms = []
    events = []

    for filename in files:
        com_vec = find_com(filename, mapdir)
        coms.append(com_vec)
        event = filename.split('_')[-1][:-4]
        #print(event)
        events.append(event)

    df = pd.DataFrame(coms)
    df['event'] = events
    df.to_csv(f"com_positions_{folder_name}.csv", index=False)

    print('Done!')


