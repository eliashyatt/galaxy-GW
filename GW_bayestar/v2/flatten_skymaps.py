import subprocess
import glob
import os
import shutil
import argparse


def flatten(mapdir, nside):

    #files = glob.glob(os.path.join("./skymaps/1000_plus/1000-2000/", "*.fits"))
    os.chdir(mapdir)

    os.makedirs('./flattened_maps', exist_ok=True)

    files = glob.glob("*.fits")

    print(f"Found {len(files)} files")

    for f in files:
        out = f.replace(".fits", "_flat.fits")
        subprocess.run(["ligo-skymap-flatten", "--nside", nside, f, out])
        dest = f'./flattened_maps/{out}'
        shutil.move(out, dest)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Flattens skymaps in a directory')
    parser.add_argument('--mapdir', help='Directory of skymaps in the form of fits files', required=True)
    parser.add_argument('--nside', help='nside for flattened maps', required=True)
    args = parser.parse_args()

    mapdir = args.mapdir
    nside = args.nside

    print('Starting...')

    flatten(mapdir, nside)

    print('Done!')