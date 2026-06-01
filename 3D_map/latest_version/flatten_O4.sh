#!/bin/bash

# Usage:
# ./flatten_all.sh --in INPUT_DIR --out OUTPUT_DIR [--nside 64]

INDIR=""
OUTDIR=""
NSIDE=64  #Default

while [[ $# -gt 0 ]]; do
    case $1 in
        --in)
            INDIR="$2"
            shift 2
            ;;
        --out)
            OUTDIR="$2"
            shift 2
            ;;
        --nside)
            NSIDE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$INDIR" || -z "$OUTDIR" ]]; then
    echo "Usage: $0 --in INPUT_DIR --out OUTPUT_DIR [--nside NSIDE]"
    exit 1
fi

mkdir -p "$OUTDIR"

for file in "$INDIR"/*.fits.gz; do

    base=$(basename "$file")
    outfile="$OUTDIR/${base%_flat.fits}"

    echo "Flattening $file -> $outfile (nside=$NSIDE)"

    ligo-skymap-flatten --nside "$NSIDE" "$file" "$outfile"
done
