#!/bin/bash

# Baysor installation:
# conda install -c conda-forge julia
# julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.build()'

if [ -n "$5" ]; then
    export JULIA_NUM_THREADS=$5
else
    export JULIA_NUM_THREADS=8
fi

MOLSPATH=$1
SCALE=$2
MINMOLS=$3
OUTPATH=$4

echo "Running Baysor on $MOLSPATH with scale $SCALE and min molecules per cell $MINMOLS"

/home/hartmana/miniconda3/share/julia/bin/baysor run $MOLSPATH --scale $SCALE --min-molecules-per-cell $MINMOLS --output $OUTPATH \
    --count-matrix-format tsv --save-polygons GeoJSON --plot
