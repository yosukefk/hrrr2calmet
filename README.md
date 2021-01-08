# hrrr2calmet
ingests HRRR datasset into CALMET/CALPUFF model

## install

miniconda and conda-forge recommended

If you are completely new to Anaconda/Python, I recommend install "Miniconda". That should provide you either (1) dedicated way to start command line terminal with conda enabled, or (2) your default setting is changed to use with conda or (3) you have to always type conda activate before using conda based tools. In either way, command line prompt will be change to have (base) PS C:\Users\Yosuke> or (base) ~ $ something like that, having "(base)" included in the prompt. You are ready to install.

Grab code from git repo

git clone https://github.com/yosukefk/hrrr2calmet.git

Grab required pacages

conda install -c conda-forge --file requirements.txt

## usage

`python3 hrrr2calmet.py input.grib2 output.m3d i0 j0 i1 j0`

where

input.grib2:  hrrr dataset in grib2 format
output.m3d:  "3D.DAT" format text file to be used in CALMET
i0: 1-based easting index of ll corner of extracted domain
j0: 1-based northin index of ll corner of extracted domain
i1: 1-based easting index of ur corner of extracted domain
j1: 1-based northin index of ur corner of extracted domain

## disclaimer

don't sue me.

currently, indices (i0, j0, i1, j1) i assume you somehow look it up, knowing where the domain is.  shouldn't be that hard to do this by specifying corner lat/lon should be fairly easy, as i now projection (using `rasterio` to read that from grib2 file), but not implemented yet.

right now one file with one time step (hourly analysis file, for example) got translated into one 3D.DAT format file.  i may add capability to process series of grib2 files into one 3D.DAT file.  But this can be done flexibily by CALMET, so i didnt do ot.
