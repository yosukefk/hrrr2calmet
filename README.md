# hrrr2calmet
Ingests [High-Resolution Rapid Refresh (HRRR)](https://rapidrefresh.noaa.gov/hrrr/) meteorological model dataset into [CALMET/CALPUFF dispersion model](http://www.src.com/).

## Install

Use of Miniconda and conda-forge recommended

If you are completely new to Anaconda/Python, I recommend install "Miniconda". That should provide you either (1) dedicated way to start command line terminal with conda enabled, 
or (2) your default setting is changed to use with conda or (3) you have to always type conda activate before using conda based tools. 
In either way, command line prompt will be change to have "(base) PS C:\Users\Yosuke> " or "(base) ~ $ " something like that, having "(base)" included in the prompt. 


Grab code from git repo

`git clone https://github.com/yosukefk/hrrr2calmet.git`

Grab required packages

`conda install -c conda-forge --file requirements.txt`

## Usage

`python3 hrrr2calmet.py outfile i0 j0 i1 j0 infile [...]`

where `outfile` is output "3D.DAT" format text file for CALMET, `i0`, `j0`, `i1` and `j1` are 1-based indices of lower-left and upper-right grid cells (i for easting, j for northing).  `infile [...]` are series of input grib2 format files.

## Disclaimer

Please use the tool at your own risk.

Currently, indices (i0, j0, i1, j1) I assume you can somehow look these up, knowing where the domain of your interest is.  

