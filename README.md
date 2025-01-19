# SaptACT
Scripts to perform SAPT calculations and generate input files for the Alexandria Chemistry Toolkit

## Prerequisites
You need to have a recent and working Psi4 installation, which you can obtain, for instance, using
[conda](https://anaconda.org/conda-forge/psi4). For conda to work, you need to install 
[miniconda](https://docs.anaconda.com/miniconda/).

## Performing calculations
You need to prepare a selection file that looks like
```
guanine#guanine
guanine#cytosine
guanine#water
cytosine#cytosine
cytosine#water
water#water
```
and so on. You also need to make sure that either the monomers are in the xyz/monomers directory
(in correct xyz format). 

Then, you can run the script
```
./run_calcs.py -v -ncore 32 -hours 48 -dimers sample.dat -ndist 3 -mindist 0.95 -maxdist 1.05 -norient 4 
```
This will create a subdirectory method-basisset/dimer-scans in which you find the results sorted by dimer name.
Note that you can change the SAPT level of theory by using the -method and -basis flags.


Alternatively, you prepare your own structures of dimers that you put in a "user" directory.
In that case the filenames should correspond to monomer1#monomer2xxx.pdb where xxx is a number. You have to use a slightly different command line to do the SAPT calculations:
```
./run_calcs.py -v -ncore 32 -hours 48 -dimers sample.dat -user 
```

Some more information is available directly from the script, by
```
./run_calcs.py -h
```

## Extracting inputs for ACT
Two steps are needed to extract the information from the Psi4 log files. First, run the following command:
```
nohup ./generate_json.py &
```
which will read the Psi4 log file and summarize the information in a small json file for each calculation. Debugging output will be stored in the nohup.out file.
Note that you can select the level of theory to be used and that the default is the same for all scripts.
The above script uses the multiprocessing Python library to speed up the reading of log files, and it will skip directories where there is a results.json file already. You can however use the -f flag to force the script to redo everything. See
```
./generate_json.py -h
```
for more information.

