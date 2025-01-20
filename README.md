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
The above command line will generate 4 sets of random rotations of the two molecules and place them at three distances along
the vector connecting the nearest two atoms in the dimer. The distances are determined from the sum of the Van der Waals radii
of the two closest atoms and then scaled by the factors indicated, in this case 0.95, 1.0, 1.05. The command should give you
3 x 4 calculations per dimer. Note that you can change the SAPT level of theory by using the -method and -basis flags.


Alternatively, you prepare your own structures of dimers that you put in a "user" directory.
In that case the filenames should correspond to monomer1#monomer2xxx.pdb where xxx is a number. 
The selection file remains the same, but you have to use a slightly different command line to do the SAPT calculations:
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
which will read the Psi4 log file and summarize the information in a small json file for each calculation. Debugging output will be stored in the nohup.out file and file called *generate_json.log*. Please check the content of these files after running the command.
Note that you can select the level of theory to be used and that the default is the same for all scripts.
The above script uses the multiprocessing Python library to speed up the reading of log files, and it will skip directories where there is a results.json file already. You can however use the -f flag to force the script to redo everything. See
```
./generate_json.py -h
```
for more information.

The second step is to store the information from the json files to the *molprop* format that the Alexandria Chemistry Toolkit can use.
To create a molprop file, run:
```
./write_molprop.py -o new.xml
```
The script will create a log file with diagnostics as well.
