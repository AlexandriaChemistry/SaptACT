# SaptACT
Scripts to perform SAPT calculations and generate input files for the Alexandria Chemistry Toolkit

## Prerequisites
You need to have a recent and working Psi4 installation, which you can obtain, for instance, using
(conda)[https://anaconda.org/conda-forge/psi4]. For conda to work, you need to install 
(miniconda)[https://docs.anaconda.com/miniconda/].

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
(in correct xyz format), or you prepare user structures of dimers that you put in the "user" directory.

Then, you can run the script
```
./run_calcs.py -v -ncore 32 -hours 48 -dimers sample.dat -ndist 3 -mindist 0.95 -maxdist 1.05 -norient 4 
```
This will create a subdirectory method-basisset/dimer-scans in which you find the results sorted by
dimer name.

## Extracting inputs for ACT
Will be added soon.
