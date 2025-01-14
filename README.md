# SaptACT
Scripts to perform SAPT calculations and generate input files for the Alexandria Chemistry Toolkit

You need to prepare a selection file that looks like
```
water#sodium
sodium#chloride
```
and so on. You also need to make sure that either the monomers are in the xyz/monomers directory
(in correct xyz format), or you prepare user structures of dimers that you put in the "user" directory.

Then, you can run the script
```
./run_calcs.py -v -ncore 32 -hours 48 -dimers na.dat -ndist 3 -mindist 0.95 -maxdist 1.05 -norient 4 
```
