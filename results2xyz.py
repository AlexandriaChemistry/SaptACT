#!/usr/bin/env python3

import glob, os, json, sys

def read_json(dataFile:str)->dict:
    if not os.path.exists(dataFile):
        print("Cannot find %s/%s" % ( os.getcwd(), dataFile ))
        return {}
    with open(dataFile, "r") as data:
        myjson = json.load(data)
        return myjson

def write_xyz(filenm:str, dimer:dict)->int:
    if not "mols" in dimer:
        print("No molecules in %s" % os.getcwd())
        return 0
    with open(filenm, "w") as outf:
        natom = 0
        for mol in range(len(dimer["mols"])):
            natom += len(dimer["mols"][mol]["atoms"])
        outf.write("%5d\n" % natom)
        energies = "energies"
        iener    = "InteractionEnergy"
        if energies in dimer and iener in dimer[energies]:
            outf.write(" %s Energy %g\n" % ( filenm, dimer[energies][iener] ) )
        else:
            outf.write(" %s\n" % filenm)
        for mol in range(len(dimer["mols"])):
            for nat in range(len(dimer["mols"][mol]["atoms"])):
                outf.write(" %3s  %12f  %12f  %12f\n" %
                           (dimer["mols"][mol]["atoms"][nat]["elem"],
                            dimer["mols"][mol]["atoms"][nat]["coords"][0],
                            dimer["mols"][mol]["atoms"][nat]["coords"][1],
                            dimer["mols"][mol]["atoms"][nat]["coords"][2]))
    return 1

if __name__ == "__main__":
    nxyz = 0
    with open("clean.sh", "w") as outf:
#        os.chdir("sapt2+(ccd)dmp2-aug-cc-pvtz/dimer-scans")
        os.chdir("sapt2+-aug-cc-pvdz/dimer-scans")
#        os.chdir("MP2-aug-cc-pvtz/monomer-sp")
        #os.chdir("MP2-aug-cc-pvtz/monomer-opt")
        for mydir in glob.glob("*"):
            if not os.path.isdir(mydir):
                continue
            os.chdir(mydir)
            for index in glob.glob("[0-9]*"):
                if not os.path.isdir(index):
                    continue
                os.chdir(index)
                myjson = read_json("results.json")
                if len(myjson.keys()) > 0:
                    nxyz += write_xyz(index+".xyz", myjson)
                else:
                    outf.write("cd '%s/..'; git rm -rf --ignore-unmatch %s\n" % ( os.getcwd(), index )) 
                os.chdir("..")
            os.chdir("..")
    print("Generated %d xyz files for SAPT calcs" % nxyz)
