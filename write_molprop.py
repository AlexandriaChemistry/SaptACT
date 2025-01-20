#!/usr/bin/env python3

import os, sys, glob, argparse

import openbabel as ob
import molprops
from gaff_to_alexandria import *
from run_calcs  import special_basis, get_monoq, get_dimer_selection
from elements import *
from psi4files  import *

write_xml_debug   = False

def get_dimerselection(selfile:str)->list:
    # Will return a list of dimer data including specific output files
    dimsel = {}
    with open(selfile, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if 2 == len(words):
                if not words[0] in dimsel:
                    dimsel[words[0]] = []
                dimsel[words[0]].append(words[1])
    # Write temporary file
    mytmp = "tmp123456789.dat"
    with open(mytmp, "w") as outf:
        for dim in dimsel.keys():
            outf.write("%s\n" % dim)
    dims = get_dimer_selection(mytmp)
    os.unlink(mytmp)
    for d in dims:
        if not d["pair"] in dimsel:
            sys.exit("Cannot find dimer %s in selection. What's up with that?" % d["pair"])
        d["selection"] = dimsel[d["pair"]]
    return dims

def add_dimers(mm:molprops.Molprop, Mol:Molecules, args:list, logf):
    lot = args.method + "-" + args.basis
    if not os.path.exists(lot):
        logf.write("No directory corresponding to Level of Theory %s\n" % lot)
        return
    monoq       = get_monoq()
    temperature = 0
    dimerlist   = []
    if None != args.selection:
        if None != args.dimerselection:
            sys.exit("Use only one of the options -selection and -dimsel")
        dimerlist = get_dimer_selection(args.selection)
    elif None != args.dimerselection:
        dimerlist = get_dimerselection(args.dimerselection)
    os.chdir(lot)
    for calc in [ "scans" ]:
        dimers = "dimer-" + calc
        if not os.path.exists(dimers):
            continue
        os.chdir(dimers)
        logf.write("Will try to add dimers for %s\n" % lot)
        for mycomplex in glob.glob("*"):
            if mycomplex == "fluoride#bromide":
                logf.write("Found %s\n" % mycomplex)
            if not os.path.isdir(mycomplex):
                continue
            myskip = len(dimerlist) != 0
            for dk in dimerlist:
                if mycomplex == dk["pair"]:
                    myskip = False
            if myskip:
                continue
            # Find compounds
            molnames      = mycomplex.split("#")
            if len(molnames) != 2:
                logf.write("Skipping incomprehensible complex %s\n" % mycomplex)
                continue
            os.chdir(mycomplex)
            outfile = Psi4Files(mycomplex, molnames, args.json)
            # Two compounds in a complex
            # First time around we need to add fragments
            if "opt" == calc:
                read_opt(mycomplex, outfile.molprop())
            else:
                outfile.find_energy_minimum()
                if write_xml_debug:
                    logf.write("Energy minimum for %s = %g\n" % ( mycomplex, outfile.edimerMin ))
                for index in glob.glob("*"):
                    if not os.path.isdir(index):
                        continue
                    # If there is a user selection with -dimsel, honour that.
                    if dimerlist:
                        selkey = "selection"
                        found = False
                        for ddd in dimerlist:
                            if ddd["pair"] == mycomplex:
                                if (selkey in ddd and index in ddd[selkey]) or not selkey in ddd:
                                    found = True
                        if not found:
                            continue
                    os.chdir(index)
                    status = outfile.read(args, g2a, molnames, monoq, calc, mycomplex, index,
                                          Mols, temperature, logf)
                    if not status == Psi4Error.OK:
                        logf.write("%s:  %s\n" % (os.getcwd(), psi4msg(status) ) )
                    os.chdir("..")
            mm.add_molecule(outfile.molprop(), True)
            os.chdir("..")
        os.chdir("..")
    os.chdir("..")

def parse_args():
    desc = "Extract data from calculations and store them in ACT molprop xml files."
    parser  = argparse.ArgumentParser(description=desc)
    defbasis = "aug-cc-pvdz"
    parser.add_argument("-basis", "--basis", help="Basis set, default is "+defbasis, type=str, default=defbasis)
    defmethod = "sap2+"
    parser.add_argument("-method","--method", help="QM method default "+defmethod, type=str, default=defmethod)
    defname = "molprop.xml"
    deltaEmax = 0.02
    parser.add_argument("-dEmax", "--deltaEmax", help="Highest energy above the minimum to include. For SAPT calcs the Exchange energy will be used. If one or both of the compounds in a dimer is charged, the value provided here will be multiplied by (1+sum q) to allow for higher energies. Default "+str(deltaEmax)+" Hartree", type=float, default=deltaEmax)
    rmax = 6
    parser.add_argument("-rmax", "--rmax", help="Largest relative distance to include, in Angstrom. Default", type=float, default=rmax)
    parser.add_argument("-sel", "--selection", help="Extract dimers based on compounds in a selection file, please provide file name with this flag", type=str, default=None)
    parser.add_argument("-dimsel", "--dimerselection", help="Extract dimer interactions based on particular calculations of compounds in a file, according to 'dimer/0xxx'. Please provide file name with this flag", type=str, default=None)
    parser.add_argument("-o", "--output", help="Name of the output molprop file, default "+defname, type=str, default=defname)
    parser.add_argument("-v", "--verbose", help="Write debugging output", action="store_true")
    parser.add_argument("-deltaHF", "--deltaHF", help="Move the HF and MP2 correction from Induction to a specific InductionCorrrection", action="store_true")
    parser.add_argument("-json", "--json", help="Read json files only", action="store_true")
    logfn = "write_molprop.log"
    parser.add_argument("-g", "--logfn", help="Debugging output file, default "+logfn, type=str, default=logfn)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    if args.verbose:
        write_xml_debug = True
        molprops.set_molprops_debug(True)
    with open(args.logfn, "w") as logf:
        mm   = molprops.Molprops()
        mm.open(args.output)
        g2a       = GaffToAlexandria()
        get_atomprops()
        Mols = Molecules()
        Mols.read_default()
        add_dimers(mm, Mols, args, logf)
        mm.close()

    print("Please check %s for diagnostic information and the output %s" % ( args.logfn, args.output ) )
