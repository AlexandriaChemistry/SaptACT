#!/usr/bin/env python3

import os, sys, glob, argparse, shutil

import molprops
import get_mol_dict as gmd
from atomic_heat_of_formation import *
from run_calcs  import special_basis, get_monoq, get_dimer_selection, get_monomer_selection, get_diatomics
from elements import *
from psi4files  import *
from mol_csv_api import *

write_xml_debug   = False

def add_one_diatomic(molname:str, diatomic:dict,
                     Mols:Molecules,
                     method:str, basis:str, deltaEmax:float, forceMax:float, 
                     rmax_factor:float, ahof) -> molprops.Molprop:
    scandir = "diatomic-scans"
    mydir = ("%s-%s/%s/%s" % ( method, basis, scandir, molname ))
    if not os.path.exists(mydir):
        #print("No such directory %s" % mydir)
        return None
    olddir = os.getcwd()
    os.chdir(mydir)
    ener = {}
    rmin = None
    emin = 1000000
    exvg = "energy.xvg"
    if not os.path.exists(exvg):
        print("No such file %s/%s" % ( os.getcwd(), exvg ))
        os.chdir(olddir)
        return None
    with open(exvg, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if len(words) == 2:
                r = float(words[0])
                ener[words[0]] = float(words[1])
                if ener[words[0]] < emin:
                    rmin = r
                    emin = float(words[1])
    force = {}
    if os.path.exists("force.xvg"):
        with open("force.xvg", "r") as inf:
            for line in inf:
                words = line.strip().split()
                if len(words) == 2:
                    r = float(words[0])
                    force[words[0]] = float(words[1])

    if not rmin:
        print("Something wrong with files in %s, no minimum found" % os.getcwd())
        os.chdir(olddir)
        return None
    elements = [ diatomic["ai"], diatomic["aj"] ]
    coords   = [ [ 0, 0, 0 ], [ 0, 0, rmin ] ]
    symmetry = 2
    mol      = Mols.find_mol(molname)
    if mol:
        symmetry = mol.symmetry_number

    mp1  = molprops.Molprop(molname)
    MD   = gmd.MoleculeDict()
    qtot = diatomic["chargei"] + diatomic["chargej"]
    if not MD.from_coords_elements(elements, coords, qtot):
        print("Cannot analyze coordinates")
        os.chdir(olddir)
        return None
    if diatomic["covalent"]:
        frag = molprops.Fragment(MD.inchi, diatomic["charge"], diatomic["mult"],
                                 symmetry, MD.atoms, MD.mol_weight, MD.formula)
        mp1.add_fragment(frag)
        for b in MD.bonds:
            mp1.add_bond(b[0], b[1], MD.bonds[b])
    else:
        MD1   = gmd.MoleculeDict()
        if not MD1.from_coords_elements(elements[:1], coords[:1], diatomic["chargei"]):
            print("Cannot analyze coordinates")
            os.chdir(olddir)
            return None
        mp1.add_fragment(molprops.Fragment(MD1.inchi, diatomic["chargei"], diatomic["mult"],
                                           symmetry, [ 1 ], 
                                           atomprops[elements[0]]["mass"],
                                           MD1.formula))
        MD2   = gmd.MoleculeDict()
        if not MD2.from_coords_elements(elements[1:], coords[1:], diatomic["chargej"]):
            print("Cannot analyze coordinates")
            os.chdir(olddir)
            return None
        mp1.add_fragment(molprops.Fragment(MD2.inchi, diatomic["chargej"], diatomic["mult"],
                                           symmetry, [ 2 ],
                                           atomprops[elements[1]]["mass"],
                                           MD2.formula))
    temp    = 0.0
    lot     = ("%s/%s" % (method, basis))
    lots    = []
    for elem in elements:
        if basis in special_basis and elem in special_basis[basis]:
            lots.append("%s/%s" % ( method, special_basis[basis][elem]) )
        else:
            lots.append("%s/%s" % ( method, basis) )
    charges = [ 0, 0 ]
    for my_r in ener:
        jobtype = "SP"
        if float(my_r) == rmin:
            jobtype = "Opt"
        # We only use energies within a certain range from the minimum
        fz = 0
        if my_r in force:
            fz = float(force[my_r])
        if (ener[my_r] - emin <= deltaEmax and 
            float(my_r) <= float(rmin)*rmax_factor and abs(fz) <= forceMax):
            myexp = molprops.Experiment("Theory", "Spoel2023a", "Psi4", method, basis,
                                        "conformation", jobtype, "none", True)
            deltaE0 = compute_dhform(ener[my_r], elements, ahof,
                                     lots, charges, temp)
            myexp.add_energy("DeltaE0", "Hartree", 0, "gas", deltaE0)
            for atom in range(2):
                obtype = MD.atoms[atom]["obtype"]
                myexp.add_atom(elements[atom], obtype, atom+1,
                               "Angstrom", 0, 0, atom*float(my_r), "Hartree/Bohr",
                               0, 0, fz)
                fz = -fz
            mp1.add_experiment(myexp)
    os.chdir(olddir)
    return mp1

def read_opt(comp:str, mp:molprops.Molprop):
    print("Do not know how to read opt file for %s. Forget about it." % comp)

def read_esp(comp:str)->list:
    grid_file  = "grid.dat"
    esp_data   = "grid_esp.dat"
    pot = []
    if os.path.exists(grid_file) and os.path.exists(esp_data):
        grid = []
        with open(grid_file, "r") as inf:
            for line in inf:
                words = line.strip().split()
                if len(words) == 3:
                    try:
                        grid.append([ float(words[0]), float(words[1]), float(words[2]) ])
                    except ValueError:
                        print("Strange line '%s' in %s/%s" % ( line.strip(), comp, grid_file ))
        with open(esp_data, "r") as inf:
            index = 0
            for line in inf:
                try:
                    val = float(line.strip())
                    if index < len(grid):
                        pot.append({ "grid": grid[index], "value": val})
                    else:
                        print("Inconsistency: %s has more lines than %s for %s" % ( esp_data, grid_file, comp ) )
                        return []
                except ValueError:
                    print("Strange line '%s' in %s/%s" % ( line.strip(), comp, grid_file ))
                index += 1
    return pot

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

def add_dimers(mm:molprops.Molprop, Mol:Molecules, args:list, logf, ahof):
    lot = args.method + "-" + args.basis
    if not os.path.exists(lot):
        logf.write("No such LoT %s\n" % lot)
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
    for calc in [ "opt", "scans" ]:
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
                    potential = []
                    status = outfile.read(args, molnames, monoq, calc, mycomplex, index,
                                          Mols, ahof, temperature, potential,logf)
                    if not status == Psi4Error.OK:
                        logf.write("%s:  %s\n" % (os.getcwd(), psi4msg(status) ) )
                    os.chdir("..")
            mm.add_molecule(outfile.molprop(), True)
            os.chdir("..")
        os.chdir("..")
    os.chdir("..")

def add_monomers(mm:molprops.Molprop, Mol:Molecules, args:list, logf, ahof):
    lot = args.method + "-" + args.basis
    if not os.path.exists(lot):
        logf.write("No such LoT %s\n" % lot)
        return
    monoq       = get_monoq()
    monomerlist = []
    if None != args.selection:
        monomerlist = get_monomer_selection(args.selection)
    temperature = 0
    ahof        = AtomicHOF(None, temperature, False)
    os.chdir(lot)

    for calc in [ "esp", "opt", "scans", "sp" ]:
        monomers = "monomer-" + calc
        if not os.path.exists(monomers):
            continue
        os.chdir(monomers)
        logf.write("Will try to add monomers for %s jobtype %s\n" % (lot, calc))
        for mymol in glob.glob("*"):
            if not os.path.isdir(mymol):
                continue
            myskip = len(monomerlist) != 0
            for km in monomerlist:
                if mymol == km["mon1"]:
                    myskip = False
            if myskip:
                continue
            os.chdir(mymol)
            # Find compounds
            # First time around we need to add fragments
            outfile = Psi4Files(mymol, [ mymol ], args.json)
            for mysubdir in glob.glob("*"):
                if not os.path.isdir(mysubdir):
                    continue
                os.chdir(mysubdir)
                potential = []
                if "esp" == calc and not args.skipESP:
                    potential = read_esp(mymol)
                if write_xml_debug:
                    logf.write("Will try to read data from %s\n" % mymol)
                status = outfile.read(args, [ mymol ], monoq, calc, mymol, mysubdir,
                                      Mols, ahof, temperature, potential, logf)
                if Psi4Error.OK != status:
                    logf.write("%s in %s\n" % ( psi4msg(status), os.getcwd() ))
                os.chdir("..")
            mm.add_molecule(outfile.molprop(), True)
            os.chdir("..")
        os.chdir("..")
    os.chdir("..")

def parse_args():
    desc = "Extract data from calculations and store them in ACT molprop xml files."
    parser  = argparse.ArgumentParser(description=desc)
    defbasis = "aug-cc-pvtz"
    parser.add_argument("-basis", "--basis", help="Basis set, default is "+defbasis, type=str, default=defbasis)
    defmethod = "MP2"
    parser.add_argument("-method","--method", help="QM method default "+defmethod, type=str, default=defmethod)
    defname = "molprop.xml"
    deltaEmax = 0.02
    parser.add_argument("-dEmax", "--deltaEmax", help="Highest energy above the minimum to include. For SAPT calcs the Exchange energy will be used. If one or both of the compounds in a dimer is charged, the value provided here will be multiplied by (1+sum q) to allow for higher energies. Default "+str(deltaEmax)+" Hartree", type=float, default=deltaEmax)
    fMax = 1
    parser.add_argument("-fmax", "--forceMax", help="Highest absolute force (Hartree/Bohr) to include, default "+str(fMax), type=float, default=fMax)
    rmax = 6
    parser.add_argument("-rmax", "--rmax", help="Largest relative distance to include, obtained by multiplying the distance at which a minimum is found by this number, applies to diatomic scans and dimer scans, but in that case absolute distance in Angstrom. Default "+str(rmax), type=float, default=rmax)
    rmin = 1.5
    parser.add_argument("-rmin", "--rmin", help="Shortest absolute distance in a dimer scan, distance in Angstrom. Default "+str(rmin), type=float, default=rmin)
    parser.add_argument("-sel", "--selection", help="Extract monomer and dimers based on compounds in a selection file, please provide file name with this flag", type=str, default=None)
    parser.add_argument("-dimsel", "--dimerselection", help="Extract dimer interactions based on particular calculations of compounds in a file, according to 'dimer/0xxx'. Please provide file name with this flag", type=str, default=None)
    parser.add_argument("-o", "--output", help="Name of the output molprop file, default "+defname, type=str, default=defname)
    parser.add_argument("-v", "--verbose", help="Write debugging output", action="store_true")
    parser.add_argument("-deltaHF", "--deltaHF", help="Move the HF and MP2 correction from Induction to a specific InductionCorrrection", action="store_true")
    parser.add_argument("-sf", "--skipforces", help="Do not write forces, e.g. for SAPT methods", action="store_true")
    parser.add_argument("-se", "--skipESP", help="Do not write ESP data to make smaller files", action="store_true")
    parser.add_argument("-json", "--json", help="Read json files only", action="store_true")
    logfn = "write_molprop.log"
    parser.add_argument("-g", "--logfn", help="Debugging output file, default "+logfn, type=str, default=logfn)
    ref = "Spoel2026a"
    parser.add_argument("-ref", "--reference", help="Reference to insert in molprop file, default  "+ref, type=str, default=ref)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    if args.verbose:
        write_xml_debug = True
        molprops.set_molprops_debug(True)
    tempxml = "temptemp.xml"
    with open(args.logfn, "w") as logf:
        mm   = molprops.Molprops()
        mm.open(tempxml)
        get_atomprops()
        Mols = Molecules()
        Mols.read_default()
        ahof = AtomicHOF(None, 0, False)
        if False:
            diatomics = get_diatomics()
            for dim in diatomics.keys():
                mp = add_one_diatomic(dim, diatomics[dim], Mols,
                                      args.method, args.basis, args.deltaEmax, args.forceMax, args.rmax, ahof)
                if mp:
                    mm.add_molecule(mp, True)
        add_monomers(mm, Mols, args, logf, ahof)
        add_dimers(mm, Mols, args, logf, ahof)
        mm.close()
        alex = "alexandria"
        if shutil.which(alex):
            os.system("%s edit_mp -mp %s -o %s -v 3" % ( alex, tempxml, args.output ) )
        else:
            print("Cannot find %s executable" % alex)
            os.system("mv %s %s" % ( tempxml, args.output ) )
