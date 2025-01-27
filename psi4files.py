#!/usr/bin/env python3

import glob, json, os, sys, math
import openbabel as ob
import molprops
import get_mol_dict       as gmd
from gaff_to_alexandria       import *
from run_calcs  import special_basis
from elements  import *
from enum import Enum
from mol_csv_api import *

class Psi4Error(Enum):
    OK = 0
    IncompleteJson = 1
    MissingJson = 2
    NoMolnames = 3
    NotADimer = 4
    Coordinates = 5
    HighEnergy = 6
    LargeDistance = 7
    MonomerCharge = 8
    InconsistentJson = 9

def psi4msg(p4:Psi4Error):
    msgs = {
        Psi4Error.OK: "",
        Psi4Error.IncompleteJson: "Psi4Error: Incomplete json file",
        Psi4Error.InconsistentJson: "Psi4Error: Inconsistent json file",
        Psi4Error.MissingJson: "Psi4Error: Missing json file",
        Psi4Error.NoMolnames: "Psi4Error: No molnames provided",
        Psi4Error.NotADimer: "Psi4Error: Not a dimer in the json file",
        Psi4Error.Coordinates: "Psi4Error: Coordinates incorrect",
        Psi4Error.HighEnergy: "Psi4Error: Energy too high",
        Psi4Error.LargeDistance: "Psi4Error: Distance too large",
        Psi4Error.MonomerCharge: "Psi4Error: Compound missing in monomer_charge.csv"
    }
    return msgs[p4]

def check_dist(allcoords:list)->float:
    if len(allcoords) != 2:
        print("allcoords has len %d" % ( len(allcoords) ) )
        return 0
    rmin = None
    for i in range(len(allcoords[0])):
        for j in range(len(allcoords[1])):
            r2 = 0
            for m in range(3):
                r2 += (allcoords[0][i][m]-allcoords[1][j][m])**2
            if not rmin or r2 < rmin:
                rmin = r2
    if rmin:
        return math.sqrt(rmin)
    else:
        return 0

verbose              = False
hasPrintedIncomplete = False
edimerMinDefault     = 1e8
class Psi4Files:
    def __init__(self, molpropname:str, molnames:list, json:bool):
        self.status    = Psi4Error.OK
        self.nmols     = len(molnames)
        if self.nmols == 0:
            self.status = Psi4Error.NoMolnames
            return
        self.MD        = []
        for i in range(self.nmols):
            self.MD.append(gmd.MoleculeDict())
        self.edimerMin = edimerMinDefault
        # Create new molecule
        self.mp1       = molprops.Molprop(molpropname)
        # Default SP, but in the best case this is part of the json
        self.jobtype   = "SP"
        self.json      = json
        self.firstDone = False

    def do_read_json(self, dataFile:str)->dict:
        with open(dataFile, "r") as data:
            myjson = json.load(data)
        jt = "jobtype"
        if jt in myjson:
            self.jobtype = myjson[jt]
        return myjson

    def read(self, args, g2a,
             molnames:list, monoq:dict, calc:str, mycomplex:str, filename:str,
             Mols:Molecules,
             temperature:float, logf)->Psi4Error:
        global hasPrintedIncomplete
        self.status = Psi4Error.OK
        nfrag = len(mycomplex)
        self.dataFile = "results.json"
        energies = "energies"
        mols     = "mols"
        atoms    = "atoms"
        if not os.path.exists(self.dataFile):
            self.status = Psi4Error.MissingJson
            return self.status
        else:
            mydict = self.do_read_json(self.dataFile)
            if (len(mydict.keys()) == 0 or 
                (energies in mydict and len(mydict[energies]) == 0) or
                not mols in mydict or
                (mols in mydict and len(mydict[mols]) == 0)):
                self.status = Psi4Error.IncompleteJson
            if (len(molnames) != self.nmols or 
                len(molnames) != len(mydict[mols])):
                self.status = Psi4Error.InconsistentJson
        if self.status != Psi4Error.OK:
            return self.status

        fAbsMax = 0
        self.lots      = []
        self.charges   = []
        # Create new experiment
        myexp = molprops.Experiment("Theory", "Spoel2025a", "Psi4", args.method, args.basis,
                                    "conformation", self.jobtype, filename, True)
        frag_charges = []
        sum_q_abs    = 0
        # TODO Fix charges!
        for j in range(self.nmols):
            if molnames[j] in monoq:
                mq = monoq[molnames[j]]["charge"]
                frag_charges.append(mq)
                sum_q_abs += abs(mq)
            else:
                self.status = Psi4Error.MonomerCharge
                return self.status

        # This is for all atoms
        elements  = []
        offset    = 0
        allcoords = []
        # Interpret data from json
        for k in range(len(mydict[mols])):
            coords   = []
            forces   = []
            myatoms  = []
            natom    = len(mydict[mols][k][atoms])
            for j in range(natom):
                elements.append(mydict[mols][k][atoms][j]["elem"])
                coords.append(mydict[mols][k][atoms][j]["coords"])
                if "forces" in mydict[mols][k][atoms][j]:
                    fj = mydict[mols][k][atoms][j]["forces"]
                    forces.append(fj)
                    fAbsMax = max(fAbsMax, math.sqrt(fj[0]**2+fj[1]**2+fj[2]**2))
                else:
                    forces.append([0,0,0])
                lot  = args.method + "/"
                if (args.basis in special_basis and
                    elements[offset+j] in special_basis[args.basis]):
                    lot += special_basis[args.basis][elements[offset+j]]
                else:
                    lot += args.basis
                self.lots.append(lot)
                # Atom counter for fragments
                myatoms.append(offset+j+1)
                self.charges.append(0)
            if k >= len(self.MD):
                sys.exit("k = %d len(self.MD) = %d" % ( k, len(self.MD)))
            if (self.status == Psi4Error.OK and
                not self.MD[k].from_coords_elements(elements[offset:], coords, charge=frag_charges[k])):
                logf.write("Cannot analyze coordinates using OpenBabel\n")
                self.status = Psi4Error.Coordinates
            if self.status == Psi4Error.OK and not self.firstDone:
                symmetry = 1
                mult     = 1
                mol      = Mols.find_mol(molnames[k])
                if mol:
                    symmetry = mol.symmetry_number
                    mult     = mol.mult
                frag     = molprops.Fragment(self.MD[k].inchi, frag_charges[k], mult, symmetry,
                                             myatoms, self.MD[k].mol_weight, self.MD[k].formula)
                self.mp1.add_fragment(frag)
                for b in self.MD[k].bonds:
                    self.mp1.add_bond(offset+b[0], offset+b[1], self.MD[k].bonds[b])
            if self.status == Psi4Error.OK:
                # Now add the atoms to the experiment
                if 0 == len(self.MD[k].atoms):
                    sys.exit("No atoms in compound %d" % k )
                for atom in range(natom):
                    obtype = g2a.rename(self.MD[k].atoms[1+atom]["obtype"])
                    # Atom numbers should match bonds
                    myexp.add_atom(elements[offset+atom], obtype, offset+atom+1,
                                   "Angstrom",
                                   coords[atom][0], coords[atom][1], coords[atom][2],
                                   "Hartree/Bohr",
                                   forces[atom][0], forces[atom][1], forces[atom][2])
                offset   += natom
                allcoords.append(coords)
        # This has to come after the atoms because we need the elements here
        if not energies in mydict:
            self.status = Psi4Error.IncompleteJson
            logf.write("No key %s in %s/%s\n" % ( energies, os.getcwd(), self.dataFile ) )
        else:
            if len(molnames) == 2:
                # We have a dimer
                if self.edimerMin == edimerMinDefault:
                    logf.write("edimerMin has not been set\n")
                    self.status = Psi4Error.IncompleteJson
                exch    = "Exchange"
                intener = "InteractionEnergy"

                if self.edimerMin == edimerMinDefault:
                    logf.write("In %s\n" % os.getcwd())
                    logf.write(f"Minimum energy not available for {molnames[0]}-{molnames[1]}. Considering exchange energy.\n")
                    
                    if exch in mydict[energies]:
                        if mydict[energies][exch] > args.deltaEmax * (1 + sum_q_abs):
                            logf.write(f"Skipping high {exch} energy {mydict[energies][exch]} for {molnames[0]}-{molnames[1]}\n")
                            self.status = Psi4Error.HighEnergy
                        else:
                            logf.write(f"Processing but ignoring exchange energy for {molnames[0]}-{molnames[1]}\n")
                    else:
                        logf.write(f"No exchange energy available for {molnames[0]}-{molnames[1]}\n")
                        self.status = Psi4Error.MissingEnergy

                if exch in mydict[energies]:
                    if mydict[energies][exch] > args.deltaEmax*(1+sum_q_abs):
                        logf.write("Skipping high %s energy %g for %s-%s\n" %
                                   ( exch, mydict[energies][exch], molnames[0], molnames[1] ))
                        self.status = Psi4Error.HighEnergy
                elif intener in mydict[energies]:
                    if (mydict[energies][intener] - self.edimerMin) > args.deltaEmax*(1+sum_q_abs):
                        logf.write("Skipping high %s energy %g for %s-%s\n" %
                                   ( intener, mydict[energies][intener], molnames[0], molnames[1] ))
                        self.status = Psi4Error.HighEnergy
                if args.deltaHF:
                    dhf = "delta HF,r (2)"
                    mp2 = "delta MP2,r (2)"
                    induc = "Induction"
                    if (dhf in mydict[energies] and
                        induc in mydict[energies]):
                        indcorr = mydict[energies][dhf]
                        if mp2 in mydict[energies]:
                            indcorr += mydict[energies][mp2]
                        mydict[energies][induc] -= indcorr
                        mydict[energies]["InductionCorrection"] = indcorr
                
                if self.status == Psi4Error.OK:
                    for energy in mydict[energies].keys():
                        myexp.add_energy(energy, "Hartree", 0, "gas", mydict[energies][energy])
                else:
                    logf.write("Skipping energies in %s because of %s" %
                               ( os.getcwd(), psi4msg(self.status) ) )
            else:
                self.status = Psi4Error.InconsistentJson
            
        
        if self.status == Psi4Error.OK:
            if self.nmols == 2:
                mydist = check_dist(allcoords)
                if mydist > args.rmax:
                    self.status = Psi4Error.LargeDistance
                    logf.write("Distance (%g) too large (max %g) between compounds in %s\n" %
                               ( mydist, args.rmax, os.getcwd() ) )
        if self.status == Psi4Error.OK:
            self.mp1.add_experiment(myexp)
        self.firstDone = True
        return self.status

    def molprop(self)->molprops.Molprop:
        return self.mp1

    def status(self)->Psi4Error:
        return self.status

    def find_energy_minimum(self):
        # First extract the lowest energy for this complex
        for index in glob.glob("*"):
            if not os.path.isdir(index):
                continue
            os.chdir(index)
            rj = "results.json"
            if os.path.exists(rj):
                mydict = self.do_read_json(rj)
                ee = "energies"
                ie = "InteractionEnergy"
                if ee in mydict and ie in mydict[ee]:
                    ener = mydict[ee][ie]
                    if ener < self.edimerMin:
                        self.edimerMin = ener;
            os.chdir("..")

if __name__ == "__main__":
    # Do nothing
    if verbose:
        print("psi4files.py")

