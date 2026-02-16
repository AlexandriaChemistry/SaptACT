#!/usr/bin/env python3

import glob, json, os, sys, math
import molprops
import get_mol_dict       as gmd
from atomic_heat_of_formation import *
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
    ShortDistance = 10
    InconsistentNumberOfAtoms = 11
    Unknown = 12

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
        Psi4Error.ShortDistance: "Psi4Error: Distance too short",
        Psi4Error.MonomerCharge: "Psi4Error: Compound missing in monomer_charge.csv",
        Psi4Error.InconsistentNumberOfAtoms: "Psi4Error: Inconsistent number of atoms",
        Psi4Error.Unknown: "Psi4Error: unknown"
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
debug                = False
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

    def add_props(self, myexp, mydict, temperature:float, logf):
        # Add electrostatic stuff
        mytype  = "Electronic"
        if "alpha" in mydict:
            alp     = mydict["alpha"]
            if None == alp:
                logf.write("Empty alpha in %s/%s\n" % ( os.getcwd(), self.dataFile ))
            else:
                if verbose:
                    logf.write("In %s\n" % os.getcwd())
                    logf.write("Found alpha {}\n".format(alp))
                average = str((alp["XX"]+alp["YY"]+alp["ZZ"])/3.0)
                myexp.add_polarisability(mytype, "Angstrom3", temperature, average, 0,
                                         str(alp["XX"]), str(alp["YY"]), str(alp["ZZ"]),
                                         str(alp["XY"]), str(alp["XY"]), str(alp["YZ"]))
        if "Dipole" in mydict:
            dip = mydict["Dipole"]
            if None == dip:
                logf.write("Empty dipole in %s/%s\n" % ( os.getcwd(), self.dataFile ))
            else:
                average = str(math.sqrt(dip["X"]**2 + dip["Y"]**2 + dip["Z"]**2))
                myexp.add_dipole(mytype, dip["unit"], temperature, average, 0,
                                 str(dip["X"]), str(dip["Y"]), str(dip["Z"]))
        if "Quadrupole" in mydict:
            quad = mydict["Quadrupole"]
            if None == quad:
                logf.write("Empty quadrupole in %s/%s\n" % ( os.getcwd(), self.dataFile ))
            else:
                quad_trace = (quad["XX"]+quad["YY"]+quad["ZZ"])/3.0
                myexp.add_quadrupole(mytype, "B", temperature,
                                     str(quad["XX"]-quad_trace), str(quad["YY"]-quad_trace), str(quad["ZZ"]-quad_trace), 0, 0, 0)
        freq = "frequencies"
        inten = "intensities"
        if freq in mydict and inten in mydict:
            for f in mydict[freq]:
                myexp.frequencies.append(str(f))
            for i in mydict[inten]:
                myexp.intensities.append(str(i))
        # Special stuff
        mytemp  = 298.15
        myphase = "gas"
        if "entropy" in mydict:
            myexp.add_cv_entropy(mydict["entropy"], mytemp, myphase)
        if "Cv" in mydict:
            myexp.add_cv_entropy(mydict["Cv"], mytemp, myphase)
        
    def read(self, args,
             molnames:list, monoq:dict, calc:str, mycomplex:str, filename:str,
             Mols:Molecules, ahof,
             temperature:float, potential:list, logf)->Psi4Error:
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
        # Check whether there are forces in json
        useForces = True
        for k in range(len(mydict[mols])):
            nForces = 0
            if atoms in mydict[mols][k]:
                for j in range(len(mydict[mols][k][atoms])):
                    if "forces" in mydict[mols][k][atoms][j]:
                        nForces += 1
                useForces = useForces and (nForces == len(mydict[mols][k][atoms]))
            else:
                useForces = False
        if verbose:
            logf.write("Found forces in %s: %s\n" % ( os.getcwd(), useForces ) )
        # Create new experiment
        myexp = molprops.Experiment("Theory", args.reference, "Psi4", args.method, args.basis,
                                    "conformation", self.jobtype, filename, useForces)
        # Add molecular properties
        self.add_props(myexp, mydict, temperature, logf)
        # Atoms etc.
        frag_charges = []
        sum_q_abs    = 0
        positive     = 1
        # TODO Fix charges!
        for j in range(self.nmols):
            if molnames[j] in monoq:
                mq = monoq[molnames[j]]["charge"]
                frag_charges.append(mq)
                positive = positive*mq
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
            if not atoms in mydict[mols][k]:
                logf.write("Incomplete json in %s/%s\n" % ( os.getcwd(), self.dataFile) )
                continue
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
                # These are atomic charges for use in heat of formation calcs
                self.charges.append(0)
            if (self.status == Psi4Error.OK and
                not self.MD[k].from_coords_elements(molnames[k], elements[offset:], coords, frag_charges[k])):
                logf.write("Cannot analyze coordinates using RDKit\n")
                self.status = Psi4Error.Coordinates
            if debug:
                logf.write("1. There are %d atoms for %s k = %d\n" % ( len(self.MD[k].atoms), molnames[k], k ) )
            if self.status == Psi4Error.OK and not self.firstDone:
                symmetry = 1
                mult     = 1
                mymol    = Mols.find_mol(molnames[k])
                if mymol:
                    symmetry = mymol.symmetry_number
                    mult     = mymol.mult
                    if mymol.stdinchi != self.MD[k].inchi:
                        logf.write("Warning: overriding generated InChi (%s) for %s with database InChi (%s)\n" % ( self.MD[k].inchi, molnames[k], mymol.stdinchi ) )
                        self.MD[k].inchi = mymol.stdinchi
                frag = molprops.Fragment(self.MD[k].inchi, frag_charges[k], mult, symmetry,
                                         myatoms, self.MD[k].mol_weight, self.MD[k].formula)
                if len(self.mp1.fragments) > 1:
                    logf.write("There are already %d fragments, skipping this one\n" % len(self.mp1.fragments))
                else:
                    self.mp1.add_fragment(frag)
                    for b in self.MD[k].bonds:
                        self.mp1.add_bond(offset+b[0], offset+b[1], self.MD[k].bonds[b])
            if debug:
                logf.write("2. There are %d atoms for %s k = %d\n" % ( len(self.MD[k].atoms), molnames[k], k ) )
            # Now add the atoms to the experiment
            for atom in range(len(self.MD[k].atoms)):
                obtype = self.MD[k].atoms[atom]["obtype"]
                # Atom numbers should match bonds
                if offset + atom >= len(elements):
                    logf.write(f"Inconsistent number of atoms for {molnames[k]}. There are {len(elements)} elements, offset = {offset}, natom = {natom}. len(self.MD[{k}].atoms) = {len(self.MD[k].atoms)}\n" )
                    self.status = Psi4Error.InconsistentNumberOfAtoms
                else:
                    # Note that atom numbers from zero, but MolDict from 1
                    myexp.add_atom(elements[offset+atom], obtype, offset+atom,
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
                deltaEmax = args.deltaEmax*(1+sum_q_abs)
                exch      = "Exchange"
                intener   = "InteractionEnergy"
                if self.edimerMin == edimerMinDefault:
                    logf.write("edimerMin has not been set\n")
                    self.status = Psi4Error.IncompleteJson
                if self.status == Psi4Error.OK:
                    if positive > 0:
                        # Compounds with equal signed charge
                        if exch in mydict[energies]:
                            if mydict[energies][exch] > deltaEmax:
                                logf.write(f"Skipping high {exch} energy {mydict[energies][exch]} for {molnames[0]}-{molnames[1]}\n")
                                self.status = Psi4Error.HighEnergy
                        else:
                            logf.write(f"No exchange energy available for {molnames[0]}-{molnames[1]}\n")
                            self.status = Psi4Error.IncompleteJson

                    else:
                        if exch in mydict[energies]:
                            if mydict[energies][exch] > deltaEmax:
                                logf.write(f"Skipping high {exch} energy {mydict[energies][exch]} for {molnames[0]}-{molnames[1]}\n")
                                self.status = Psi4Error.HighEnergy
                        else:
                            logf.write(f"No exchange energy available for {molnames[0]}-{molnames[1]}\n")
                            self.status = Psi4Error.IncompleteJson
#                        if (mydict[energies][intener] - self.edimerMin) > deltaEmax:
#                            logf.write("Skipping high %s energy %g for %s-%s\n" %
#                                       ( intener, mydict[energies][intener], molnames[0], molnames[1] ))
#                            self.status = Psi4Error.HighEnergy

                if self.status == Psi4Error.OK:
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
                
                    for energy in mydict[energies].keys():
                        myexp.add_energy(energy, "Hartree", 0, "gas", mydict[energies][energy])
                else:
                    logf.write("Skipping energies in %s because of %s" %
                               ( os.getcwd(), psi4msg(self.status) ) )
            else:
                # Monomers
                unit = "Hartree"
                if "unit" in mydict[energies]:
                    unit = mydict[energies]["unit"]
                deltaE0 = compute_dhform(mydict[energies]["energy"],
                                         elements, ahof,
                                         self.lots, self.charges, temperature)
                myexp.add_energy("DeltaE0", unit, 0, "gas", deltaE0)
                # Zero point vibrational energy
                ZPVE = "ZPVE"
                if ZPVE in mydict[energies]:
                    myexp.add_energy("ZPE", unit, 0, "gas", mydict[energies][ZPVE])

        if self.status == Psi4Error.OK:
            if self.nmols == 2:
                mydist = check_dist(allcoords)
                if mydist > args.rmax:
                    self.status = Psi4Error.LargeDistance
                    logf.write("Distance (%g) too large (max %g) between compounds in %s\n" %
                               ( mydist, args.rmax, os.getcwd() ) )
                elif mydist < args.rmin:
                    self.status = Psi4Error.ShortDistance
                    logf.write("Distance (%g) too short (max %g) between compounds in %s\n" %
                               ( mydist, args.rmin, os.getcwd() ) )
        if self.status ==  Psi4Error.OK and len(potential) > 0:
            self.add_potential(potential, myexp)
        if self.status == Psi4Error.OK:
            self.mp1.add_experiment(myexp)
        #if self.status == Psi4Error.OK:
        self.firstDone = True
        return self.status

    def add_potential(self, potential:list, myexp:molprops.Experiment):
        for np in range(len(potential)):
            p = potential[np]
            myexp.add_potential(str(np), "Hartree/e", "pm",
                                100*p["grid"][0], 100*p["grid"][1], 100*p["grid"][2],
                                p["value"])

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

