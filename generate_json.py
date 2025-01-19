#!/usr/bin/env python3

import argparse, glob, os, json, shutil, sys, gzip, io, time
from run_calcs import get_monoq
from multiprocessing import Pool, cpu_count

energies   = "energies"
monoq_data = "data/monomer_charges.csv"
Bohr       = 0.529177
au2Debye   = 2.5417464519
au2DebAng  = 1.3450342976

def get_sapt_map(method:str)->list:
    mydict = {}
    if method.find("sapt0") >= 0:
        mymap = [ { "key": "Total sSAPT0", "name": "InteractionEnergy", "n": 2 },
                  { "key": "Electrostatics sSAPT0", "name": "Electrostatics", "n": 1},
                  { "key": "Dispersion sSAPT0", "name": "Dispersion", "n": 1},
                  { "key": "Exchange sSAPT0", "name": "Exchange", "n": 1},
                  { "key": "Induction sSAPT0", "name": "Induction",  "n": 1} ]
        mydict = mymap
    else:
        comps = [ "Electrostatics   ", "Elst10,r    ", "Elst12,r    ",
                  "Exchange   ", "Exch10         ", "Exch10(S^2)    ",
                  "Exch11(S^2)    ", "Exch12(S^2)    ",
                  "Exch-Ind20,r  ", "Exch-Ind22  ",
                  "Induction   ", "Ind20,r     ", "Ind22        ",
                  "delta HF,r (2)  ",
                  "Exch-Disp20   ",
                  "Dispersion   ",
                  "Disp20  ", "Disp21  ", "Disp22 (SDQ)  ", "Disp22 (T)   ",
                  "Est. Disp22 (T) " ]
        if method == "sapt2+-aug-cc-pvdz":
            mydict = [ { "key": "Total SAPT2+", "name": "InteractionEnergy", "n": 2 } ]
        elif method == "sapt2+(ccd)dmp2-aug-cc-pvtz":
            mydict = [ { "key": "Total SAPT2+(CCD)dMP2", "name": "InteractionEnergy", "n": 2 } ]
            comps += [  "delta MP2,r (2)    ", "Disp2 (CCD)   ", "Disp22 (S) (CCD)  ", "Disp22 (T) (CCD)  ", "Est. Disp22 (T) (CCD)  " ]
        if len(mydict) > 0:
            for comp in comps:
                # Find location of energy on line
                cstrip = comp.strip()
                n = len(cstrip.split())
                mydict.append({ "key": comp, "name": cstrip, "n": n})
    return mydict

def to_element(elem:str)->str:
    # Convert element name to uppercase for first character only.
    myelem = elem.upper()[:1]
    if len(elem) > 1:
        myelem += elem[1:].lower()
    return myelem

def lineToAtom(line:str, nwords:int)->dict:
    myatom = {}
    words = line.split()
    if len(words) == nwords:
        try:
            myatom = {"elem": to_element(words[0]),
                      "coords": [ float(words[1]),
                                  float(words[2]),
                                  float(words[3]) ] }
        except ValueError:
            sys.exit("Incomprehensible line '%s' in logfile in %s" %
                     ( line, os.getcwd() ))
    return myatom

class Psi4Reader:
    def __init__(self, method:str, index:str, system:str,
                 verbose:bool, root:str, debug:bool=False):
        self.mydict = { "method": method,
                        "mols": [],
                        "energies": {},
                        "system": system,
                        "index": index }
        self.root    = root
        self.verbose = verbose
        self.debug   = debug
        self.msgs    = []
        self.sapt_map = {}
        if method.find("sapt") >= 0:
            if self.debug:
                self.msgs.append("SAPT detected in %s" % ( os.getcwd()))
            self.sapt_map = get_sapt_map(method)
            if len(self.sapt_map) == 0:
                self.msgs.append("Unknown SAPT method %s" % method)
        self.read_log_lines(index)
        self.read_out_lines(index)

    def read_log_lines(self, index:str):
        self.logLines = []
        self.logFile  = index + ".log"
        if os.path.exists(self.logFile):
            with open(self.logFile, "r") as log:
                self.logLines = log.readlines()
        else:
            self.logFile += ".gz"
            if os.path.exists(self.logFile):
                with gzip.open(self.logFile, 'rb') as ip:
                    with io.TextIOWrapper(ip, encoding='utf-8') as decoder:
                        self.logLines = decoder.readlines()
        if len(self.logLines) > 0:
            self.mydict["logFile"] = self.logFile

    def read_out_lines(self, index:str):
        self.outLines = []
        self.outFile  = index+".out"
        if not os.path.exists(self.outFile):
            self.outFile = self.mydict["system"] + ".out"
        if os.path.exists(self.outFile):
            with open(self.outFile, "r") as inf:
                self.outLines = inf.readlines()
        if len(self.outLines) > 0:
            self.mydict["outFile"] = self.outFile

    def read_electrostatics(self)->int:
        # Interpret the lines in order
        outmap = {
            "DIPOLE POLARIZABILITY": { "short": "alpha",      "nwords": 4, "subtype": 2, "value": 3, "factor": Bohr**3 },
            "DIPOLE":                { "short": "dipole",     "nwords": 3, "subtype": 1, "value": 2, "factor": au2Debye  },
            "QUADRUPOLE":            { "short": "quadrupole", "nwords": 3, "subtype": 1, "value": 2, "factor": au2DebAng  }
        }
        nlines    = 0
        foundELEC = False
        for line in self.outLines:
            words = line.strip().split()
            try:
                for ok in outmap.keys():
                    if line.find(ok) >= 0 and len(words) == outmap[ok]["nwords"]:
                        short = outmap[ok]["short"]
                        if not short in self.mydict:
                            self.mydict[short] = {}
                        self.mydict[short][words[outmap[ok]["subtype"]]] = float(words[outmap[ok]["value"] ]) * outmap[ok]["factor"]
                        foundELEC = True
                        nlines += 1
            except ValueError:
                self.msgs.append("Cannot read line '%s' in %s/%s\n" %
                                 ( line.strip(), os.getcwd(), self.outFile ) )
                break
        if foundELEC:
            self.mydict["jobtype"] = "Topology"
        return nlines
        
    def add_mol(self, compound:str, monoq:dict, atoms:list):
        self.mydict["mols"].append({ "compound": compound,
                                     "charge": monoq[compound]["charge"],
                                     "mult": monoq[compound]["mult"],
                                     "atoms": atoms })

    def out_error(self, nline: int, label:str):
        self.msgs.append("Error %s reading line '%s' in %s/%s" % ( label, self.outLines[nline].strip(), os.getcwd(), self.outFile ) )

    def read_out(self, compounds:list, monoq:dict)->bool:
        if len(self.outLines) == 0:
            self.msgs.append("Empty or missing outfile %s/%s" % ( os.getcwd(), self.outFile ))
            return False
        nline    = 0
        natomtot = 0
        for nm in range(len(compounds)):
            if not compounds[nm] in monoq:
                self.msgs.append("Compound %s missing from %s\n" % ( compounds[nm], monoq_data ))
                return
            try:
                natom = int(self.outLines[nline])
            except ValueError:
                self.out_error(nline, "reading natom")
                return False
            if len(self.outLines) < natom+2:
                self.msgs.append("File %s/%s incomplete" % ( os.getcwd(), self.outFile ) )
                return
            nline += 1
            atoms = []
            for j in range(natom):
                try:
                    words = self.outLines[nline+j].strip().split()
                    if len(words) == 4:
                        atoms.append( { "elem": to_element(words[0]),
                                        "coords": [  float(words[1]),
                                                     float(words[2]),
                                                     float(words[3]) ] } )
                except ValueError:
                    self.out_error(nline, "reading atom")
                    return False
            self.add_mol(compounds[nm], monoq, atoms)
            nline += natom
            natomtot += natom
        # Check for energy
        if nline < len(self.outLines):
            words = self.outLines[nline].strip().split()
            if len(words) >= 2 and words[0] == "Energy":
                try:
                    self.mydict["energies"]["energy"] = float(words[1])
                except ValueError:
                    self.out_error(nline, "reading energy")
                    return False
            nline += 1
        else:
            self.out_error(nline, "missing energy")
            return False

        # Routine will return the number of lines read for electrostatics
        oldnline = nline
        nline += self.read_electrostatics()
        self.msgs.append("In %s len(outLines) = %d, oldnline = %d nline = %d natomtot = %d" % ( os.getcwd(), len(self.outLines), oldnline, nline, natomtot ))
        if len(self.outLines)-nline == natomtot:
            # We likely have forces here
            for nm in range(len(compounds)):
                for j in range(len(self.mydict["mols"][nm]["atoms"])):
                    words = self.outLines[nline].strip().split()
                    try:
                        self.msgs.append("Adding a force")
                        force = [ float(words[0]), float(words[1]), float(words[2]) ]
                        self.mydict["mols"][nm]["atoms"][j]["forces"] = force
                    except ValueError:
                        self.out_error(nline, "reading forces")
                        return False
                    nline += 1
        return True

    def logcoordsToMydict(self, logcoords:list, monoq:dict, compounds:list)->str:
        n0     = 0
        elem   = "elem"
        coords = "coords"
        forces = "forces"
        self.mydict["mols"] = []
        
        msgs   = []
        for n in range(len(compounds)):
            atoms = []
            comp  = compounds[n]
            if not comp in monoq:
                self.msgs.append("Compound not present in %s" % monoq_data)
                return
            if monoq[comp]["natom"] == 0:
                self.msgs("Incorrect number of atoms for %s in %s" % (comp, monoq_data) )
                return
            for k in range(monoq[comp]["natom"]):
                atom           = {}
                if len(logcoords) <= n0+k:
                    msgs += ("Not enough atoms (%d) in logcoords, expected at least %d in %s" % ( len(logcoords), n0+k, os.getcwd() ))
                    return
                if elem in logcoords[n0+k] and coords in logcoords[n0+k]:
                    atom["elem"]   = logcoords[n0+k]["elem"]
                    atom["coords"] = logcoords[n0+k]["coords"]
                    if forces in logcoords[n0+k]:
                        atom[forces] = logcoords[n0+k][forces]
                    atoms.append(atom)
                else:
                    msgs += ("Incorrect logcoords '%s' in %s" % ( str(logcoords[n0+k]), os.getcwd()))
                    return
            self.add_mol(comp, monoq, atoms)
            n0 += monoq[comp]["natom"]
        return msgs

    def read_opt_log(self, compounds:list, monoq:dict)->bool:
        if len(self.logLines) == 0:
            self.msgs.append("Empty or missing logfile %s/%s\n" % ( os.getcwd(), self.logFile))
            return False
        close = 0
        logcoords = []
        for i in range(len(self.logLines)):
            line = self.logLines[i].strip()
            ener = None
            if line.find("Step    Total Energy") >= 0 and i+4 < len(self.logLines):
                words = self.logLines[i+4].strip().split()
                if len(words) >= 10:
                    try:
                        self.mydict["energies"]["energy"] = float(words[1])
                    except ValueError:
                        self.msgs.append("Cannot read energy on line %s in %s" % ( i+4, self.logFile))
            if line.find("Final optimized geometry and variables") >= 0:
                close = 1
            if close == 1 and line.find("Geometry") >= 0:
                close = 2
            if close == 2:
                atom = lineToAtom(line, 4)
                if len(atom.keys()) == 2:
                    logcoords.append(atom)

        self.logcoordsToMydict(logcoords, monoq, compounds)
        self.mydict["jobtype"] = "Opt"
        return True

    def read_sp_log(self, compounds:list, monoq:dict)->bool:
        if len(self.logLines) == 0:
            self.msgs.append("Empty logfile %s/%s\n" % ( os.getcwd(), self.logFile))
            return False
        centerDone = False
        energy     = None
        close = 0
        logcoords = []
        for ii in range(len(self.logLines)):
            line = self.logLines[ii].strip()
            if not centerDone and line.find("Center") >= 0 and line.find("Mass") >= 0:
                k = ii + 2
                atom = lineToAtom(self.logLines[k].strip(), 5)
                while (len(atom.keys()) == 2):
                    logcoords.append(atom)
                    k += 1
                    atom = lineToAtom(self.logLines[k].strip(), 5)
                
                centerDone = True
            if line.find("-Total Gradient:") >= 0:
                forces = []
                for k in range(ii+3, ii+3 + len(logcoords)):
                    if len(self.logLines) > k:
                        words = self.logLines[k].strip().split()
                        if len(words) == 4:
                            try:
                                forces.append( [ -float(words[1]), -float(words[2]), -float(words[3]) ] )
                            except ValueError:
                                self.msgs.append("Incorrect forces on line '%s'" % self.logLines[k].strip())
                if len(forces) == len(logcoords):
                    for k in range(len(logcoords)):
                        logcoords[k]["forces"] = forces[k]
            if line.find("Total Energy") >= 0:
                words = line.split()
                if len(words) >= 4:
                    try:
                        self.mydict["energies"]["energy"] = float(words[3])

                    except ValueError:
                        self.msgs.append("Cannot read energy from line '%s' in %s" % ( line.strip(), os.getcwd() ) )

        self.logcoordsToMydict(logcoords, monoq, compounds)
        self.mydict["jobtype"] = "SP"
        return True
         
    def read_sapt_log(self, compounds:list, monoq:dict)->bool:
        if len(self.logLines) == 0:
            self.msgs.append("Empty logfile %s/%s\n" % ( os.getcwd(), self.logFile))
            return False

        self.mydict[energies] = {}
        centerDone  = False
        logcoords   = []
        for ii in range(len(self.logLines)):
            line = self.logLines[ii].strip()
            if not centerDone and line.find("Center") >= 0 and line.find("Mass") >= 0:
                k = ii + 2
                atom = lineToAtom(self.logLines[k].strip(), 5)
                while (len(atom.keys()) == 2):
                    logcoords.append(atom)
                    k += 1
                    atom = lineToAtom(self.logLines[k].strip(), 5)

                centerDone = True
                # Skip to next line
                continue

            # Now check for SAPT energies, reset energy from outfile if present
            nmax = 0
            mkey = None
            for mm in self.sapt_map:
                words = line.split()
                if len(words) > 6:
                    if line.find(mm["key"]) >= 0:
                        n = mm["n"]
                        if n > nmax:
                            nmax = n
                            mkey = mm["name"]
            # Check whether we found something
            if mkey and nmax > 0:
                try:
                    if mkey in self.mydict[energies]:
                        self.msgs.append("Warning: duplicate key '%s' in mydict\n" % mkey)
                    self.mydict[energies][mkey] = float(words[nmax])/1000
                except ValueError:
                    self.msgs.append("Incomprehensible line '%s' n = %d\n" % (line.strip(), nmax))
                    return False

        self.logcoordsToMydict(logcoords, monoq, compounds)
        if self.debug:
            self.msgs.append("Found %d energies in %s/%s\n" %
                             ( len(self.mydict[energies]),
                               os.getcwd(), self.logFile ) )
        self.mydict["jobtype"] = "SP"
        return True

    def save_json(self, outFile:str):
        with open(outFile, "w") as jf:
            json.dump(self.mydict, jf, sort_keys=True, indent=4)
        if self.verbose:
            self.msgs.append("Stored %s/%s with %d energies\n" %
                             ( os.getcwd(), outFile, len(self.mydict[energies]) ) )

def check_json(rj:str)->bool:
    json_ok = False
    try:
        with open(rj, "r") as inf:
            mydict = json.load(inf)
            if ("energies" in mydict and len(mydict["energies"]) > 0 and
                "mols" in mydict and len(mydict["mols"]) > 0):
                json_ok = True
    except:
        json_ok = False
    return json_ok

def do_one(workdir:str, method:str, calc:str, compounds, monoq:dict, index:str,
           system:str, verbose:bool, debug:bool, force:bool, root:str):
    os.chdir(workdir)
    msgs = []
    if debug:
        msgs.append("Hello from %s" % os.getcwd())
    reader = Psi4Reader(method, index, system, verbose, root, debug)
    resjson = "results.json"
    json_ok = os.path.exists(resjson) and check_json(resjson)
    if not json_ok or force:
        # Make a backup of the existing file.
        if os.path.exists(resjson):
            os.rename(resjson, resjson+".bak")
        readOK = False
        # Read stuff from outfile first ...
        # readOK = reader.read_out(compounds, monoq)
        # ... but if we have a correct log file, rewrite the coordinates
        # and the energies
        if len(reader.logLines) > 0:
            if method.find("sapt") >= 0:
                readOK = reader.read_sapt_log(compounds, monoq)
            elif calc == 'opt':
                readOK = reader.read_opt_log(compounds, monoq)
            elif calc == 'sp':
                readOK = reader.read_sp_log(compounds, monoq)
        else:
            msgs += ( "No log file in %s" % os.getcwd())
        if readOK:
            reader.save_json(resjson)
        msgs += reader.msgs

    return msgs

    
def parse_args():
    desc = """Convert (gzipped) log file(s) and/or out files to results.json.
    If a results.json file exists already, the directory will be skipped, unless
    the --force flag is used. 
    """
    parser  = argparse.ArgumentParser(description=desc)
    defbasis = "aug-cc-pvdz"
    parser.add_argument("-basis", "--basis", help="Basis set, default is "+defbasis, type=str, default=defbasis)
    defmethod = "sapt2+"
    parser.add_argument("-method","--method", help="QM method default "+defmethod, type=str, default=defmethod)
    parser.add_argument("-f", "--force", help="If a results.json file exists no action will be taken, unless this flag is set", action="store_true")
    parser.add_argument("-v", "--verbose", help="Write verbose output", action="store_true")
    parser.add_argument("-dbg", "--debug", help="Write debugging output", action="store_true")
    maxcore = 16
    parser.add_argument("-ncore", "--ncore", help="Max number of cores to use (should be less than the number of available cores), default "+str(maxcore), type=int, default=maxcore)
    args = parser.parse_args()

    return args

if __name__ == "__main__":                
    args = parse_args()
    method = args.method + "-" + args.basis
    if not os.path.isdir(method):
        sys.exit("No such dir '%s'" % method)

    # Monomer properties
    monoq = get_monoq()
    logfnm = "generate_json.log"
    outf = open(logfnm, "w")
    os.chdir(method)
    root = os.getcwd()

    # Output
    all_msgs = []

    # Parallel processing
    ncpu     = 1
    num_jobs = 0
    pool     = None
    if args.ncore > 1:
        ncpu         = min(args.ncore, cpu_count())-1
        pool         = Pool(ncpu)
        pool_results = []
        print("Will run the script on %d cores" % ( ncpu ))
        args.verbose = False
        args.debug   = False
    complexes = { "1": "monomer", "2": "dimer" }
    for ckey in complexes.keys():
        for calc in [ "opt", "sp", "scans", "esp" ]:
            subdir = root + "/" + complexes[ckey] + "-" + calc
            if not os.path.isdir(subdir):
                continue
            os.chdir(subdir)
            for system in sorted(glob.glob("*")):
                if not os.path.isdir(system):
                    continue
                if args.verbose:
                    print("%s/%s" % ( subdir, system ) )
                compounds = system.split("#")
                for cc in compounds:
                    if not cc in monoq:
                        outf.write("Compound '%s' not found in %s, ignoring '%s'\n" % ( cc, monoq_data, system ))
                        continue
                sysdir = subdir + "/" + system
                os.chdir(sysdir)
                for index in glob.glob("*"):
                    if not os.path.isdir(index):
                        continue
                    inddir = sysdir + "/" + index
                    if args.debug:
                        outf.write("Now in dir %s\n" % os.getcwd())
                    if pool:
                        pool_results.append(pool.apply_async(do_one,
                                                             [ inddir, method, calc, compounds, monoq, index, system, args.verbose, args.debug, args.force, root ] ) )
                    else:
                        all_msgs += do_one( inddir, method, calc, compounds, monoq, index, system, args.verbose, args.debug, args.force )
                    num_jobs += 1

                os.chdir("..")
            os.chdir("..")
    os.chdir(root)

    if pool:
        # Check output from parallel runs
        print("\nWaiting for the %d workers and %d jobs to finish ..." % ( ncpu, num_jobs ) )
        iter  = 0
        sleep = 2
        all_finished = False
        while iter < 10*num_jobs and not all_finished:
            all_finished = True
            for res in pool_results:
                if not res.ready():
                    all_finished = False
                    
            time.sleep(sleep)
            iter += sleep
        for res in pool_results:
            all_msgs += res.get()
        
    for msg in all_msgs:
        outf.write("%s\n" % msg)
    outf.close()
    print("Please check your output and %s for more information" % logfnm)
