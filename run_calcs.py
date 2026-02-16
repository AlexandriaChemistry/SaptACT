#!/usr/bin/env python3

import argparse, copy, os, glob, math, shutil, sys
import random as rnd
from get_csv_rows import *
from diatomics import get_diatomics
import elements as elements

BOHR    = 0.529177
psi4ACT = os.getcwd()

special_basis = { "aug-cc-pvtz": { "K": "def2-TZVPP", "Ca": "def2-TZVPP", "Rb": "def2-TZVPP", "Cs": "def2-TZVPP",
                                   "I": "aug-cc-pvtz-pp", "Kr": "aug-cc-pwcvtz-pp", "Xe": "aug-cc-pwcvtz-pp" },
#                                   "I": "aug-cc-pwcvtz-pp", "Kr": "aug-cc-pwcvtz-pp", "Xe": "aug-cc-pwcvtz-pp" },
                  "aug-cc-pvqz": { "K": "def2-QZVPP", "Ca": "def2-QZVPP", "Rb": "def2-QZVPP", "Cs": "def2-QZVPP",
                                   "I": "aug-cc-pwcvqz-pp", "Kr": "aug-cc-pwcvqz-pp", "Xe": "aug-cc-pwcvqz-pp" },
                  "aug-cc-pv5z": { "I": "aug-cc-pwcv5z-pp", "Kr": "aug-cc-pwcv5z-pp", "Xe": "aug-cc-pwcv5z-pp" },
                  "aug-cc-pvtz-pp": { "I": "aug-cc-pvtz-pp", "Kr": "aug-cc-pvtz-pp", "Xe": "aug-cc-pvtz-pp" },
                  "6-311G": { "I": "koko" }
              }

def make_form(elem:str, charge:int)->str:
    form = elem
    if 0 > charge:
        form += str(charge) + "-"
    elif 0 < charge:
        form += str(charge) + "+"
    return form

def get_monoq()->dict:
    monoq  = {}
    for line in get_csv_rows(psi4ACT+"/data/monomer_charge.csv", 4, delim=","):
        monoq[line[0]] = { "charge": int(line[1]), "mult": int(line[2]), "natom": int(line[3]) }
    return monoq

def get_vdwradii(vdwr = "data/vdwradii.dat")->dict:
    vdwdict = {}
    if os.path.exists(vdwr):
        with open(vdwr, "r") as inf:
            for line in inf:
                words = line.strip().split()
                if len(words) == 3:
                    try:
                        # Input is in nanometer, but DO NOT move to Angstrom here
                        # since this is used in GROMACS tools as well.
                        vdwdict[words[1]] = float(words[2])
                    except ValueError:
                        print("Invalid line '%s' in %s" % ( line.strip(), vdwr ))
    return vdwdict

def get_dimer_selection(selection:str)->list:
    monoq  = get_monoq()
    print("Read monomer data for %d compounds" % len(monoq.keys()))
    dimers = []
    for line in get_csv_rows(selection, 2, delim="#"):
        mon1 = line[0]
        mon2 = line[1]
        # Check for ACT formatted files
        if mon2.find("|") > 0:
            mon2 = mon2.split("|")[0]
        if mon1 in monoq and mon2 in monoq:
            dimers.append({ "mon1": mon1, "q1": monoq[mon1]["charge"], "m1": monoq[mon1]["mult"], "nat1": monoq[mon1]["natom"],
                            "mon2": mon2, "q2": monoq[mon2]["charge"], "m2": monoq[mon1]["mult"], "nat2": monoq[mon2]["natom"],
                            "pair": ("%s#%s" % (mon1, mon2))})
        else:
            print("Either '%s' or '%s' missing from the data/monomer_charge.csv file" % ( mon1, mon2 ))
    return dimers

def get_monomer_selection(selection:str)->list:
    monoq    = get_monoq()
    monomers = []
    for line in get_csv_rows(selection, 1):
        mon1 = line[0]
        q1   = 0
        m1   = 1
        if mon1 in monoq:
            monomers.append({ "mon1": mon1, "q1": monoq[mon1]["charge"], "m1": monoq[mon1]["mult"], "nat1": monoq[mon1]["natom"]})
        else:
            print("%s missing from the data/monomer_charge.csv file" % ( mon1 ))
    return monomers

def get_xyz(basename:str, resetCOM:bool, atomprops:dict):
    if basename.endswith("xyz") and os.path.exists(basename):
        xyz = basename
    else:
        xyz = (f"{psi4ACT}/xyz/monomers/{basename}.xyz")
    if not os.path.exists(xyz):
        print("Cannot find %s, I am in %s" % (xyz, os.getcwd()))
        return None, None
    with open(xyz, "r") as inf:
        lines  = inf.readlines()
        if len(lines) == 0:
            print("Empty file %s" % xyz)
            return None, None
        natom  = int(lines[0])
        elem   = []
        coords = []
        com    = [ 0.0, 0.0, 0.0 ]
        tmass  = 0
        for i in range(natom):
            www = lines[i+2].strip().split()
            elem.append(www[0])
            thiscoord = [ float(www[1]), float(www[2]), float(www[3]) ]
            coords.append(thiscoord)
            mymass = 0
            if www[0] in atomprops:
                mymass = atomprops[www[0]]["mass"]
                tmass += mymass
            for m in range(3):
                com[m] += mymass*thiscoord[m]
        if resetCOM and tmass > 0:
            for m in range(3):
                com[m] /= tmass
            for xx in coords:
                for m in range(3):
                    xx[m] -= com[m]
        return elem, coords
    sys.exit("Cannot interpret %s" % xyz)
    return None, None

def update_vdwradii(factor:float):
    # get vdw radii
    vdw = "vdwradii.dat"
    vdw_dict = get_vdwradii(vdw)
    # replace all occurrences of the required string
    with open(vdw, "w") as outf:
        for i in vdw_dict:
            outf.write("??? %s %g\n" % (i, vdw_dict[i]*factor) )

def orient(labels:list, coords:list, alpha:float, beta:float, gamma:float):
    # The rotation matrix, Maple output
    #
    # Rot :=
    #
    #    [cos(beta) cos(gamma) , -cos(beta) sin(gamma) , sin(beta)]
    #
    #    [sin(alpha) sin(beta) cos(gamma) + cos(alpha) sin(gamma) ,
    #
    #    -sin(alpha) sin(beta) sin(gamma) + cos(alpha) cos(gamma) ,
    #
    #    -sin(alpha) cos(beta)]
    #
    #    [-cos(alpha) sin(beta) cos(gamma) + sin(alpha) sin(gamma) ,
    #
    #    cos(alpha) sin(beta) sin(gamma) + sin(alpha) cos(gamma) ,
    #
    #    cos(alpha) cos(beta)]

    cosa = math.cos(alpha)
    sina = math.sin(alpha)
    cosb = math.cos(beta)
    sinb = math.sin(beta)
    cosg = math.cos(gamma)
    sing = math.sin(gamma)
    A = [ [ cosb * cosg, -cosb * sing, sinb ],
          [ sina * sinb * cosg + cosa * sing, -sina * sinb * sing + cosa * cosg, -sina * cosb ],
          [-cosa * sinb * cosg + sina * sing, cosa * sinb * sing + sina * cosg, cosa * cosb ] ]

    myx = []
    for m in range(len(coords)):
        newx = [ 0, 0, 0 ]
        for i in range(3):
            for j in range(3):
                newx[i] += A[j][i] * coords[m][j]
        myx.append(newx)

    return myx
    
def compute_mindist(x1:list, x2:list):
    mindist2 = 1000
    mdvec    = [ 0, 0, 0 ]
    ij       = ( 0, 0 )
    for i in range(len(x1)):
        for j in range(len(x2)):
            md2 = 0
            newvec = [ 0, 0, 0 ]
            for m in range(3):
                newvec[m] = x2[j][m]-x1[i][m]
                md2 += (newvec[m])**2
            if md2 < mindist2:
                mindist2 = md2
                mdvec = newvec
                ij = (i, j)
    if 1000 == mindist2:
        sys.exit("Internal error computing mindist")
    return math.sqrt(mindist2), mdvec, ij

def compute_coms(x1:list, x2:list):
    com1 = [ 0, 0, 0 ]
    for i in range(len(x1)):
        for m in range(3):
            com1[m] += x1[i][m]
    com2 = [ 0, 0, 0 ]
    for i in range(len(x2)):
        for m in range(3):
            com2[m] += x2[i][m]
    comvec = [ 0, 0, 0 ]
    for m in range(3):
        com1[m] /= len(x1)
        com2[m] /= len(x2)
    return com1, com2

def subtract_com(xyz:list, com:list)->list:
    xyzcom = []
    for i in range(len(xyz)):
        xyzc = []
        for m in range(3):
            xyzc.append(xyz[i][m]-com[m])
        xyzcom.append(xyzc)
    return xyzcom

def get_monomers(root:str, userfile:str, mon1:str, nat1:int, mon2:str, nat2:int):
    pdbfile = None
    swapped = False
    if len(userfile) > 0:
        pdbfile = userfile
    else:
        pdbfile = ( "%s/xyz/dimers/%s#%s.pdb" % ( root, mon1, mon2 ))
        if not os.path.exists(pdbfile):
            print("Cannot find %s, giving up" % pdbfile)
            return None, None, None, None
            #pdbfile = ( "%s/xyz/dimers/%s#%s.pdb" % ( root, mon2, mon1 ))
            #swapped = True
    if not os.path.exists(pdbfile):
        print("No such file %s" % pdbfile )
        return None, None, None, None
    if (pdbfile.find("%s#%s" % ( mon1, mon2 )) < 0 and
        pdbfile.find("%s#%s" % ( mon2, mon1 )) < 0):
        print("Mismatch between filename %s and monomers %s resp. %s" % ( pdbfile, mon1, mon2 ))
        return None, None, None, None
    label = []
    xyz   = []
    with open(pdbfile, "r") as inf:
        atom = 0
        for line in inf:
            if line.find("ATOM") >= 0 or line.find("HETATM") >= 0:
                words = line.split()
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    coord  = [ x, y, z ]
                    elem   = words[2]
                    if atom == 0 or atom == nat1:
                        label.append([])
                        xyz.append([])
                    label[-1].append(elem)
                    xyz[-1].append(coord)
                    atom += 1
                except ValueError:
                   print(" incomprehensible line '%s' in %s" % ( line.strip(), pdbfile))
                   return None, None, None, None
    if len(label) == 2 and len(xyz) == 2:
        if swapped:
            return label[1],  xyz[1], label[0], xyz[0]
        else:
            return label[0],  xyz[0], label[1], xyz[1]
    else:
        print("len(label) = %d len(xyz) = %d trying to get_monomers from %s" % ( len(label), len(xyz), pdbfile ) )
        return None, None, None, None

def extract_grid_points_pdb(surface_files:list):
    x_values = []
    y_values = []
    z_values = []
    for surface_file_path in surface_files:
        print(surface_file_path)
        if os.path.isfile(surface_file_path):
            inFileStream = open(surface_file_path, 'r')
            for line in inFileStream:
                if (len(line.split()) == 10):
                    Atom, number1, dot, dot, number2, x, y, z, number3, number4 = line.split()
                    if Atom == "ATOM":
                        try:
                            x=float(x)
                            y=float(y)
                            z=float(z)
                            x_values.append(x)
                            y_values.append(y)
                            z_values.append(z)
                        except:
                            continue
            inFileStream.close()
        else:
            print("{0} is not available.".format(surface_file))	        
    return x_values, y_values, z_values
    
def write_grid_points(x, y, z, labels:list, coords: list):
  x_grid = []
  y_grid = []
  z_grid = []	

  output_file = "grid.dat"
  output = open(output_file, 'w')
  
  # atom coordinates
  for value in range(len(labels)):
    output.write("{0} {1} {2}\n".format(coords[value][0], coords[value][1], coords[value][2])) 
    x_grid.append(coords[value][0])
    y_grid.append(coords[value][1])
    z_grid.append(coords[value][2])

  # grid points
  for value in range(len(x)):
    output.write("{0} {1} {2}\n".format(x[value],y[value],z[value]))
    x_grid.append(x[value])
    y_grid.append(y[value])
    z_grid.append(z[value])
  output.close()  
  
  # this writes grid.xyz that can be opened in e.g. pymol
  if False:
      output_file = "grid.xyz"
      output = open(output_file, 'w')
      output.write("{0}\n".format(len(x_grid)))
      output.write("grid\n")
      for value in range(len(x_values)):
          output.write("H   {0} {1} {2}\n".format(x_values[value],y_values[value],z_values[value]))

      for value in range(len(x)):
          output.write("H   {0} {1} {2}\n".format(x[value],y[value],z[value]))
      output.close() 
  return x_grid, y_grid, z_grid   

def next_idir(idir:int):
    my_idir = ("%04d" % idir)
    while os.path.isdir(my_idir):
        idir += 1
        my_idir = ("%04d" % idir)
    return idir, my_idir
    

class Psi4jobs:
    '''Class to manage jobs run using the Psi4 software on HPC cluster'''
    def __init__(self, args):
        self.method      = None
        self.basis       = None
        self.memory      = args.memory
        self.ncores      = args.ncores
        self.hours       = args.hours
        self.frozen_core = args.frozen_core
        self.dryrun      = args.dryrun
        self.submit      = args.submit
        self.scratch     = args.scratch_dir
        self.optimize    = args.optimize
        self.esp         = args.esp
        self.frequency   = args.frequency
        self.nconf       = args.nconf
        self.gromacs     = args.gromacs
        self.root        = os.getcwd()
        self.verbose     = args.verbose
        self.vdw0        = args.vdw0
        self.nlayer      = args.nlayer
        self.dlayer      = args.layer

    def set_method(self, method):
        self.method = method
        
    def set_basis(self, basis):
        self.basis = basis
        
    def lot(self)->str:
        lot = self.method+"-"+self.basis
        if self.frozen_core:
            lot += "-fc"
        return lot

    def run_one_job(self, job:str):
        if self.dryrun:
            if self.verbose:
                print("Not submitting %s" % job)
            return
        if shutil.which(self.submit):
            os.system(self.submit + " " + job)
        else:
            os.chmod(job, 0o755)
            # TODO: Check for path
            os.system("./%s 2>&1 output.dat &" % job)
    
    def write_job_header(self, outf, outfile:str, mult:int):
        outf.write("#!/usr/bin/env python3\n")
        outf.write("#SBATCH -t %d:00:00\n" % self.hours)
        outf.write("#SBATCH -c %d\n" % self.ncores)
        sn = "SNIC_RESOURCE"
        if sn in os.environ:
            outf.write("#SBATCH -A naiss2024-3-13\n")
            if os.environ[sn] == "dardel":
                outf.write("#SBATCH -p shared\n")
        else:
            outf.write("#SBATCH -p CLUSTER,CLUSTER-AMD\n")
        outf.write("import os, sys\n")
        outf.write("import numpy as np\n")
        outf.write("import psi4 as psi4\n")
        if self.method == "MP2":
            self.write_calculate_polarizabilities(outf)
        outf.write("psi4.core.set_num_threads(%d)\n" % self.ncores)
        p4iopt = { "cachelevel": 1, "print": 1, 'damping_percentage': 20 }
        p4sopt = { "guess": "read" }
        reference = 'rhf'
        if mult > 1:
            reference = 'rohf'
        p4sopt['reference'] = reference
            
        if self.frozen_core:
            outf.write(", 'freeze_core': 'true'")
        if self.optimize:
            p4sopt["opt_coordinates"] = "cartesian"
            p4sopt["fail_on_maxiter"] = "false"
#            p4sopt["scf_initial_accelerator"] =  "none"
#            p4sopt["G_CONVERGENCE"] = "GAU_VERYTIGHT"
            p4iopt["geom_maxiter"] = 10000
            p4iopt["soscf_max_iter"] = 35

        outf.write("psi4.set_options({")
        comma = False
        for pp in p4sopt:
            if comma:
                outf.write(", ")
            outf.write("'%s': '%s'" % ( pp, p4sopt[pp] ))
            comma = True
        for pp in p4iopt:
            if comma:
                outf.write(", ")
            outf.write("'%s': %d" % ( pp, p4iopt[pp] ))
        outf.write("})\n")

        outf.write("psi4.set_memory(%.0f)\n" % (self.memory*self.ncores*1.0e6))
        outf.write("psi4_io = psi4.core.IOManager.shared_object()\n")
        if None != self.scratch:
            outf.write("tmpdir = \"%s\"\n" % self.scratch)
        else:
            outf.write("tmpdir = \"PDC_TMP\"\n")
#            outf.write("tmpdir = \"TMPDIR\"\n")
            outf.write("if tmpdir in os.environ:\n")
            outf.write("    myhead,mytail = os.path.split(os.getcwd())\n")
            outf.write("    myhead2,mytail2 = os.path.split(myhead)\n")
            outf.write("    tmpdir = os.environ[tmpdir]+\"/\" + mytail2 + \"/\" + mytail\n")
            outf.write("    os.makedirs(tmpdir, exist_ok=True)\n")
            outf.write("else:\n")
            outf.write("    tmpdir = \".\"\n")
        outf.write("psi4_io.set_default_path(tmpdir)\n")
        outf.write("psi4.core.set_output_file('%s', False)\n" % outfile)

    def gen_diatomic_scan_script(self, name:str, ai:str, aj:str, dist:float, charge:int, mult:int)->str:
        jobname = name + "-scan.py"
        rmin    = round (0.84*dist, 2)
        rmax    = 1.2*dist
        npoints = int((rmax - rmin)*200)+1
        with open(jobname, "w") as outf:
            outfile = name + ".out"
            self.write_job_header(outf, outfile, mult)
            outf.write("geometry= \"\"\"\n")
            outf.write(" %d %d\n" % ( charge, mult ) )
            outf.write(" %s 0 0 0\n" % ai)
            outf.write(" %s 0 0 {0}\n" % aj)
            outf.write("\"\"\"\n")
            outf.write("table = []\n")
            outf.write("npoints    = %d\n" % npoints) 
            outf.write("for i in range(npoints):\n")
            outf.write("    r         = %g+0.005*i\n" % ( rmin ))
            outf.write("    geom      = psi4.geometry(geometry.format(r))\n")
            outf.write("    psi4.basis_helper(\"\"\"\n")
            outf.write("assign %s\n" % basis)
            numeric = False
            if self.basis in special_basis:
                if ai in special_basis[self.basis]:
                    outf.write("assign %s %s\n" % ( ai, special_basis[self.basis][ai] ) )
                    numeric = True
                    shutil.copy(("%s/basis/%s.gbs" % ( self.root, special_basis[self.basis][ai])), ".")
                if aj in special_basis[self.basis] and aj != ai:
                    outf.write("assign %s %s\n" % ( aj, special_basis[self.basis][aj] ) )
                    numeric = True
                    shutil.copy(("%s/basis/%s.gbs" % ( self.root, special_basis[self.basis][aj])), ".")
            if self.method.find("ccsd") >= 0:
                numeric = True
            outf.write("\"\"\", name=\'dz_PLUS\')\n")
            forces = False
            outf.write("    try:\n")
            if forces:
                outf.write("        grad, wfn = psi4.gradient(\"%s\", molecule=geom, return_wfn=True" % self.method)
                if numeric:
                    outf.write(", dertype=0")
                outf.write(")\n")
                outf.write("        myener    = wfn.energy()\n")
                outf.write("        forces    = grad.to_array()\n")
                outf.write("        fz        = -forces[0][2]\n")
            else:
                outf.write("        myener = psi4.energy(\"%s\", molecule=geom)\n" % self.method)
                outf.write("        fz     = 0\n")
            outf.write("        table.append([r, myener, fz])\n")
            outf.write("    except Exception as ex:\n")
            outf.write("        print(ex)\n")
            outf.write("        continue\n")
            outf.write("""# Now find minimum through bisection
npoints = len(table)
if npoints < 3:
    sys.exit(\"Not enough data points!\")
ia = 0
ib = npoints-1
# Find minimum
for i in range(npoints):
    if table[i][1] < table[ia][1]:
        ia = i
if ia == npoints-1:
    sys.exit(\"Incomplete scan\")
if table[ia-1][1] < table[ia+1][1]:
    ia = ia -1
ib = ia + 1
toler = 0.00001 # Angstrom
while (table[ib][0]-table[ia][0]) > toler:
    rx     = (table[ia][0]+table[ib][0])/2
    geom   = psi4.geometry(geometry.format(rx))
    myener = psi4.energy(\"%s\", molecule=geom)
    fz     = 0
    table.append([rx, myener, fz])
    if table[ia][1] < table[ib][1]:
        ib = len(table)-1
    else:
        ia = len(table)-1
""" % self.method)
            outf.write("if len(table) > 0:\n")
            outf.write("    with open(\"energy.xvg\", \"w\") as outf:\n")
            outf.write("        for tab in sorted(table):\n")
            outf.write("            outf.write(\"%20.15f  %20.15f\\n\" % ( tab[0], tab[1] ))\n")
            outf.write("    with open(\"force.xvg\", \"w\") as outf:\n")
            outf.write("        for tab in sorted(table):\n")
            outf.write("            outf.write(\"%20.15f  %20.15f\\n\" % ( tab[0], tab[2] ))\n")
        return jobname

    def gen_atomization_script(self, ai:str, charge:int, mult:int, args)->str:
        jobname = ai + ".py"
        with open(jobname, "w") as outf:
            outfile = ai + ".out"
            self.write_job_header(outf, outfile, mult)
            outf.write("geometry= \"\"\"\n")
            outf.write("%d %d\n" % (charge, mult ))
            outf.write(" %s 0 0 0\n" % ai)
            basis = self.basis
            os.popen("ls")
            if self.basis in special_basis:
                if ai in special_basis[self.basis]:
                    basis = special_basis[self.basis][ai]
                    shutil.copy(("%s/basis/%s.gbs" % ( self.root, basis)), ".")
            outf.write("\"\"\"\n")
            outf.write("geom    = psi4.geometry(geometry)\n")
            outf.write("myener  = psi4.energy(\"%s/%s\", molecule=geom)\n" % ( self.method, basis ))
            if args.frozen_core:
                outf.write("mybasis = \"%s-fc\"\n" % basis)
            else:
                outf.write("mybasis = \"%s\"\n" % basis)
            outf.write("with open('atomization_energy.dat', 'w') as aefile:\n")
            outf.write("    aefile.write('Atomization energy %f %s\\n' % ( myener, mybasis ))\n")

        return jobname

    def run_diatomic_scan(self):
        diatomics = get_diatomics()
        lot = self.lot()
        os.makedirs(lot, exist_ok=True)
        os.chdir(lot)
        scandir = "diatomic-scans"
        os.makedirs(scandir, exist_ok=True)
        os.chdir(scandir)
        for dim in diatomics:
            os.makedirs(dim, exist_ok=True)
            os.chdir(dim)
            doit = True
            for fn in [ "energy.xvg" ]:
                if os.path.exists(fn) and os.path.getsize(fn) > 0:
                    doit = False
            if doit:
                job = self.gen_diatomic_scan_script(dim,
                                                    diatomics[dim]["ai"],
                                                    diatomics[dim]["aj"],
                                                    diatomics[dim]["distance"],
                                                    diatomics[dim]["charge"],
                                                    diatomics[dim]["mult"])
                print("Generated new script in %s" % job)
                self.run_one_job(job)
            else:
                print("There is output already in %s" % os.getcwd())
            os.chdir("..")
        os.chdir("../..")

    def run_atoms(self):
        lot = self.lot()
        os.makedirs(lot, exist_ok=True)
        os.chdir(lot)
        monomers = "atomization"
        os.makedirs(monomers, exist_ok=True)
        os.chdir(monomers)
        for elem in elements.atomprops.keys():
            os.makedirs(elem, exist_ok=True)
            os.chdir(elem)
            if not os.path.exists("atomization_energy.dat"):
                job = self.gen_atomization_script(elements.atomprops[elem]["symbol"], 
                                                  elements.atomprops[elem]["charge"],
                                                  elements.atomprops[elem]["mult"], args)
                self.run_one_job(job)
            os.chdir("..")
        os.chdir("../..")

    def write_calculate_polarizabilities(self,outf):
        outf.write("def calculate_polarizabilities():\n")
        outf.write("    pert = 1.8897261250e-5\n")
        outf.write("    lambdas = [pert, -pert, 2.0*pert, -2.0*pert]\n")
        outf.write("    method = \"MP2\"\n")
        outf.write("    perturbed_energies = np.zeros((len(lambdas),3))\n")
        outf.write("    perturbed_dipoles = np.zeros((len(lambdas),3,3))\n")
        outf.write("    psi4.set_options({\"perturb_h\":True,\n")
        outf.write("                    \"perturb_with\":\"dipole\",})\n")
        outf.write("    for step, l in enumerate(lambdas):\n")
        outf.write("        for d in range(3):\n")
        outf.write("            perturb_dipole = list(l * np.eye(3)[d])\n")
        outf.write("            psi4.set_options({\"perturb_dipole\": perturb_dipole,})\n")
        outf.write("            perturbed_energies[step,d] = psi4.properties(method, properties=['dipole'])\n")
        outf.write("            dip = psi4.core.variable(method + \" dipole\")\n")
        outf.write("            for m in range(3):\n")
        outf.write("                perturbed_dipoles[step,m,d] = dip[m]\n")
        outf.write("    alpha_5pt = np.zeros((3,3))\n")
        outf.write("    for i in range(3):\n")
        outf.write("        alpha_5pt[i] = -(8.0*perturbed_dipoles[0,:,i] - 8.0*perturbed_dipoles[1,:,i] - perturbed_dipoles[2,:,i] + perturbed_dipoles[3,:,i]) / (12.0*pert)\n")
        outf.write("    return alpha_5pt\n")
    
    def write_esp_input(self, outf):
        if not method.strip() in [ "MP2", "ccsd(t)" ]:
            outf.write("_, opt_wfn = psi4.optimize(\"%s\", molecule=geom, return_wfn=True)\n" % method)
            outf.write("mydict[\"energies\"][\"energy\"], wfn = psi4.energy(\"%s\", molecule=opt_wfn.molecule(), properties_origin=[\"NUCLEAR_CHARGE\"], properties=[\"GRID_ESP\", \"DIPOLE_POLARIZABILITIES\"], return_wfn=True)\n" % method)
            outf.write("mydict[\"polar\"] = psi4.properties(\"%s\", properties=[\"DIPOLE_POLARIZABILITIES\"])\n" % method)
        else:
            outf.write("_, opt_wfn = psi4.optimize(\"%s\", molecule=geom, return_wfn=True)\n" % method)            
            outf.write("mydict[\"energies\"][\"energy\"], wfn = psi4.energy(\"%s\", molecule=opt_wfn.molecule(), properties_origin=[\"NUCLEAR_CHARGE\"], properties=[\"GRID_ESP\"], return_wfn=True)\n" % method)            
        outf.write("oeprops = psi4.core.OEProp(wfn)\n")
        outf.write("oeprops.add(\"GRID_ESP\")\n")
        outf.write("oeprops.add(\"MULTIPOLE(4)\")\n")
        outf.write("oeprops.compute()\n")
        outf.write("psi4.core.print_variables()\n")
        
    def write_esp_output(self, outf):
        outf.write("    compvars = [ \"XX\", \"XY\", \"XZ\", \"YY\", \"YZ\", \"ZZ\" ]\n")
        outf.write("    wfnvars = wfn.variables(\"DIPOLE\")\n")
        outf.write("    qvars = [ \"X\", \"Y\", \"Z\" ]\n")
        outf.write("    dip = wfnvars[\"DIPOLE\"]\n")
        outf.write("    mydict[\"dipole\"] = {}\n")
        outf.write("    for m in range(3):\n")
        outf.write("        dm = dip[m]\n")
        outf.write("        result.write(\"DIPOLE %s %10g\\n\" % ( qvars[m], dm ))\n")
        outf.write("        mydict[\"dipole\"][qvars[m]] = dm\n")
        outf.write("    quad = wfnvars[\"QUADRUPOLE\"]\n")
        outf.write("    mydict[\"quadrupole\"] = {}\n")
        outf.write("    qindex = 0\n")
        outf.write("    for m in range(3):\n")
        outf.write("        for n in range(m,3):\n")
        outf.write("            qmn = quad[m][n]\n")
        outf.write("            result.write(\"QUADRUPOLE %s %g\\n\" % ( compvars[qindex], qmn ))\n")
        outf.write("            mydict[\"quadrupole\"][compvars[qindex]] = qmn\n")
        outf.write("            qindex += 1\n")
        if self.method.strip() != "MP2":
            outf.write("    for comp in compvars:\n")
            outf.write("        pvar    = \"DIPOLE POLARIZABILITY \" + comp\n")
            outf.write("        result.write(\"%s %s\\n\" % ( pvar, str(psi4.variable(pvar))))\n")
        else:
            outf.write("    new_geometry = \"\\n symmetry c1\" + geometry\n")
            outf.write("    new_geom = psi4.geometry(new_geometry)\n")
            outf.write("    polarizabilities = calculate_polarizabilities()\n")
            outf.write("    qindex = 0\n")
            outf.write("    mydict[\"alpha\"] = {}\n")
            outf.write("    for i in range(3):\n")
            outf.write("        for j in range(i,3):\n")
            outf.write("            alphaij = polarizabilities[i][j]\n")
            outf.write("            result.write(\"DIPOLE POLARIZABILITY %s %g\\n\" % ( compvars[qindex], alphaij ))\n")
            outf.write("            mydict[\"alpha\"][compvars[qindex]] = alphaij\n")
            outf.write("            qindex += 1\n")

    def write_dimer_input(self, myname:str, remark:str,
                          q1:int, m1:int, label1:list, xyz1ori:list,
                          q2:int, m2:int, label2:list, xyz2ori:list)->str:
        if self.esp:
            labels = []
            for lll in [ label1, label2 ]:
                for l in lll:
                    labels.append(l)
            xyzs   = []
            for xyz in [ xyz1ori, xyz2ori ]:
                for xxx in xyz:
                    xyzs.append(xxx)
            self.create_grid(labels, xyzs)
        myjob = myname + ".py"
        with open(myjob, "w") as outf:
            outfile = myname + ".log"
            self.write_job_header(outf, outfile, max(m1, m2))
            extrabasis = {}
            outf.write("geometry= \"\"\"\n")
            outf.write(" %d %d\n" % (q1, m1))
            for i in range(len(label1)):
                outf.write(" %s %s %s %s\n" % ( label1[i], xyz1ori[i][0],
                                                xyz1ori[i][1], xyz1ori[i][2] ))
                if self.basis in special_basis and label1[i] in special_basis[self.basis]:
                    extrabasis[label1[i]] = 1
            outf.write(" --\n")
            outf.write(" %d %d\n" % (q2, m2))
            for i in range(len(label2)):
                outf.write(" %s %s %s %s\n" % ( label2[i], xyz2ori[i][0],
                                                xyz2ori[i][1], xyz2ori[i][2] ))
                if self.basis in special_basis and label2[i] in special_basis[self.basis]:
                    extrabasis[label2[i]] = 1
            if self.esp:
                # See https://github.com/psi4/psi4/issues/2964
                outf.write(" noreorient\n")
                outf.write(" nocom\n")

            outf.write("\"\"\"\n")
            outf.write("geom = psi4.geometry(geometry)\n")
            outf.write("psi4.basis_helper(\"\"\"\n")
            outf.write("assign %s\n" % self.basis)
            for eb in extrabasis.keys():
                outf.write("assign %s %s\n" % ( eb, special_basis[self.basis][eb] ))
                shutil.copy(("%s/basis/%s.gbs" % (self.root, special_basis[basis][eb])), ".")
            outf.write("\"\"\")\n")
            useSAPT = method.find("sapt") >= 0
            outf.write("forces = None\n")
            outf.write("mydict = {}\n")
            outf.write("mydict[\"energies\"] = {}\n")
            if len(remark) > 0:
                outf.write("mydict[\"remark\"] = \"%s\"\n" % remark)
            saptENER = {}
            if useSAPT:
                outf.write("mydict[\"energies\"][\"energy\"] = psi4.energy(\"%s\")\n" % self.method)
                saptENER = { "Electrostatics": 'SAPT ELST ENERGY',
                             "Exchange": 'SAPT EXCH ENERGY',
                             "Induction": 'SAPT IND ENERGY',
                             "Dispersion": 'SAPT DISP ENERGY',
                             "InteractionEnergy": 'SAPT TOTAL ENERGY' }
                for enm,esapt in saptENER.items():
                    outf.write("mydict[\"energies\"][\"%s\"] = psi4.variable('%s')\n" % ( enm, esapt))

            elif self.esp:
                self.write_esp_input(outf)
            elif self.optimize:
                outf.write("try:\n")
                outf.write("    mydict[\"energies\"][\"energy\"], wfn = psi4.gradient(\"%s\", molecule=geom, return_wfn=True)\n" % self.method)
                outf.write("    geom = np.asarray(wfn.molecule().geometry())\n")
                output = myname + ".xyz"
                outf.write("    with open(\"%s\", \"w\") as result:\n" % output)
                labels = label1+label2
                outf.write("        result.write(\"%5d\\n\")\n" % (len(labels)))
                outf.write("        dist = %f*((geom[0][0]-geom[1][0])**2 + (geom[0][1]-geom[1][1])**2 + (geom[0][2]-geom[1][2])**2 )**0.5\n" % BOHR)
                outf.write("        result.write(\"%12.8f  %12.8f\\n\" % ( mydict[\"energies\"][\"energy\"], dist))\n")
                for i in range(len(labels)):
                    outf.write("        xyz = [ geom[%d][0]*%f, geom[%d][1]*%f, geom[%d][2]*%f ]\n" %
                               ( i, BOHR, i, BOHR, i, BOHR ))
                    outf.write("        result.write(\"%5s  %%12.8f  %%12.8f  %%12.8f\\n\" %% (xyz[0], xyz[1], xyz[2]))\n" % (labels[i]))
                outf.write("except psi4.ConvergenceError:\n")
                outf.write("    print(\"Did not converge\\n\")\n")
                if self.frequency:
                    outf.write("mydict[\"freq_en\"], mydict[\"freq_wfn\"] = psi4.frequency(\"%s/%s\", molecule=geom, return_wfn=True, dertype='gradient')\n" % (self.method, basis))

            else:
                outf.write("grad, wfn = psi4.gradient(\"%s\", molecule=geom, return_wfn=True" % self.method)
                if len(extrabasis.keys()) > 0 or self.method.find("ccsd") >= 0:
                    outf.write(", dertype=0")
                outf.write(")\n")
                outf.write("mydict[\"energies\"][\"energy\"] = wfn.energy()\n")
                outf.write("forces    = grad.to_array()\n")
            if not self.optimize:
                if self.frequency:
                    outf.write("psi4.driver.gradient(\"%s/%s\", return_wfn=True)\n" % (self.method, basis))
                    outf.write("freq_en, freq_wfn = psi4.frequency(\"%s/%s\", molecule=geom, return_wfn=True, dertype='gradient')\n" % (self.method, basis))
                output = myname + ".out"
                outf.write("with open(\"%s\", \"w\") as result:\n" % output)
                outf.write("    result.write(\"%5d\\n\")\n" % (len(label1)))

                for i in range(len(label1)):
                    outf.write("    result.write(\" %s %s %s %s\\n\")\n" % ( label1[i], xyz1ori[i][0],
                                                                             xyz1ori[i][1], xyz1ori[i][2] ))
                outf.write("    result.write(\"%5d\\n\")\n" % (len(label2)))
                for i in range(len(label2)):
                    outf.write("    result.write(\" %s %s %s %s\\n\")\n" % ( label2[i], xyz2ori[i][0],
                                                                             xyz2ori[i][1], xyz2ori[i][2] ))
                # Formatting weird, sorry!
                outf.write("    result.write(\"Energy %s %s\\n\" %s mydict[\"energies\"][\"energy\"])\n" % ( "%12.8f", remark, "%" ))
                if self.esp:
                    self.write_esp_output(outf)
                outf.write("    if None != forces and None != forces.all():\n")
                outf.write("        for f in forces:\n")
                outf.write("            result.write(\"%12.8f  %12.8f  %12.8f\\n\" % (f[0], f[1], f[2]))\n")

        return myjob

    def create_grid(self, labels:list, coords:list):
        pdb_file = "mymol.pdb"
        with open(pdb_file, "w") as outf:
            outf.write("MODEL  1\n")
            for i in range(len(labels)):
                outf.write("%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" %
                           ( "HETATM", i+1, labels[i], " ", "RES", 'A', 1, " ", coords[i][0], coords[i][1], coords[i][2],
                             0.0, 0.0, labels[i], "0") )
            outf.write("ENDMDL\n")
        
        vdW_value = []
        for i in range(self.nlayer):
            vdW_value.append(self.vdw0+self.dlayer*i)

        ndots = 200 
        surface_files = []
        folders = []
        for distance in vdW_value:
            folder = "vdW_{0}".format(distance)
            surfpdb = ("surface_%d.pdb" % (int(distance*10) ) )
            surface_files.append("%s/%s" % ( folder, surfpdb))
            os.makedirs(folder, exist_ok=True)
            os.chdir(folder)
            folders.append(folder)
            vdwr = "vdwradii.dat"
            os.system("cp %s/data/%s ." % ( psi4ACT, vdwr))
            os.system("cp ../%s ." % pdb_file)
            update_vdwradii(distance)
            os.system('export GMXDATA=.; echo 0 | %s sasa -f %s -s %s -q %s -probe 0.001 -ndots %d -noprot -o area_%g' %
                      ( self.gromacs, pdb_file, pdb_file, surfpdb, ndots, distance) )
            #if False:
            #    os.unlink(vdwr)
            os.chdir("..")

        x,y,z = extract_grid_points_pdb(surface_files)
        x_grid, y_grid, z_grid = write_grid_points(x, y, z, labels, coords)
        if False:
            for sf in folders:
                shutil.rmtree(sf)

    def write_monomer_input(self, myname:str, q1:int, m1:int, label1:list, xyz1ori:list)->str:
        if self.esp:
            self.create_grid(label1, xyz1ori)
        if not label1 or len(label1) == 0:
            print("Incorrect input to write_monomer_input")
            return None
        myjob = myname + ".py"
        with open(myjob, "w") as outf:
            outfile = myname + ".log"
            self.write_job_header(outf, outfile, m1)
            extrabasis = {}
            outf.write("geometry= \"\"\"\n")
            outf.write(" %d %d\n" % (q1, m1))
            for i in range(len(label1)):
                outf.write(" %s %s %s %s\n" % ( label1[i], xyz1ori[i][0],
                                                xyz1ori[i][1], xyz1ori[i][2] ))
                if basis in special_basis and label1[i] in special_basis[basis]:
                    extrabasis[label1[i]] = 1
                    shutil.copy(("%s/basis/%s.gbs" % (self.root, special_basis[basis][label1[i]])), ".")
            if self.esp:
                # See https://github.com/psi4/psi4/issues/2964
                outf.write(" noreorient\n")
                outf.write(" nocom\n")
                outf.write(" symmetry c1\n")
            outf.write("\"\"\"\n")
            outf.write("geom = psi4.geometry(geometry)\n")
            outf.write("mydict = {}\n")
            outf.write("mydict[\"energies\"] = {}\n")
            outf.write("psi4.basis_helper(\"\"\"\n")
            outf.write("assign %s\n" % basis)
            numeric = False
            for eb in extrabasis.keys():
                outf.write("assign %s %s\n" % ( eb, special_basis[basis][eb] ))
                numeric = True
            outf.write("\"\"\")\n")
            outf.write("forces = None\n")
            if self.esp:
                self.write_esp_input(outf)
            elif self.optimize:
                outf.write("mydict[\"energies\"][\"energy\"], wfn = psi4.optimize(\"%s\", molecule=geom, return_wfn=True)\n" % self.method)
                outf.write("geom = np.asarray(wfn.molecule().geometry())\n")
                if self.frequency:
                    outf.write("psi4.frequency(\"%s\")\n" % self.method)
                # After optimization write output xyz file for later processing
                output = myname + "-output.xyz"
                outf.write("with open(\"%s\", \"w\") as result:\n" % output)
                outf.write("    result.write(\"%5d\\n\")\n" % (len(label1)))
                if len(label1) == 2:
                    outf.write("    dist = %f*((geom[0][0]-geom[1][0])**2 + (geom[0][1]-geom[1][1])**2 + (geom[0][2]-geom[1][2])**2 )**0.5\n" % BOHR)
                    outf.write("    result.write(\"%12.8f  %12.8f\\n\" % ( mydict[\"energies\"][\"energy\"], dist))\n")
                else:
                    outf.write("    result.write(\"%12.8f\\n\" % ( mydict[\"energies\"][\"energy\"]))\n")
                for i in range(len(label1)):
                    outf.write("    xyz = []\n")

                    for m in range(3):
                        outf.write("    xyz.append(geom[%d][%d]*%f)\n" % ( i, m, BOHR ))
                    outf.write("    result.write(\"%5s  %%12.8f  %%12.8f  %%12.8f\\n\" %% (xyz[0], xyz[1], xyz[2]))\n" % (label1[i] ))
            else:
                outf.write("grad, wfn = psi4.gradient(\"%s\", molecule=geom, return_wfn=True" % method)
                if numeric:
                    outf.write(", dertype=0")
                outf.write(")\n")
                outf.write("mydict[\"energies\"][\"energy\"]    = wfn.energy()\n")
                outf.write("forces    = grad.to_array()\n")
            if not self.optimize:
                output = myname + ".out"
                outf.write("with open(\"%s\", \"w\") as result:\n" % output)
                outf.write("    result.write(\"%5d\\n\")\n" % (len(label1)))
                for i in range(len(label1)):
                    outf.write("    result.write(\" %s %s %s %s\\n\")\n" % ( label1[i], xyz1ori[i][0],
                                                                             xyz1ori[i][1], xyz1ori[i][2] ))
                outf.write("    result.write(\"Energy %12.8f\\n\" % ( mydict[\"energies\"][\"energy\"] ) )\n")
                outf.write("    polstr = \"\"\n")
                if self.esp:
                    self.write_esp_output(outf)

        return myjob
    
    def get_dist_fractions(self, mindist:float, maxdist:float, ndist:int)->list:
        df = []
        for idist in range(ndist):
            if ndist > 1:
                df.append((mindist + (idist*(maxdist-mindist))/(ndist-1)))
            else:
                df.append(0)
        return df

    def run_one_dimer(self, dimer, userfile:str, vdwdict:dict):
        ndist    = args.ndist
        if len(userfile) > 0:
            norient = 0
        else:
            norient  = args.norient
        ncores   = args.ncores
        memory   = args.memory
        if norient == 0:
            label1, xyz1, label2, xyz2 = get_monomers(self.root, userfile, dimer["mon1"], dimer["nat1"], dimer["mon2"], dimer["nat2"])
            if None == label1 or None == label2:
                print("Cannot find monomer information for %s or %s" % (dimer["mon1"], dimer["mon2"])) 
                return
        else:
            label1, xyz1 = get_xyz(dimer["mon1"], True, elements.atomprops)
            label2, xyz2 = get_xyz(dimer["mon2"], True, elements.atomprops)
        if not xyz1 or not label1 or 0 == len(xyz1) or len(label1) != len(xyz1):
            print("Could not find coordinates for %s" % ( dimer["mon1"] ))
            return
        if not xyz2 or not label2 or 0 == len(xyz2) or len(label2) != len(xyz2):
            print("Could not find coordinates for %s" % ( dimer["mon2"] ))
            return
        com1, com2 = compute_coms(xyz1, xyz2)
        Pi  = math.pi
        Pi2 = 2*Pi
        mindist = args.mindist # Angstrom or fraction
        maxdist = args.maxdist # id.
        idir    = 200
        mydir   = ( "%s#%s" % ( dimer["mon1"], dimer["mon2"] ))
        mydir2  = ( "%s#%s" % ( dimer["mon2"], dimer["mon1"] ))
        if os.path.exists(mydir2) and not os.path.exists(mydir):
            mydir = mydir2
        os.makedirs(mydir, exist_ok=True)
        os.chdir(mydir)
        if self.verbose:
            print("There are %d orientations and %d distances for %s" % ( norient, ndist, mydir ))
            
        delta = 0.02
        maxorient = norient
        df = self.get_dist_fractions(mindist, maxdist, ndist)
        if 0 == norient:
            xyz1ori0 = xyz1[:]
            xyz2ori0 = xyz2[:]
            ori_dist, dist_vec, ij = compute_mindist(xyz1ori0, xyz2ori0)
            print(dist_vec)
            if ori_dist == 0:
                sys.exit("Distance between monomers is zero, check your input in %s" % ( os.getcwd()) )
                return
            # Loop over distances
            xyz1ori = xyz1ori0[:]
            for dist_frac in df:
                xyz2ori = copy.deepcopy(xyz2ori0)
                for i in range(len(xyz2)):
                    for m in range(3):
                        xyz2ori[i][m] = xyz2ori0[i][m] + dist_frac*dist_vec[m]/ori_dist
                # print("dist_frac %f dist %f xyz2ori %f xyz2ori0 %f" % ( dist_frac, xyz2ori[0][2]-xyz1ori[0][2], xyz2ori[0][2], xyz2ori0[0][2] ))
                # Make directory
                idir, my_idir = next_idir(idir)
                os.makedirs(my_idir, exist_ok=False)
                os.chdir(my_idir)
                myjob = self.write_dimer_input(my_idir, userfile,
                                               dimer["q1"], dimer["m1"], label1, xyz1ori,
                                               dimer["q2"], dimer["m2"], label2, xyz2ori)
                self.run_one_job(myjob)
                                  
                os.chdir("..")
                idir += 1
            
        else:
            # Loop over orientations
            for iorient in range(norient):
                # Put both compounds on the origin
                xyz1com  = subtract_com(xyz1, com1)
                xyz2com  = subtract_com(xyz2, com2)
                xyz1ori0 = orient(label1, xyz1com, rnd.uniform(0, Pi2), rnd.uniform(0, Pi2), rnd.uniform(0, Pi))
                xyz2ori0 = orient(label2, xyz2com, rnd.uniform(0, Pi2), rnd.uniform(0, Pi2), rnd.uniform(0, Pi))
                dist_vec = [ 0, 0, 1 ]
                
                # Loop over distances
                for dist_frac in df:
                    # Make directory
                    idir, my_idir = next_idir(idir)
                    os.makedirs(my_idir, exist_ok=False)
                    os.chdir(my_idir)
                    xyz1ori = copy.deepcopy(xyz1ori0)
                    xyz2ori = copy.deepcopy(xyz2ori0)

                    if args.absdist:
                        for i in range(len(xyz2ori)):
                            xyz2ori[i][2] = xyz2ori0[i][2] + dist_frac
                        if self.verbose:
                            print(f"xyz1ori {xyz1ori} xyz2ori {xyz2ori} xyz2ori0 {xyz2ori0}")
                    else:
                        #xyz1ori = orient(label1, xyz1com, rnd.uniform(0, Pi2), rnd.uniform(0, Pi2), rnd.uniform(0, Pi))
                        #xyz2ori = orient(label2, xyz2com, rnd.uniform(0, Pi2), rnd.uniform(0, Pi2), rnd.uniform(0, Pi))
                        rel_dist = 0
                        while rel_dist < dist_frac:
                            # Shift all atoms in compound 2 by delta in the Z direction
                            for m in range(len(xyz2ori)):
                                xyz2ori[m][2] += delta
                            # Compute relative distance between the compounds
                            new_ori_dist, new_dist_vec, new_ij = compute_mindist(xyz1ori, xyz2ori)
                            # Need to convert to Angstrom here
                            vdw_sum  = 10*(vdwdict[label1[new_ij[0]]] + vdwdict[label2[new_ij[1]]])
                            rel_dist = math.sqrt(new_dist_vec[0]**2+new_dist_vec[1]**2+new_dist_vec[2]**2)/vdw_sum
                            if self.verbose:
                                print("new_ori_dist %g vdw_sum %g rel_dist %g dist_frac %g" % ( new_ori_dist, vdw_sum, rel_dist, dist_frac ))

                    if self.verbose:
                        print("dist_frac %f dist %f xyz2ori %f xyz2ori0 %f" %
                              ( dist_frac, xyz2ori[0][2]-xyz1ori[0][2], xyz2ori[0][2], xyz2ori0[0][2] ))
                    myjob = self.write_dimer_input(my_idir, userfile,
                                                   dimer["q1"], dimer["m1"], label1, xyz1ori,
                                                   dimer["q2"], dimer["m2"], label2, xyz2ori)
                    self.run_one_job(myjob)
                                  
                    os.chdir("..")
                    idir += 1
        os.chdir("..")
    
    def opt_dimer(self, dimer, userfile:str, force:bool):
        label1, xyz1, label2, xyz2 = get_monomers(self.root, userfile, dimer["mon1"], dimer["nat1"], dimer["mon2"], dimer["nat2"])
        if None == label1 or None == label2 or not xyz1 or not xyz2:
            print("Could not read dimer file for %s#%s" % ( dimer["mon1"], dimer["mon2"] ))
            return
        if self.verbose:
            print("Will try and optimize dimer " + dimer)
        mydir   = ( "%s#%s" % ( dimer["mon1"], dimer["mon2"] ))
        os.makedirs(mydir, exist_ok=True)
        os.chdir(mydir)
        output = mydir + ".xyz"
        if self.optimize and (force or not os.path.exists(output)):
            myjob = self.write_dimer_input(mydir, "optimization",
                                           dimer["q1"], dimer["m1"], label1, xyz1, 
                                           dimer["q2"], dimer["m2"], label2, xyz2)
            self.run_one_job(myjob)
        elif self.esp:
            idir    = 0
            my_idir = ("%04d" % idir)
            while os.path.isdir(my_idir):
                idir += 1
                my_idir = ("%04d" % idir)
            os.makedirs(my_idir, exist_ok=True)
            os.chdir(my_idir)
            myjob = self.write_dimer_input(mydir, userfile,
                                           dimer["q1"], dimer["m1"], label1, xyz1, 
                                           dimer["q2"], dimer["m2"], label2, xyz2)
            self.run_one_job(myjob)
            os.chdir("..")
            
        os.chdir("..")
    
    def run_one_monomer(self, monomer, monoq):
        if len(monomer.strip()) == 0:
            print("Empty monomer name '%s'" % monomer)
            return
        nconf = self.nconf
        xyzs  = []
        curdir = os.getcwd()
        if nconf == 1:
            xyzs.append(monomer)
        else:
            # Check whether we have single molecule files
            os.chdir(self.root + "/gromacs")
            if os.path.exists(monomer):
                os.chdir(monomer)
                for index in glob.glob("*"):
                    xyz = f"{index}/{index}.xyz"
                    if os.path.isdir(index) and os.path.exists(xyz):
                        xyzs.append(os.path.abspath(xyz))
                os.chdir("..")
            os.chdir(curdir)
            if nconf == 0:
                nconf = len(xyzs)
            else:
                nconf = min(nconf, len(xyzs))
        mydir   = monomer
        os.makedirs(mydir, exist_ok=True)
        os.chdir(mydir)
        idir = 0
        charge = 0
        mult   = 1
        if monomer in monoq:
            charge = monoq[monomer]["charge"]
            mult   = monoq[monomer]["mult"]
        for xyz in xyzs:
            label1, xyz1 = get_xyz(xyz, True, elements.atomprops)
            if label1 and len(label1) > 0:
                idir, my_idir = next_idir(idir)
                os.makedirs(my_idir, exist_ok=False)
                os.chdir(my_idir)
                myjob = self.write_monomer_input(my_idir, charge, mult, label1, xyz1)       
                self.run_one_job(myjob)
                os.chdir("..")
            else:
                print(f"Problem reading {xyz}. Got label {label1}")
        # Back to where we started before this job
        os.chdir(curdir)
    
    def run_dimers(self, dimerfile:str, userfiles:bool,
                   vdwdict:dict, force:bool):
        if not os.path.exists(dimerfile):
            sys.exit("Selection file %s does not exist" % dimerfile)
        dimers  = get_dimer_selection(dimerfile)
        print("There are %d dimers in the selection %s" % ( len(dimers), dimerfile ))
        if self.optimize:
            scandir = "dimer-opt"
        elif self.esp:
            scandir = "dimer-esp"
        else:
            scandir = "dimer-scans"
        lot     = self.lot()
        mydir   = ("%s/%s" % ( lot, scandir ))
        os.makedirs(mydir, exist_ok=True)
        os.chdir(mydir)
        if self.verbose:
            print("mydir = %s" % os.getcwd())
        
        for d in dimers:
            if userfiles:
                template = d["mon1"] + "#" + d["mon2"]
                for pdbfile in glob.glob(("%s/user/%s*.pdb" % ( psi4ACT, template ))):
                    if self.optimize or self.esp:
                        self.opt_dimer(d, pdbfile, force)
                    else:
                        self.run_one_dimer(d, pdbfile, vdwdict)
            else:
                if self.optimize or self.esp:
                    self.opt_dimer(d, "", force)
                else:
                    self.run_one_dimer(d, "", vdwdict)
        os.chdir("../../")
    
    def run_monomers(self, monomerfile:str):
        monomers = []
        with open(monomerfile, "r") as inf:
            for line in inf:
                monomers.append(line.strip())
        print("There are %d monomers in the selection %s" % ( len(monomers), monomerfile ))
        monoq  = get_monoq()
        scandir = "monomer"
        if self.optimize:
            scandir += "-opt"
        if self.esp:
            scandir += "-esp"
        if not self.optimize and not self.esp:
            scandir += "-sp"
        lot = self.lot()
        mydir   = ("%s/%s" % ( lot, scandir ))
        os.makedirs(mydir, exist_ok=True)
        os.chdir(mydir)
        for m in monomers:
            print("Will run %s" % m)
            self.run_one_monomer(m, monoq)
        os.chdir("../../")
    
def add_lot_args(parser):
    defbasis = [ "aug-cc-pvtz", "aug-cc-pvtz-pp", "aug-cc-pvqz", "def2-TZVPP", "def2-QZVPP" ]
    parser.add_argument("-basis", "--basis", nargs="+", help="Basis set(s)", type=str, default=defbasis)
    defmethod = [ "B3LYP", "WB97X", "HF", "B3LYP-D3BJ", "WB97X-D3BJ", "MP2", "ccsd(t)", "PWPB95-D3BJ" ]
    parser.add_argument("-method","--method", nargs="+", help="QM method(s)", type=str, default=defmethod)

def parse_args():
    desc = """
Run Psi4 calculations for use in ACT. When running dimers, the product of the ndist and norient and
the number of dimers in the selection file equals the total number of calculations. It can be a lot!
If the number of orientations is zero, the script will look for minimized structures for the dimers
and make a distance scan based on that. In this case the minimum and maximum distance are interpreted
as a fraction of the distance of the closest atoms in the minimized input structure, e.g. from 0.9 to 1.5
times that distance. If norient is larger than zero, monomers will be randomly rotated and
a distance scan done along the z-coordinate. The mindist and maxdist given by the user are interpreted
as being relative to the sum of the Van der Waals radii of the closest atoms.
    """
    parser  = argparse.ArgumentParser(description=desc)
    add_lot_args(parser)
    parser.add_argument("-fc","--frozen_core", help="Use frozen core orbitals, useful for expensive methods. Note that this is the default in other QM codes, such as Gaussian, but not in Psi4. If used, the basis set name will be extended with '-fc'", action="store_true")
    parser.add_argument("-v", "--verbose", help="Write debugging output", action="store_true")
    parser.add_argument("-diatomics", "--diatomics", help="Run distance scans for diatomic compounds", action="store_true")
    parser.add_argument("-atomization", "--atomization", help="Compute atomization energies", action="store_true")
    parser.add_argument("-dimers", "--dimers", help="Compute dimer interactions based on compounds in a selection file, please provide file name with this flag", type=str, default=None)
    parser.add_argument("-monomers", "--monomers", help="Compute energy for monomeric compounds in a selection file, please provide file name with this flag", type=str, default=None)
    parser.add_argument("-user", "--userfiles", help="Employ input coordinates supplied in the 'user' directory with name similar to those in xyz/dimers, namely mol1#mol2#number.pdb. You still need to provide the -dimer flag", action="store_true")
    gromacs = "gmx"
    parser.add_argument("-gmx", "--gromacs", help="Command to run GROMACS, only needed for ESP calculations. Default "+gromacs, type=str, default=gromacs)
    parser.add_argument("-esp", "--esp", help="Compute electrostatic potential for monomers", action="store_true")
    parser.add_argument("-opt", "--optimize", help="Optimize geometry, only for monomers", action="store_true")
    vdw0 = 1.4
    parser.add_argument("-vdw0", "--vdw0", help="Factor on Van der Waals radii for creating the first grid layer around atoms when doing ESP calculations. Default "+str(vdw0), type=float, default=vdw0)
    nlayer = 4
    parser.add_argument("-nlayer", "--nlayer", help="Number of grid layers around atoms when doing ESP calculations. Default "+str(nlayer), type=int, default=nlayer)
    dlayer = 0.2
    parser.add_argument("-layer", "--layer", help="Thickness of grid layers relative to the Van der Waals radius around atoms when doing ESP calculations. Default "+str(dlayer), type=float, default=dlayer)
    parser.add_argument("-freq", "--frequency", help="Compute frequencies. This will also optimize the geometry, only for monomers", action="store_true")
    nconf = 1
    parser.add_argument("-nconf", "--nconf", help="Number of conformations for each monomer in single point calculations, zero means all available conformation in gromacs/compound. default "+str(nconf), type=int, default=nconf)
    ndist   = 1
    parser.add_argument("-ndist", "--ndist", help="Number of dimer distances to compute, see helptext. Default "+str(ndist), type=int, default=ndist)
    mindist = 2.5
    parser.add_argument("-mindist", "--mindist", help="Minimum dimer distances to compute, default "+str(mindist), type=float, default=mindist)
    maxdist = 5.0
    parser.add_argument("-maxdist", "--maxdist", help="Maximum dimer distances to compute, default "+str(maxdist), type=float, default=maxdist)
    parser.add_argument("-absdist", "--absdist", help="Use the distances as absolute numbers", action="store_true")
    norient = 1
    parser.add_argument("-norient", "--norient", help="Number of dimer orientation to compute, see help text. Default "+str(norient), type=int, default=norient)
    ncores = 4
    parser.add_argument("-ncores", "--ncores", help="Number of cores to use for each calculation, default "+str(ncores), type=int, default=ncores)
    memory = 3000
    parser.add_argument("-mem", "--memory", help="Memory per core to use for each calculation, default "+str(memory)+" Mb", type=int, default=memory)
    hours = 24
    parser.add_argument("-hours", "--hours", help="Amount of hours to queue for, default "+str(hours), type=int, default=hours)
    parser.add_argument("-dryrun", "--dryrun", help="Write scripts but do not perform the runs", action="store_true")
    parser.add_argument("-force" , "--force", help="Run optimization even though it may have been done before", action="store_true")
    host = "HOST"
    scr  = None
    if host in os.environ and os.environ[host] == "csb.bmc.uu.se":
        scr = "/scratch"
    parser.add_argument("-scratch", "--scratch_dir", help="Directory to write scratch files. If not given, the TMPDIR environment variable will be checked first, and if it does not exist the current working directory will be used. ", type=str, default=scr)
    submit = "sbatch"
    parser.add_argument("-submit", "--submit", help="Command to submit to batch system, default "+submit, type=str, default=submit)
    args = parser.parse_args()
    if args.verbose:
        print("Turning on debugging.")
    return args

if __name__ == "__main__":
    args    = parse_args()
    if args.esp and None == shutil.which(args.gromacs):
        sys.exit("GROMACS command %s not in search path" % args.gromacs)
    if not args.dryrun and None == shutil.which("psi4"):
        sys.exit("Psi4 command not found in search path")
    psi4    = Psi4jobs(args)
    methods = args.method
    bases   = args.basis
    elements.get_atomprops()
    if type(methods) == str:
        methods = [ args.method ]
    if type(bases) == str:
        bases = [ args.basis ]

    vdwdict = get_vdwradii()
    for method in methods:
        psi4.set_method(method)
        for basis in bases:
            psi4.set_basis(basis)
            if args.diatomics:
                psi4.run_diatomic_scan()
            if args.atomization:
                psi4.run_atoms()
            if args.dimers:
                psi4.run_dimers(args.dimers, args.userfiles, vdwdict, args.force)
            if args.monomers:
                psi4.run_monomers(args.monomers)
