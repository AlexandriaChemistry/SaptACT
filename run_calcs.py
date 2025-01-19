#!/usr/bin/env python3

import argparse, copy, os, glob, math, shutil, sys
import random as rnd
from get_csv_rows import *
import elements as elements

BOHR    = 0.529177
psi4ACT = os.getcwd()

special_basis = { "aug-cc-pvtz": { "K": "def2-TZVPP", "Ca": "def2-TZVPP", "Rb": "def2-TZVPP", "Cs": "def2-TZVPP",
                                   "I": "aug-cc-pvtz-pp", "Kr": "aug-cc-pwcvtz-pp", "Xe": "aug-cc-pwcvtz-pp" },
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

def get_vdwradii()->dict:
    vdwr = (psi4ACT+"/data/vdwradii.dat")
    vdwdict = {}
    with open(vdwr, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if len(words) == 3:
                try:
                    # Input is in nanometer, move to Angstrom
                    vdwdict[words[1]] = 10.0*float(words[2])
                except ValueError:
                    print("Invalid line '%s' in %s" % ( line.strip(), vdwr ))
    return vdwdict

def get_dimer_selection(selection:str)->list:
    monoq  = get_monoq()
    dimers = []
    for line in get_csv_rows(selection, 2, delim="#"):
        mon1 = line[0]
        mon2 = line[1]
        q1 = 0
        q2 = 0
        m1 = 1
        m2 = 1
        if mon1 in monoq and mon2 in monoq:
            dimers.append({ "mon1": mon1, "q1": monoq[mon1]["charge"], "m1": monoq[mon1]["mult"], "nat1": monoq[mon1]["natom"],
                            "mon2": mon2, "q2": monoq[mon2]["charge"], "m2": monoq[mon1]["mult"], "nat2": monoq[mon2]["natom"],
                            "pair": ("%s#%s" % (mon1, mon2))})
        else:
            print("Either %s or %s missing from the data/monomer_charge.csv file" % ( mon1, mon2 ))
    return dimers

def get_xyz(basename:str, resetCOM:bool, atomprops:dict):
    if basename.endswith("xyz") and os.path.exists(basename):
        xyz = basename
    else:
        xyz = ("../../xyz/monomers/%s.xyz" % basename)
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
        self.root        = os.getcwd()
        self.verbose     = args.verbose

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
            outf.write("#SBATCH -A naiss2023-5-531\n")
            if os.environ[sn] == "dardel":
                outf.write("#SBATCH -p shared\n")
        else:
            outf.write("#SBATCH -p CLUSTER,CLUSTER-AMD\n")
        outf.write("import os, sys\n")
        outf.write("import numpy as np\n")
        outf.write("import psi4 as psi4\n")
        outf.write("psi4.core.set_num_threads(%d)\n" % self.ncores)
        p4iopt = { "cachelevel": 1, "print": 1, 'damping_percentage': 20 }
        p4sopt = { "guess": "read" }
        reference = 'rhf'
        if mult > 1:
            reference = 'rohf'
        p4sopt['reference'] = reference
            
        if self.frozen_core:
            outf.write(", 'freeze_core': 'true'")

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

    def write_dimer_input(self, myname:str, remark:str,
                          q1:int, m1:int, label1:list, xyz1ori:list,
                          q2:int, m2:int, label2:list, xyz2ori:list)->str:
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

            else:
                outf.write("grad, wfn = psi4.gradient(\"%s\", molecule=geom, return_wfn=True" % self.method)
                if len(extrabasis.keys()) > 0 or self.method.find("ccsd") >= 0:
                    outf.write(", dertype=0")
                outf.write(")\n")
                outf.write("mydict[\"energies\"][\"energy\"] = wfn.energy()\n")
                outf.write("forces    = grad.to_array()\n")
            if True:
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
                outf.write("    if None != forces and None != forces.all():\n")
                outf.write("        for f in forces:\n")
                outf.write("            result.write(\"%12.8f  %12.8f  %12.8f\\n\" % (f[0], f[1], f[2]))\n")

        return myjob
    
    def next_idir(self, idir:int)->str:
        my_idir = ("%04d" % idir)
        while os.path.isdir(my_idir):
            idir += 1
            my_idir = ("%04d" % idir)
        return idir, my_idir
    
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
                idir, my_idir = self.next_idir(idir)
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
                dist_vec = [ 0, 0, 1 ]
                
                # Loop over distances
                for dist_frac in df:
                    # Make directory
                    idir, my_idir = self.next_idir(idir)
                    os.makedirs(my_idir, exist_ok=False)
                    os.chdir(my_idir)

                    xyz1ori = orient(label1, xyz1com, rnd.uniform(0, Pi2), rnd.uniform(0, Pi2), rnd.uniform(0, Pi))
                    xyz2ori = orient(label2, xyz2com, rnd.uniform(0, Pi2), rnd.uniform(0, Pi2), rnd.uniform(0, Pi))
                    rel_dist = 0
                    while rel_dist < dist_frac:
                        # Shift all atoms in compound 2 by delta in the Z direction
                        for m in range(len(xyz2ori)):
                            xyz2ori[m][2] += delta
                        # Compute relative distance between the compounds
                        if self.verbose:
                            print("dist_frac %f dist %f" % ( dist_frac, xyz2ori[0][2]-xyz1ori[0][2] ))
                        new_ori_dist, new_dist_vec, new_ij = compute_mindist(xyz1ori, xyz2ori)
                        vdw_sum  = vdwdict[label1[new_ij[0]]] + vdwdict[label2[new_ij[1]]]
                        rel_dist = math.sqrt(new_dist_vec[0]**2+new_dist_vec[1]**2+new_dist_vec[2]**2)/vdw_sum
                        if self.verbose:
                            print("new_ori_dist %g vdw_sum %g rel_dist %g dist_frac %g" % ( new_ori_dist, vdw_sum, rel_dist, dist_frac ))

                    myjob = self.write_dimer_input(my_idir, userfile,
                                                   dimer["q1"], dimer["m1"], label1, xyz1ori,
                                                   dimer["q2"], dimer["m2"], label2, xyz2ori)
                    self.run_one_job(myjob)
                                  
                    os.chdir("..")
                    idir += 1
        os.chdir("..")
    
    def run_dimers(self, dimerfile:str, userfiles:bool,
                   vdwdict:dict, force:bool):
        if not os.path.exists(dimerfile):
            sys.exit("Selection file %s does not exist" % dimerfile)
        dimers  = get_dimer_selection(dimerfile)
        print("There are %d dimers in the selection %s" % ( len(dimers), dimerfile ))
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
                    self.run_one_dimer(d, pdbfile, vdwdict)
            else:
                self.run_one_dimer(d, "", vdwdict)
        os.chdir("../../")
    
def parse_args():
    desc = """
Run Psi4 calculations for use in ACT. When running dimers, the product of the ndist and norient and
the number of dimers in the selection file equals the total number of calculations. It can be a lot!
If the number of orientations is zero, the script will look for minimized structures for the dimers
and make a distance scan based on that. In this case the minimum and maximum distance are interpreted
as a fraction of the center of mass distance of the minimized input structure, e.g. from 0.9 to 1.5
times the center of mass distance. If norient is larger than zero, monomers will be randomly rotated and
a distance scan done along the z-coordinate. The mindist and maxdist given by the user are interpreted
as being relative to the sum of the Van der Waals radii of the closest atoms.
    """
    parser  = argparse.ArgumentParser(description=desc)
    defbasis = "aug-cc-pvdz"
    parser.add_argument("-basis", "--basis", help="Basis set, default "+defbasis, type=str, default=defbasis)
    defmethod = "sapt2+"
    parser.add_argument("-method","--method", help="QM method, default "+defmethod, type=str, default=defmethod)

    parser.add_argument("-fc","--frozen_core", help="Use frozen core orbitals, useful for expensive methods. Note that this is the default in other QM codes, such as Gaussian, but not in Psi4. If used, the basis set name will be extended with '-fc'", action="store_true")
    parser.add_argument("-v", "--verbose", help="Write debugging output", action="store_true")
    parser.add_argument("-dimers", "--dimers", help="Compute dimer interactions based on compounds in a selection file, please provide file name with this flag", type=str, default=None)
    parser.add_argument("-user", "--userfiles", help="Employ input coordinates supplied in the 'user' directory with name similar to those in xyz/dimers, namely mol1#mol2#number.pdb. You still need to provide the -dimer flag", action="store_true")
    ndist   = 1
    parser.add_argument("-ndist", "--ndist", help="Number of dimer distances to compute, see helptext. Default "+str(ndist), type=int, default=ndist)
    mindist = 2.5
    parser.add_argument("-mindist", "--mindist", help="Minimum dimer distances to compute, default "+str(mindist), type=float, default=mindist)
    maxdist = 5.0
    parser.add_argument("-maxdist", "--maxdist", help="Maximum dimer distances to compute, default "+str(maxdist), type=float, default=maxdist)
    norient = 1
    parser.add_argument("-norient", "--norient", help="Number of dimer orientation to compute, see help text. Default "+str(norient), type=int, default=norient)
    ncores = 16
    parser.add_argument("-ncores", "--ncores", help="Number of cores to use for each calculation, default "+str(ncores), type=int, default=ncores)
    memory = 3000
    parser.add_argument("-mem", "--memory", help="Memory per core to use for each calculation, default "+str(memory)+" Mb", type=int, default=memory)
    hours = 24
    parser.add_argument("-hours", "--hours", help="Amount of hours to queue for, default "+str(hours), type=int, default=hours)
    parser.add_argument("-dryrun", "--dryrun", help="Write scripts but do not perform the runs", action="store_true")
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
            if args.dimers:
                psi4.run_dimers(args.dimers, args.userfiles, vdwdict, False)
