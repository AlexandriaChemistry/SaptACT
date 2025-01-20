#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import datetime
import string,sys,csv,os
import getpass
import time
from pathlib     import Path

# ACT python code
from get_csv_rows import *
from get_mol_dict import *
from molutils import *

try:
    from pubchempy import *
except ImportError:
    print("Warning: pubchempy is not installed. Proceed at your own risk.")

wtol    = 0.5 # a.m.u.
debug   = 0

def alex_datetime():
    mydate = datetime.datetime.utcnow()
    return str(mydate)[0:-7]

def set_debug_level(level):
    global debug
    print("Setting debug to %d" % ( int(level) ) )
    debug = level

def get_debug_level():
    return debug

def mol_csv_die(reason):
    print("Death horror: %s" % ( reason ) )
    exit(1)

def filename_to_dir(filename):
    ddd = { 'zmat': "ZMAT", 'xyz': "XYZ", 'sdf': "SDF" }
    ext = filename[filename.rfind('.')+1:]
    if ext in ddd:
        return ddd[ext]
    print("Incorrect filename %s" % filename)
    exit(1)

def remove_charge_from_formula(formula):
    if (-1 != formula.find("+")):
        formula = formula.replace("+", " ").split()[0]
    elif (-1 != formula.find("-")):
        formula = formula.replace("-", " ").split()[0]
    return formula

def formula_equal(our, their, remove_charge = False):
    if (their.find('.') != -1):
        return True
    if (remove_charge):
        return (remove_charge_from_formula(our) ==
                remove_charge_from_formula(their))
    else:
        return (their == our)

def compare_formula(logfile, their, our):
    if (their.find('.') != -1):
        return
    tt = remove_charge_from_formula(their)
    oo = remove_charge_from_formula(our)
    if (tt != oo):
        fff_tt = parse_formula(tt)
        fff_oo = parse_formula(oo)
        for ft in fff_tt:
            if (not ft in fff_oo or (fff_tt[ft] != fff_oo[ft])):
                logfile.write("Formula mismatch expected '%s' pybel says '%s'\n" % ( oo, tt ) )
        for fo in fff_oo:
            if (not fo in fff_tt or (fff_tt[fo] != fff_oo[fo])):
                logfile.write("Formula mismatch expected '%s' pybel says '%s'\n" % ( oo, tt ) )

def compare_mol_formula(a, b):
    aa = remove_charge_from_formula(a)
    bb = remove_charge_from_formula(b)
    delta = 0
    if (aa != bb):
        fff_aa = parse_formula(aa)
        fff_bb = parse_formula(bb)
        for organic in [ "C", "H" ]:
            if organic in fff_aa:
                if organic in fff_bb:
                    if delta == 0:
                        delta = fff_aa[organic] - fff_bb[organic]
                else:
                    delta = 1
            elif organic in fff_bb:
                delta = -1
    return delta

def inchi_equal(their:str, our:str)->bool:
    if their == our:
        # All good
        return True
    else:
        return False
    if None == our:
        return False
    if None == their:
        return False
    # Check the layers. TODO: fix this to be correct.
    in1 = their.split('/')
    in2 = our.split('/')
    minlen = min(len(in1), len(in2))
    for k in range(minlen):
        if (not ((in1[k] == in2[k]) or
                 (k >= 3 and ((in1[k][0] in [ 'p', 'q' ]) and
                              (in2[k][0] in [ 'p', 'q' ]) and
                              (in1[k][1:] == in2[k][1:]))))):
            return False
    return True
    
def compare_inchi(logfile, their, our):
    if (their == our):
        # All good
        return

    # Check the layers
    in1 = their.split('/')
    in2 = our.split('/')
    minlen = min(len(in1), len(in2))
    for k in range(minlen):
        if (not ((in1[k] == in2[k]) or
                 (k >= 3 and ((in1[k][0] in [ 'p', 'q' ]) and
                              (in2[k][0] in [ 'p', 'q' ]) and
                              (in1[k][1:] == in2[k][1:]))))):
            logfile.write("InChI mismatch %d ours '%s' pybel '%s'\n" % ( k, our, their ) )
            return
    # So far so good, first bit is the same, now check the remaining
    if (len(in1) > len(in2)):
        logfile.write("Missing layers for our InChI ")
        for k in range(len(in2), len(in1)):
            logfile.write("/%s" % in1[k])
        logfile.write("\n")
    elif (len(in1) < len(in2)):
        logfile.write("Missing layers in pybel InChI ")
        for k in range(len(in1), len(in2)):
            logfile.write("/%s" % in2[k])
        logfile.write("\n")
    else:
        logfile.write("InChI mismatch in layers %d-%d " % ( minlen, len(in1) ) )
        for k in range(minlen, len(in1)):
            logfile.write("our: %s their %s" % ( in2[k], in1[k] ) )
        logfile.write("\n")

def set_sdf_charge(sdf_fn, charge):
    tmp = sdf_fn + ".bak"
    try:
        os.rename(sdf_fn, tmp)
    except:
        print("Can not rename file %s to %s" % ( sdf_fn, tmp ) )
        return
    out = open(sdf_fn, "w")
    try:
        sdf = open(tmp, "r")
        found_chg = False
        for line in sdf:
            if (line.find("M  CHG") >= 0):
                found_chg = True
            if (line.find("M  END") >= 0 and not found_chg):
                out.write("M  CHG  1   1  %2d\n" % charge)
            out.write("%s" % line)
    finally:
        sdf.close()
    out.close()

def supported_elements() -> str:
    return ['H','He', 'C', 'N', 'O', 'F', 'Ne', 'P', 'S', 'Cl', 'Ar', 'Br', 'Kr', 'I']

class Molecule:
    """Class to handle an entry from the alexandria.csv database"""

    def __init__(self, alexandria:int):
        self.alexandria = alexandria

    def check_html(self):
        if len(self.html_name) == 0:
            self.html_name = self.iupac

    def check_latex(self):
        if (get_debug_level() > 0):
            if (self.html_name.find('$') != -1):
                print("LaTeX codes in compound %d html name '%s'" % ( self.alexandria, self.html_name ) )
            if (self.iupac.find('$') != -1):
                print("LaTeX codes in compound %d iupac name '%s'" % ( self.alexandria, self.iupac ) )
        if (self.iupac.find("[") >= 0 or self.iupac.find("]") >= 0) and (len(self.latex_name) == 0):
            self.latex_name = self.iupac.replace("[", "{[}").replace("]", "{]}")

    def from_row(self, row, version:int):
        if (version != 3):
            sys.exit("Version %d not supported" % version)

        if (len(row) < 22):
            sys.exit("row[0] %s invalid input" % ( row[0] ))
        if self.alexandria <= 0:
            self.alexandria = int(row[0])
        self.parentid        = int(row[1])
        self.iupac           = row[2]
        self.formula         = row[3]
        self.charge          = int(row[4])
        self.mult            = int(row[5])
        self.weight          = row[6]
        self.cas             = row[7]
        self.csid            = row[8]
        self.pubchem         = row[9]
        self.stdinchi        = row[10]
        self.inchikey        = row[11]
        self.classification  = row[12]
        self.synonyms        = row[13]
        self.filename        = row[14]
        self.point_group     = row[15]
        self.symmetry_number = row[16]
        self.rotbond         = row[17]
        self.html_name       = row[18]
        self.latex_name      = row[19]
        self.author          = row[20]
        self.timestamp       = row[21]
        self.comment         = None
        if len(row) >= 23:
            self.comment     = row[22]
        self.check_latex()
        self.check_html()

    def classes(self):
        c = []

        cc = self.classification.strip()
        for s in cc.split(";"):
            ss = s.strip()
            if (len(ss) > 0) and ss not in c:
                c.append(ss)
        return c

    def names(self):
        n = [ self.iupac ]

        if (None != self.cas and len(self.cas) > 0):
            n.append(str(self.cas))
        if (None != self.stdinchi and len(self.stdinchi) > 0):
            n.append(str(self.stdinchi))
        syns = self.synonyms.strip()
        if (len(syns) > 0):
            for s in syns.split(";"):
                if (len(s) > 0):
                    n.append(str(s.strip()))
        n.append(Path(self.filename).stem)
        return n

    def add_classification(self, classification):

        if (classification and len(classification) >0):
            clss = self.classes()
            if (not classification in clss):
                if (len(self.classification) == 0):
                    self.classification = classification
                else:
                    self.classification += ";" + classification
            # always replace all spaces with underscores
            self.classification = self.classification.replace(" ","_")
            # remove duplicates (https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists)
            cc = sorted(list(dict.fromkeys(self.classification.split(";"))))
            self.classification = ';'.join(str(x) for x in cc)

    def clean_classification(self):
        # Clean classification
        newclass = None
        for c in self.classes():
            if (newclass == None):
                newclass = c
            else:
                newclass = newclass + ';' + c
        if (newclass):
            self.classification = newclass

    def add_synonym(self, newsynonym):
        self.synonyms += (";%s" % newsynonym)

    def clean_synonyms(self):
        # Clean synonyms
        newsyn = None
        for c in self.names():
            if (c != self.iupac and
                self.filename != (c + '.zmat') and
                self.filename != (c + '.xyz') and
                self.cas != c and
                self.stdinchi != c and
                self.csid != c):
                if (newsyn == None):
                    newsyn = c
                else:
                    newsyn = newsyn + ';' + c
        if (newsyn):
            self.synonyms = newsyn

    def clean(self):
        self.clean_synonyms()
        self.clean_classification()

    def print_self(self, logfile):
        logfile.write("Self-checking compound %d %s fn %s\n" % ( self.alexandria, self.iupac, self.filename ) )

    def add_ion_categories(self):
        if (self.formula.find('+') != -1 or self.charge > 0):
            self.classification += ";cation"
            if (not (self.formula.find('+') != -1 and self.charge > 0)):
                print("Inconsistency: formula %s charge %d filename %s" % (self.formula, self.charge, self.filename))
        if (self.formula.find('-') != -1 or self.charge < 0):
            self.classification += ";anion"
            if (not (self.formula.find('-') != -1 and self.charge < 0)):
                print("Inconsistency: formula %s charge %d filename %s" % (self.formula, self.charge, self.filename))
        if ((int(self.mult) % 2) == 0):
            self.classification += ";radical"

    def myprint(self):
        print("alexandria_id   %d" % ( self.alexandria ))
        print("parent_id       %d" % ( self.parentid ))
        print("iupac           %s" % ( self.iupac ))
        print("formula         %s" % ( self.formula ))
        print("charge          %d" % ( self.charge ))
        print("multiplicity    %d" % ( self.mult ))
        print("weight          %s" % ( self.weight ))
        print("cas             %s" % ( self.cas ))
        print("chemspider_id   %s" % ( self.csid ))
        print("pubchem_id      %s" % ( self.pubchem ))
        print("stdinchi        %s" % ( self.stdinchi ))
        print("inchikey        %s" % ( self.inchikey ))
        print("classification  %s" % ( self.classification ))
        print("synonyms        %s" % ( self.synonyms ))
        print("filename        %s" % ( self.filename ))
        print("point_group     %s" % ( self.point_group ))
        print("symmetry_number %s" % ( self.symmetry_number ))
        print("rotbond         %s" % ( self.rotbond ))
        print("html_name       %s" % ( self.html_name ))
        print("latex_name      %s" % ( self.latex_name ))
        print("author          %s" % ( self.author ))
        print("timestamp       %s" % ( self.timestamp ))
        print("comment         %s" % ( self.comment ))

    def is_comment(self):
        return (self.iupac[0] == '#')

    def set_pubchem_id(self):
        for compound in get_compounds(self.iupac, 'name'):
            if (compound.molecular_formula == self.formula and
                compound.iupac_name == self.iupac and
                abs(float(compound.molecular_weight)-float(self.weight)) < wtol and
                compound.charge == self.charge):
                self.pubchem  = str(compound.cid)
                self.inchikey = str(compound.inchikey)

    def row(self, version):
        delim = '|'
        if (version == 1):
            myrow = self.iupac + delim + self.formula + delim + self.weight + delim + self.cas + delim + self.csid + delim
            myrow += self.stdinchi + delim + self.inchikey + delim
            myrow += self.classification + delim + self.synonyms + delim + self.filename + delim + self.point_group + delim + self.symmetry_number + delim + self.rotbond + delim + self.html_name
        elif (version == 2):
            if (None == self.alexandria):
                return None
            myrow = str(self.alexandria) + delim + self.iupac + delim + self.formula + delim
            myrow += str(self.charge) + delim + str(self.mult) + delim + self.weight + delim + self.cas + delim + self.csid + delim
            myrow += self.pubchem + delim + self.stdinchi + delim + self.inchikey + delim
            myrow += self.classification + delim + self.synonyms + delim + self.filename + delim + self.point_group + delim + self.symmetry_number + delim + self.rotbond + delim + self.html_name + delim
            myrow += self.latex_name + delim + self.author + delim + self.timestamp + delim + self.comment
        elif (version == 3):
            if (None == self.alexandria):
                return None
            myrow = str(self.alexandria) + delim + str(self.parentid) + delim + self.iupac + delim + self.formula + delim
            myrow += str(self.charge) + delim + str(self.mult) + delim + self.weight + delim + self.cas + delim + self.csid + delim
            myrow += self.pubchem + delim + self.stdinchi + delim + self.inchikey + delim
            myrow += self.classification + delim + self.synonyms + delim + self.filename + delim + self.point_group + delim + self.symmetry_number + delim + self.rotbond + delim
            if self.html_name != self.iupac:
                myrow += self.html_name
            myrow += delim
            if self.latex_name != self.iupac:
                myrow += self.latex_name
            myrow += delim + self.author + delim + self.timestamp
            if self.comment:
                myrow += delim + self.comment
        else:
            mol_csv_die("version")
        return myrow

    def unstamp(self):
        self.author     = ""
        self.timestamp  = ""

    def stamp(self):
        self.author     = getpass.getuser()
        self.timestamp  = alex_datetime()

    def get_stamp(self):
        return self.author + " " + self.timestamp

    def stamped(self):
        if self.author and self.timestamp:
            return True
        else:
            return False

    def get_iupac(self):
        return self.iupac

    def is_supported(self):
        elements = parse_formula(self.formula)
        if all(element in supported_elements() for element in elements.keys()) and (self.mult % 2 == 1):
            return True
        return False

    def is_excluded(self):
        needle = self.filename.split('.')[0]
        excluded_file = "data/excluded_molecules.csv"
        for row in get_csv_rows(excluded_file, 2):
            if len(row) == 2:
                if row[0] == needle:
                   return True
        return False

    # Couple of functions for sorting on formula
    def __lt__(self, other):
        return compare_mol_formula(self.formula, other.formula) < 0
    def __gt__(self, other):
        return compare_mol_formula(self.formula, other.formula) > 0
    def __eq__(self, other):
        return compare_mol_formula(self.formula, other.formula) == 0
    def __le__(self, other):
        return compare_mol_formula(self.formula, other.formula) <= 0
    def __ge__(self, other):
        return compare_mol_formula(self.formula, other.formula) >= 0
    def __ne__(self, other):
        if None == other:
           return True
        return compare_mol_formula(self.formula, other.formula) != 0

class Molecules:
    """Class to handle the alexandria.csv database"""
    def __init__(self):
        self.mols            = {}
        self.name_dict       = {}
        self.class_dict      = {}
        self.pubchem_dict    = {}
        self.chemspider_dict = {}

    def make_name_dict(self):
        for m in self.mols:
            mol = self.mols[m]
            if (not mol.is_comment()):
                for name in mol.names():
                    if ((name in self.name_dict) and
                        (self.name_dict[name] in self.mols) and
                        (mol.iupac != self.mols[self.name_dict[name]].iupac) ):
                        if (get_debug_level() > 0):
                            print("Identifier '%s' is used for '%s', will be ignored for '%s'" %
                                (name, mol.iupac, self.mols[self.name_dict[name]].iupac ) )
                    else:
                        self.name_dict[name] = str(mol.alexandria)

    def make_class_dict(self):
        for m in self.mols:
            mol = self.mols[m]
            if (not mol.is_comment()):
                for ccc in mol.classes():
                    if (ccc in self.class_dict):
                        self.class_dict[ccc] = (self.class_dict[ccc] + ";" + mol.iupac)
                    else:
                        self.class_dict[ccc] = mol.iupac

    def make_id_dict(self):
        for m in self.mols:
            mol = self.mols[m]
            csid = mol.csid
            if (len(csid) > 0):
                self.chemspider_dict[str(csid)] = mol.alexandria
            pcid = mol.pubchem
            if (len(pcid) > 0):
                self.pubchem_dict[str(pcid)] = mol.alexandria

    def fetch_pubchem(self, pubchem):
        if (str(pubchem) in self.pubchem_dict):
            return self.mols[self.pubchem_dict[str(pubchem)]]
        return None

    def fetch_chemspider(self, csid):
        if (str(csid) in self.chemspider_dict):
            return self.mols[self.chemspider_dict[str(csid)]]
        return None

    def next_id(self):
        max_id = 0
        for m in self.mols:
            max_id = max(max_id, self.mols[m].alexandria)
        return 1+max_id

    def delete_molecule(self, mol_id):
        self.mols.pop(mol_id, None)

    def add_molecule(self, mol):
        if (mol.alexandria in self.mols):
            print("Consistency error: alexandria ID %d already in the database. Will generate new id." % ( mol.alexandria ))
            mol.alexandria = self.next_id()
        self.mols[mol.alexandria] = mol

    def read(self, input_file:str, version, clean:bool):
        rows    = get_csv_rows(input_file, 18)
        for row in rows:
            # Default mol_id = 0 means read from file.
            mol = Molecule(0)
            mol.from_row(row, version)
            self.add_molecule(mol)

        self.make_name_dict()
        self.make_class_dict()
        self.make_id_dict()
        if (clean):
            self.clean()

    def read_default(self):
        acsv = "data/alexandria.csv"
        self.read(acsv, 3, False)

    def write(self, output_file:str, version:int):
        outputfile = open(output_file, "w", encoding='utf-8')

        try:
            for mol in sorted(self.mols):
                row = self.mols[mol].row(version)
                if (row):
                    outputfile.write(row + "\n")

        finally:
            outputfile.close()

    def set_pubchem_id(self):
        for mol in self.mols:
            time.sleep(1)
            mol.set_pubchem_id()

    def add_qmult(self, qmcsv):
        inputfile = open(qmcsv, "r", encoding='utf-8')

        csv.register_dialect('pipes', delimiter=delim)

        try:
            reader = csv.reader(inputfile, dialect='pipes')

            for row in reader:
                if (row[0] in self.name_dict):
                    mid = self.name_dict[row[0]]
                    self.mols[mid].charge = int(row[1])
                    self.mols[mid].mult   = int(row[2])

        finally:
            inputfile.close()

    def add_ion_categories(self):
        nclean_cls = 0
        for m in self.mols:
            mol = self.mols[m]
            cls = mol.classification
            mol.add_ion_categories()
            if (cls != mol.classification):
                nclean_cls += 1
        print("Just added ion classed to %d mols" % ( nclean_cls ) )

    def clean(self):
        nclean_syn = 0
        nclean_cls = 0
        for m in self.mols:
            mol = self.mols[m]
            syn = mol.synonyms
            cls = mol.classification
            mol.clean()
            if (syn != mol.synonyms):
                nclean_syn += 1
            if (cls != mol.classification):
                nclean_cls += 1
        print("Just cleaned %d synonyms and %d classifications" % ( nclean_syn, nclean_cls ) )

    def use_sdf(self):
        nupdate   = 0
        nnotfound = 0
        for m in self.mols:
            mol = self.mols[m]
            sdf = None
            if (mol.filename[-4:] == ".xyz"):
                sdf = mol.filename[:-4] + ".sdf"
            elif (mol.filename[-5:] == ".zmat"):
                sdf = mol.filename[:-5] + ".sdf"
            elif (mol.filename[-4:] != ".sdf"):
                print("Incorrect filename %s" % mol.filename)
                continue
            if (None != sdf):
                if (os.path.isfile("../MOLECULES/SDF/" + sdf)):
                    mol.filename = sdf
                    nupdate += 1
                else:
                    print("%s not found" % sdf)
                    nnotfound += 1
        print("For %d compounds the filenames were updated and entries unstamped." % nupdate )
        print("For %d compounds no sdf files were found." % nnotfound )

    def find_mol(self, molname):
        if (molname in self.name_dict):
            molindex = self.name_dict[str(molname)]
            if (debug):
                print("Found %s in name_dict, value %s" %
                    ( molname, molindex ) )
            if (int(molindex) in self.mols):
                if (debug):
                    print("Found %s in mols" % ( self.name_dict[molname] ) )
                return self.mols[int(molindex)]
        elif debug:
            print("%s not found." % molname)
        return None

    def find_iupac(self, molname):
        m = self.find_mol(molname)
        if m:
            return m.iupac
        else:
            return None

    def find_pubchem(self, pubchem):
        if (debug):
            print("Looking for %d" % pubchem)
        for m in self.mols:
            if self.mols[m].pubchem == pubchem:
                return self.mols[m]
        return None

    def find_inchi(self, inchi:str):
        if (debug):
            print("Looking for %s" % inchi)
        for m in self.mols:
            if None == self.mols[m].stdinchi:
               print("No inchi in %s" % ( self.mols[m].iupac ))
            elif inchi_equal(self.mols[m].stdinchi, inchi):
                return self.mols[m]
        return None

    def find_class(self, mclass):
        if (mclass in self.class_dict):
            return self.class_dict[mclass].split(";")
        else:
            return []

    def dump(self):
        for name in self.dict:
            print("%s|%s" % ( name, self.dict[name].iupac ) )

    def nmols(self):
        return len(self.mols)

# Utility to prettify a formula
def to_latex_formula(formula):
    myq = ""
    if formula.find("-"):
        w = formula.split("-")
        formula = w[0]
        if len(w) > 1:
            myq     = ("$^{-%s}$" % w[1])
    elif formula.find("+"):
        w = formula.split("+")
        formula = w[0]
        if len(w) > 1:
            myq     = ("$^{+%s}$" % w[1])
    
    for i in range(10):
        istr = str(i)
        jstr = ("$_%d$" % i)
        formula = formula.replace(istr, jstr)
    return formula + myq

