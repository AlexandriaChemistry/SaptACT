#!/usr/bin/env python3

import argparse, copy, os, json, glob, math, sys, xmltodict
from run_calcs import add_lot_args, get_dimer_selection

Hartree = 2625.5
minstr  = "min"
maxstr  = "max"

def parse():
    desc = "Generate plots of dimer energies as a function of distance."
    parser  = argparse.ArgumentParser(description=desc)
    add_lot_args(parser)
    parser.add_argument("-sel", "--selection", help="Plot dimer interactions based on compounds in a selection file, please provide file name with this flag. Default is to plot all dimers.", type=str, default=None)
    parser.add_argument("-mp", "--molprop", help="Use molprop file as input instead of json files", type=str, default=None)
    parser.add_argument("-cap","--caption", help="Extra text to add to the captoions", type=str, default="")
    dEmax = 1
    parser.add_argument("-dEmax", "--dEmax", help="Highest exchange energy to include. Ignored when using a molprop file. Default "+str(dEmax)+" Hartree", type=float, default=dEmax)
    parser.add_argument("-newdist", "--newdist", help="Generate csv file with minimum and maximum distance between monomers", type=str, default=None)

    args = parser.parse_args()
    return args

def doplot(xvgf:str, pdff:str):
    os.system("viewxvg -f %s -legend_x 0.5 -ls None -mk o x '*' + v  -pdf %s -noshow -tickfs 24 -lfs 30 -alfs 30" % ( xvgf, pdff ) )

class NewDist:
    def __init__(self):
        self.nd = {}

    def add(self, dimer:str, dist:float):
        if not dimer in self.nd:
            self.nd[dimer] = { minstr: dist, maxstr: dist, "N": 1 }
        else:
            self.nd[dimer][minstr] = min(dist, self.nd[dimer][minstr])
            self.nd[dimer][maxstr] = max(dist, self.nd[dimer][maxstr])
            self.nd[dimer]["N"]  += 1

    def save(self, filenm:str):
        with open(filenm, "w") as outf:
            outf.write("#Dimer,Min,Max,N\n")
            for dim in self.nd:
                outf.write("%s,%f,%f,%d\n" %
                           ( dim,
                             self.nd[dim][minstr],
                             self.nd[dim][maxstr],
                             self.nd[dim]["N"]))

def json2pdfs(dimsel:list, sapt:list, destdir:str, newdist:dict)->list:
    pdfs = []
    for mydir in dimsel:
        if not os.path.isdir(mydir):
            continue
        os.chdir(mydir)
        confs = []
        for index in glob.glob("0*"):
            if not os.path.isdir(index):
                continue
            os.chdir(index)
            jsonfile = "results.json"
            if os.path.exists(jsonfile):
                rawdata =  open(jsonfile, "r")
                confs.append(json.load(rawdata))
                rawdata.close()
            os.chdir("..")
        curve = []
        for c in confs:
            eee = "energies"
            iee = "Exchange"
            if not (eee in c and iee in c[eee] and c[eee][iee] < args.dEmax):
                continue
            if len(c["mols"]) == 2 and eee in c:
                # Dimer found
                mindist = 1e8
                for ix in c["mols"][0]["atoms"]:
                    for jx in c["mols"][1]["atoms"]:
                        dist2 = 0
                        for m in range(3):
                            dist2 += (ix["coords"][m] - jx["coords"][m])**2
                        mindist = min(mindist, math.sqrt(dist2))
                cmap = {}
                for term in sapt:
                    if term in c[eee]:
                        cmap[term] = float(c[eee][term])*Hartree
                if len(cmap.keys()) == len(sapt):
                    curve.append( ( mindist, cmap ) )
                    if newdist:
                        newdist.add(mydir, mindist)
        print("Found %d calcs in %s, extracted %d points" % ( len(confs), mydir, len(curve) ) )
        if len(curve) > 0:
            curve.sort(key=lambda a: a[0])
            xvgf = mydir + ".xvg"
            Npoints = 0
            with open(xvgf, "w") as outf:
                outf.write("@ title \"%s (N = %d)\"\n" % ( mydir.replace('#', '-' ), len(curve) ) )
                outf.write("@ xaxis label \"Minimum distance ($\\mathrm{\\AA}$)\"\n")
                outf.write("@ yaxis label \"Energy (kJ/mol)\"\n")
                iii = 0
                for ss in sapt:
                    outf.write("@ s%d legend \"%s\"\n" % ( iii, ss ))
                    iii += 1
                for kk in curve:
                    ck = kk[1]
                    outf.write("%10g" % ( kk[0] ) )
                    for ss in ck:
                        outf.write("  %10g" % ck[ss])
                    outf.write("\n")
            if os.path.exists(xvgf):
                mypdf = mydir.replace('#', '-') + ".pdf"
                pdfs.append(mypdf)
                pdff = destdir + "/" + mypdf
                doplot(xvgf, pdff)
        os.chdir("..")
    return pdfs

def calc_dist(frags:list, atoms:list)->float:
    mindist = None
    if len(frags) == 2:
        for i in frags[0].split():
            for j in frags[1].split():
                d2 = 0
                for c in [ "x", "y", "z" ]:
                    d2 += (float(atoms[int(i)-1][c]) - 
                           float(atoms[int(j)-1][c]))**2
                if not mindist:
                    mindist = math.sqrt(d2)
                else:
                    mindist = min(mindist, math.sqrt(d2))
    else:
        print(frags)
    return mindist
    
def molprop2pdfs(molprop:str, sapt:list, destdir:str, newdist:dict)->list:
    pdfs = []
    verbose = False
    debug   = False
    with open(molprop) as fd:
        doc = xmltodict.parse(fd.read())

        mols = "molecules"
        if not mols in doc:
            sys.exit("No %s in %s" % ( mols, infile ))
        mol = "molecule"
        nmols = len(doc[mols][mol])
        print("Number of molecules %d" % nmols)
        exp  = "experiment"
        frag = "fragments"
        atom = "atom"
        molatoms = {}
        for m in range(nmols):
            molname = doc[mols][mol][m]["@molname"]
            nelem   = len(doc[mols][mol][m])
            myfrags = []
            plotdict = []
            for k in doc[mols][mol][m].keys():
                if debug:
                    print(doc[mols][mol][m])
                    print("k = {}".format(k))
                if frag == k:
                    if debug:
                        print("FF {}".format(doc[mols][mol][m][k]))
                    for ff in range(len(doc[mols][mol][m][k]['fragment'])):
                        myfrags.append(doc[mols][mol][m][k]['fragment'][ff]['#text'])
                        if debug:
                            print("Found fragment {}".format(myfrags[len(myfrags)-1]))
                elif exp == k:
                    myexps = doc[mols][mol][m][k]
                    if isinstance(myexps, dict):
                        myexps = [ myexps ]
                    for me in myexps:
                        if debug:
                            print("me: {}".format(me))
                        atoms = []
                        eners = {}
                        for n in me:
                            if debug:
                                print("We have {}".format(n))
                            if n == 'atom':
                                for a in me[n]:
                                    if debug:
                                        print("a {}".format(a['x']))
                                    atoms.append({ "@atomid": a["@atomid"], "@name": a["@name"], "x": a['x'], "y": a['y'], "z": a['z'] })
                            elif n == 'energy':
                                for ie in range(len(me[n])):
                                    etype = me[n][ie]["@type"]
                                    evalu = float(me[n][ie]["average"])
                                    if debug:
                                        print("%s %s" % ( etype, evalu ))
                                    eners[etype] = evalu
                            else:
                                continue
                        if not molname in molatoms:
                            molatoms[molname] = copy.deepcopy(atoms)
                            if verbose:
                                print("Found %d atoms for %s %s" % ( len(atoms), molname, str(atoms)))
                        if debug:
                            print(atoms)
                        dist = calc_dist(myfrags, atoms)
                        if newdist:
                            newdist.add(molname, dist)
                        if debug:
                            print("Distance %s" % str(dist))
                        plotdict.append( { "dist": dist, "ener": eners } )
            dxvg = destdir + "/" + molname + ".xvg"
            eleg = [ "Dispersion", "Exchange", "Electrostatics",
                     "Induction", "InteractionEnergy" ]

            with open(dxvg, "w") as outf:
                outf.write("@ xaxis label \"Distance ($\\mathrm{\\AA}$)\"\n")
                outf.write("@ yaxis label \"Energy (kJ/mol)\"\n")
                for k in range(len(eleg)):
                    outf.write("@ s%d legend \"%s\"\n" % ( k, eleg[k] ))
                for k in sorted(plotdict, key=lambda pd : pd["dist"]):
                    outf.write("%10g" % k["dist"])
                    for n in range(len(eleg)):
                        outf.write("  %12g" % ( k["ener"][eleg[n]] ) )
                    outf.write("\n")
            pdff = molname.replace("#", "-") + ".pdf"
            doplot(dxvg, destdir + "/" + pdff)
            pdfs.append(pdff)

    return pdfs

if __name__ == "__main__":
    args = parse()
    sapt = [ "Dispersion", "Exchange", "Electrostatics", "Induction", "InteractionEnergy" ]
    outdir = 'SaptPlots'
    destdir = os.getcwd() + "/" + outdir
    os.makedirs(destdir, exist_ok=True)

    firstlot = None
    lots = []
    for m in args.method:
        for b in args.basis:
            if not firstlot:
                firstlot = m + "/" + b
            lots.append(m + "-" + b)
    if len(lots) != 1:
        sys.exit("Please specify one method and one basisset")
    lot = lots[0]

    newdist = None
    if args.newdist:
        newdist = NewDist()
    if args.molprop:
        print("Will use %s to read dimers and energies from" % args.molprop)
        pdfs = molprop2pdfs(args.molprop, sapt, destdir, newdist)
    else:
        print("Will use the '%s' level of theory" % lot)
        os.chdir("%s/dimer-scans" % lot)
        if args.selection:
            dimsel = []
            for d in get_dimer_selection(args.selection):
                dimsel.append(d["pair"])
        else:
            dimsel = sorted(glob.glob("*"))
        print("There are %d dimers in the selection" % len(dimsel))
        pdfs = json2pdfs(dimsel, sapt, destdir, newdist)

    os.chdir(destdir)
    texfn = outdir + ".tex"
    pdffn = outdir + ".pdf"
    with open(texfn, "w") as outf:
        outf.write("""\\documentclass[11pt]{article}

\\usepackage{sectsty,hyperref}
\\usepackage{graphicx}

% Margins
\\topmargin=-0.45in
\\evensidemargin=0in
\\oddsidemargin=0in
\\textwidth=6.5in
\\textheight=9.0in
\\headsep=0.25in

\\title{SAPT Data}
\\author{The Alexandria Chemistry Initiative}
\\date{\\today}

\\begin{document}
\\maketitle	
\\pagebreak
\\listoffigures
\\pagebreak
%\\section*{Plots}
""")

        lot = firstlot.replace("dmp2", "$\\delta$MP2").replace("ccd", "CCD")
        n = 0
        for pdf in pdfs:
            if not pdf == pdffn:
                outf.write("\\begin{figure}[ht]\n")
                outf.write("\\includegraphics[width=14cm]{%s}\n" % pdf)

                outf.write("\\caption{Energy components at the %s level of theory as a function of distance for %s. %s}\n" % ( lot, pdf[:-4], args.caption) )
                outf.write("\\end{figure}\n")
                if n % 30 == 0:
                    outf.write("\\cleardoublepage\n")
                n += 1
        outf.write("\\end{document}\n")
    
    for m in range(3):
        os.system("pdflatex %s" % texfn)
    print("Please check %s/%s" % ( outdir, pdffn ) )

    if newdist:
        os.chdir("..")
        newdist.save(args.newdist)
        print("Please check %s" % args.newdist)
