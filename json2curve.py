#!/usr/bin/env python3

import argparse, os, json, glob, math, sys
from run_calcs import get_dimer_selection

Hartree = 2625.5

def parse():
    desc = "Generate plots of dimer energies as a function of distance"
    parser  = argparse.ArgumentParser(description=desc)

    defbasis = "aug-cc-pvdz"
    parser.add_argument("-basis", "--basis", help="Basis set, default "+defbasis, type=str, default=defbasis)
    defmethod = "sapt2+"
    parser.add_argument("-method","--method", help="QM method, default "+defmethod, type=str, default=defmethod)

    parser.add_argument("-sel", "--selection", help="Plot dimer interactions based on compounds in a selection file, please provide file name with this flag. Default is to plot all dimers.", type=str, default=None)
    dEmax = 0.1
    parser.add_argument("-dEmax", "--dEmax", help="Highest exchange energy to include. Default "+str(dEmax)+" Hartree", type=float, default=dEmax)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse()
    lots = []
    for m in args.method:
        for b in args.basis:
            lots.append(m + "-" + b)
    if len(lots) != 1:
        sys.exit("Please specify one method and one basisset")
    lot = lots[0]
    print("Will use the '%s' level of theory" % lot)
    sapt = [ "Dispersion", "Exchange", "Electrostatics", "Induction", "InteractionEnergy" ]
    outdir = 'SaptPlots'
    destdir = os.getcwd() + "/" + outdir
    os.makedirs(destdir, exist_ok=True)
    if args.selection:
        dimsel = []
        for d in get_dimer_selection(args.selection):
            dimsel.append(d["pair"])
    else:
        dimsel = sorted(glob.glob("*"))
    
    os.chdir("%s/dimer-scans" % lot)
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
        print("Found %d calcs in %s, extracted %d points" % ( len(confs), mydir, len(curve) ) )
        if len(curve) > 0:
            curve.sort(key=lambda a: a[0])
            xvgf = mydir + ".xvg"
            Npoints = 0
            with open(xvgf, "w") as outf:
                outf.write("@ title \"%s (N = %d)\"\n" % ( mydir.replace('#', '-' ), len(curve) ) )
                outf.write("@ xaxis label \"Minimum distance (Angstrom)\"\n")
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
                os.system("viewxvg -f %s -legend_x 0.5 -ls None -mk -pdf %s -noshow" % ( xvgf, pdff ) )
        os.chdir("..")
    
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

        n = 0
        for pdf in pdfs:
            if not pdf == pdffn:
                outf.write("\\begin{figure}[ht]\n")
                outf.write("\\includegraphics[width=14cm]{%s}\n" % pdf)
                outf.write("\\caption{%s}\n" % pdf)
                outf.write("\\end{figure}\n")
                if n % 30 == 0:
                    outf.write("\\cleardoublepage\n")
                n += 1
        outf.write("\\end{document}\n")
    
    for m in range(3):
        os.system("pdflatex %s" % texfn)
    print("Please check %s/%s" % ( outdir, pdffn ) )


    
