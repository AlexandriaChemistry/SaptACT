#!/usr/bin/env python3

import argparse, json, math, os, sys, copy

HARTREE = 2625

def write_atom(outf, atom:int, resnr:int, name:str, x:float, y:float, z:float):
    outf.write("HETATM %4d %2s   UNK     %d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s\n" % (atom, name, resnr, x, y, z, 0, name))

def write_pdb(mydict:dict, defdist:float, outfn:str):
    move = [ 0, 0, 0 ]
    if defdist > 0:
        cogs = []
        for i in range(len(mydict["mols"])):
            imol = mydict["mols"][i]
            for iatom in imol["atoms"]:
                cog = [ 0, 0, 0 ]
                for m in range(len(cog)):
                    cog[m] += iatom["coords"][m]/len(imol["atoms"])
            cogs.append(cog)
        dist  = [ 0, 0, 0 ]
        dist2 = 0
        for m in range(3):
            dist[m] = cogs[0][m]-cogs[1][m]
            dist2  += dist[m]**2
        if dist2 == 0:
            move = [ 0, 0, defdist ]
        else:
            dist  = math.sqrt(dist2)
            scale = (defdist-dist)/dist
            for m in range(3):
                move[m] = scale*dist
    with open(outfn, "w") as outf:
        atom  = 1
        resnr = 1
        for i in range(len(mydict["mols"])):
            imol = mydict["mols"][i]
            for iatom in imol["atoms"]:
                write_atom(outf, atom, resnr, iatom["elem"],
                           iatom["coords"][0]+i*move[0],
                           iatom["coords"][1]+i*move[1],
                           iatom["coords"][2]+i*move[2])
                atom += 1
            resnr += 1
        outf.write("END\n")

def parse_args():
    desc = "Convert out file from run_calcs for dimer to pdb file. If multiple files are given, the one with the lowest energy will be written."
    parser  = argparse.ArgumentParser(description=desc)
    parser.add_argument("-f","--infile", nargs="+", help="Input file, type results.json", type=str, default=None)
    defout = "dimer.pdb"
    parser.add_argument("-o","--outfile", help="Output file, default "+defout, type=str, default=defout)
    parser.add_argument("-minimize","--minimize", help="Minimize energy of structure", action="store_true")
    defdist = 0
    parser.add_argument("-dist", "--dist", help="Set distance between monomers to this value, if zero (default) do nothing", type=float, default=defdist)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    if args.infile and args.outfile:
        mindict = None
        emin    = 1e8
        minfile = None
        for infile in args.infile:
            with open(infile, "r") as inf:
                mydict = json.load(inf)
            energies = "energies"
            iener    = "InteractionEnergy"
            if energies in mydict and iener in mydict[energies]:
                if mydict[energies][iener] < emin:
                    mindict = copy.deepcopy(mydict)
                    emin    = mydict[energies][iener]
                    minfile = infile
        temp = "temp.pdb"
        if mindict:
            write_pdb(mindict, args.dist, temp)
        if args.minimize:
            os.system("obminimize %s > %s" % ( temp, args.outfile ))
        else:
            os.system(f"mv {temp} {args.outfile}")
        print("Lowest energy %g kJ/mol file: %s" % ( emin*HARTREE, minfile ))
