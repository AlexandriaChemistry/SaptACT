#!/usr/bin/env python3

import argparse, os, sys
import random

def parse_args():
    desc = """
Generate list of dimers from a string of monomers.
    """
    parser  = argparse.ArgumentParser(description=desc)
    parser.add_argument("-sel", "--selection", nargs="+", help="List of monomers", type=str)
    parser.add_argument("-sel2", "--selection2", nargs="+", help="Second list of monomers,if present only dimers with one compound in the first and one compound in the second selection will be generated", type=str)
    output = "selection.dat"
    parser.add_argument("-o", "--output", help="Output file name, default "+output, type=str, default=output)
    parser.add_argument("-act", "--act", help="Make selection file for ACT with random Train and Test select", action="store_true")
    tf = 0.3
    parser.add_argument("-tf", "--test_fraction", help="Fraction test compounds in the ACT compatible output, default "+str(tf), type=float, default=tf)
    args = parser.parse_args()
    return args

def get_monomers(sel:str)->list:
    monomers = []
    for aas in sel:
        xyz = ( "xyz/monomers/%s.xyz" % aas )
        if os.path.exists(xyz):
            monomers.append(aas)
        else:
            print("No such file %s" % xyz)
    return monomers
    
if __name__ == "__main__":
    args = parse_args()
    gend = 0
    if args.selection and len(args.selection) > 0:
        mono1 = get_monomers(args.selection)
        mono2 = None
        if args.selection2:
            mono2 = get_monomers(args.selection2)
        
        if len(mono1) > 0:
            selection = []
            for mm in range(len(mono1)):
                if mono2:
                    for nn in range(len(mono2)):
                        selection.append("%s#%s" % ( mono1[mm], mono2[nn] ))
                        gend += 1
                else:
                    for nn in range(mm, len(mono1)):
                        selection.append("%s#%s" % ( mono1[mm], mono1[nn] ))
                        gend += 1
            with open(args.output, "w") as outf:
                for s in selection:
                    if args.act:
                        if random.uniform(0.0, 1.0) < args.test_fraction:
                            outf.write("%s|Test\n" % s)
                        else:
                            outf.write("%s|Train\n" % s)
                    else:
                        outf.write("%s\n" % s)
    if gend == 0:
        print("No monomers selected")
