#!/usr/bin/env python
from __future__ import division, print_function
from past.builtins import xrange
import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging
# Uncomment for debugging messages
import os
SRCDIR=os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument("-d", action="store_true", help="Turn on debug mode")
parser.add_argument("gro", type=str, help="Gromacs file from which to extract sidechain dihedrals")
parser.add_argument("-ignoreres", type=str, default=["BCL", "NA", "SOL"], nargs='+', help="names of residues to ignore")
args = parser.parse_args()

if args.d:
    logging.basicConfig(level=logging.DEBUG)

GRO   = args.gro
AALIST=SRCDIR+"/aminoacid_dihedrals.txt"


def build_dihedrals(atom_list):
    all_dihedrals = []
    if len(atom_list) >= 4:
        for dih_count in xrange(len(atom_list) - 3):
            all_dihedrals.append(atom_list[dih_count:dih_count+4])
    return all_dihedrals


aa_dihedral = {}

with open(AALIST) as f:
    for l in f:
        l = l.split('#')[0]
        l = l.strip()
        if l != "":
            aa = l[:3]
            dihedrals = l[3:].split()
            aa_dihedral[aa] = dihedrals

amino_acids = aa_dihedral.keys()
for key in amino_acids:
    logging.debug(key + ": " + " ".join(aa_dihedral[key]))


# IGNORE_TYPES=["BCL", "NA", "SOL"]
IGNORE_TYPES=args.ignoreres

dihedral_list   = []
dihedral_resid  = []
dihedral_resname= []
this_resid_done = True

with open(GRO) as f:
    ## Skip first line, use second line as number of lines (total)"
    ## Note: this means that for a multi-frame file, only the first frame
    ##       will be used!
    _ = f.readline()
    n = f.readline()
    n = int(n)
    this_resid = 0
    this_resname = ""
    this_dihedral_atomtype = []
    this_dihedral_atomnum  = []
    
    ## Iterate through the number of lines specified by .gro file, line2
    for i in xrange(n):
        l = f.readline()
        label   = l[0:8]
        resid   = int(label[:5].strip())
        resname = label[-3:].strip()
        atomname= l[8:15].strip()
        atomnum = int(l[15:20].strip())

        if resname in IGNORE_TYPES:
            continue
        ## If we are in a new residue type, build local variabes
        ## Guaranteed:  this_dihedral_atomtype (str, atomtype names)
        ##              this_dihedral_atomnum (int, dummy indices)
        ##              this_resid_done = False
        if (this_resid != resid):
            this_resid_done = False
            this_resid = resid
            this_resname = resname
            logging.debug(resname.rstrip())
            badkey = (not resname in amino_acids)
            if badkey:
                raise ValueError("Amino acid {} found in GRO but not in aalist".format(resname))
            ## Look up the atomtypes involved in each dihedral
                
            this_dihedral_atomtype = build_dihedrals(aa_dihedral[resname])
            this_dihedral_atomnum  = [ [-1, -1, -1, -1] for dihedral in this_dihedral_atomtype ]
            logging.debug(this_dihedral_atomtype)
        ## For each line, iterate through all dihedrals in that residue;
        ##    if the current atomtype matches the type in the dictionary,
        ##    replace the dummy index in this_dihedral_atomnum with the
        ##    current index
        for i, dihedral in enumerate(this_dihedral_atomtype):
            if atomname in dihedral:
                atom_index = dihedral.index(atomname)
                this_dihedral_atomnum[i][atom_index] = atomnum

        ## Check that all dummy indices of "-1" have been replaced and that we have
        ##    not already written out dihedral lists
        atomnum_flat = [item for sublist in this_dihedral_atomnum for item in sublist]
        atomnum_done = [atomnum != -1 for atomnum in atomnum_flat]
        if all(atomnum_done) and not this_resid_done:
            this_resid_done = True
            for dihedral_atomnum in this_dihedral_atomnum:
                dihedral_list.append(dihedral_atomnum)
                dihedral_resid.append(this_resid)
                dihedral_resname.append(this_resname)

which_mono=1
which_resid=None
which_chi=1
for resid, dihedral, resname in zip(dihedral_resid, dihedral_list, dihedral_resname):
    if which_resid is None:
        which_resid = resid
    if resid != which_resid:
        if resid < which_resid:
            which_mono += 1
        which_resid = resid
        which_chi = 1
    dihedral = ["{:4d}".format(atomnum) for atomnum in dihedral]
    print("{}{}_M{}X{} ".format(resid, resname, which_mono, which_chi) + " ".join(dihedral))
    which_chi += 1
