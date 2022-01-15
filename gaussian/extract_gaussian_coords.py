#!/usr/bin/env python3

"""
Rangsiman Ketkaew
July 7th, 2020

Extract the cartesian coordinates of all optimized geometries and their energies 
from Gaussian output files and save all of them as an XYZ file.

Note that the XYZ output file will be written at the directory where 
this script is executed unless optional option is used.

Normal usage:
$ python extract_gaussian_coords.py /path/to/folder/containing/Gaussian_outputs/
"""

import argparse
import glob
import os
import shutil
import sys

import numpy as np


def getVector(nline, start, end):
    vector = []
    for line in nline[start + 1:end]:
        dat = line.split()
        vector.append(dat)

    return vector

def index2symbol(atom):
    atom = int(atom)
    if   atom ==   1 : atom = "H"
    elif atom ==   2 : atom = "He"
    elif atom ==   3 : atom = "Li"
    elif atom ==   4 : atom = "Be"
    elif atom ==   5 : atom = "B"
    elif atom ==   6 : atom = "C"
    elif atom ==   7 : atom = "N"
    elif atom ==   8 : atom = "O"
    elif atom ==   9 : atom = "F"
    elif atom ==  10 : atom = "Ne"
    elif atom ==  11 : atom = "Na"
    elif atom ==  12 : atom = "Mg"
    elif atom ==  13 : atom = "Al"
    elif atom ==  14 : atom = "Si"
    elif atom ==  15 : atom = "P"
    elif atom ==  16 : atom = "S"
    elif atom ==  17 : atom = "Cl"
    elif atom ==  18 : atom = "Ar"
    elif atom ==  19 : atom = "K"
    elif atom ==  20 : atom = "Ca"
    elif atom ==  21 : atom = "Sc"
    elif atom ==  22 : atom = "Ti"
    elif atom ==  23 : atom = "V"
    elif atom ==  24 : atom = "Cr"
    elif atom ==  25 : atom = "Mn"
    elif atom ==  26 : atom = "Fe"
    elif atom ==  27 : atom = "Co"
    elif atom ==  28 : atom = "Ni"
    elif atom ==  29 : atom = "Cu"
    elif atom ==  30 : atom = "Zn"
    elif atom ==  31 : atom = "Ga"
    elif atom ==  32 : atom = "Ge"
    elif atom ==  33 : atom = "As"
    elif atom ==  34 : atom = "Se"
    elif atom ==  35 : atom = "Br"
    elif atom ==  36 : atom = "Kr"
    elif atom ==  37 : atom = "Rb"
    elif atom ==  38 : atom = "Sr"
    elif atom ==  39 : atom = "Y"
    elif atom ==  40 : atom = "Zr"
    elif atom ==  41 : atom = "Nb"
    elif atom ==  42 : atom = "Mo"
    elif atom ==  43 : atom = "Tc"
    elif atom ==  44 : atom = "Ru"
    elif atom ==  45 : atom = "Rh"
    elif atom ==  46 : atom = "Pd"
    elif atom ==  47 : atom = "Ag"
    elif atom ==  48 : atom = "Cd"
    elif atom ==  49 : atom = "In"
    elif atom ==  50 : atom = "Sn"
    elif atom ==  51 : atom = "Sb"
    elif atom ==  52 : atom = "Te"
    elif atom ==  53 : atom = "I"
    elif atom ==  54 : atom = "Xe"
    elif atom ==  55 : atom = "Cs"
    elif atom ==  56 : atom = "Ba"
    elif atom ==  57 : atom = "La"
    elif atom ==  58 : atom = "Ce"
    elif atom ==  59 : atom = "Pr"
    elif atom ==  60 : atom = "Nd"
    elif atom ==  61 : atom = "Pm"
    elif atom ==  62 : atom = "Sm"
    elif atom ==  63 : atom = "Eu"
    elif atom ==  64 : atom = "Gd"
    elif atom ==  65 : atom = "Tb"
    elif atom ==  66 : atom = "Dy"
    elif atom ==  67 : atom = "Ho"
    elif atom ==  68 : atom = "Er"
    elif atom ==  69 : atom = "Tm"
    elif atom ==  70 : atom = "Yb"
    elif atom ==  71 : atom = "Lu"
    elif atom ==  72 : atom = "Hf"
    elif atom ==  73 : atom = "Ta"
    elif atom ==  74 : atom = "W"
    elif atom ==  75 : atom = "Re"
    elif atom ==  76 : atom = "Os"
    elif atom ==  77 : atom = "Ir"
    elif atom ==  78 : atom = "Pt"
    elif atom ==  79 : atom = "Au"
    elif atom ==  80 : atom = "Hg"
    elif atom ==  81 : atom = "Tl"
    elif atom ==  82 : atom = "Pb"
    elif atom ==  83 : atom = "Bi"
    elif atom ==  84 : atom = "Po"
    elif atom ==  85 : atom = "At"
    elif atom ==  86 : atom = "Rn"
    elif atom ==  87 : atom = "Fe"
    elif atom ==  88 : atom = "Ra"
    elif atom ==  89 : atom = "Ac"
    elif atom ==  90 : atom = "Th"
    elif atom ==  91 : atom = "Pa"
    elif atom ==  92 : atom = "U"
    elif atom ==  93 : atom = "Np"
    elif atom ==  94 : atom = "Pu"
    elif atom ==  95 : atom = "Am"
    elif atom ==  96 : atom = "Cm"
    elif atom ==  97 : atom = "Bk"
    elif atom ==  98 : atom = "Cf"
    elif atom ==  99 : atom = "Es"
    elif atom == 100 : atom = "Fm"
    elif atom == 101 : atom = "Md"
    elif atom == 102 : atom = "No"
    elif atom == 103 : atom = "Lr"
    elif atom == 104 : atom = "Rf"
    elif atom == 105 : atom = "Db"
    elif atom == 106 : atom = "Sg"
    elif atom == 107 : atom = "Bh"
    elif atom == 108 : atom = "Hs"
    elif atom == 109 : atom = "Mt"
    elif atom == 110 : atom = "Ds"
    elif atom == 111 : atom = "Rg"
    elif atom == 112 : atom = "Cn"
    elif atom == 113 : atom = "Uut"
    elif atom == 114 : atom = "Fl"
    elif atom == 115 : atom = "Uup"
    elif atom == 116 : atom = "Lv"
    elif atom == 117 : atom = "Uus"
    elif atom == 118 : atom = "Uuo"
    
    return atom

def getMol(f):
    gaussian_file = open(f, "r")
    nline = gaussian_file.readlines()
    gaussian_file.close()

    all_atom = []
    all_energy = []
    all_coord = []

    # <Number of atoms>
    # <Energy>
    # <Atomic symbol> <coord_x> <coord_y> <coord_z>
    # ...

    # Extract the final single point energies
    for i in range(len(nline)):
        if "SCF Done" in nline[i]: 
            energy = nline[i].split(' ')[7].rstrip("\n\r")
            all_energy.append(energy)

    for i in range(len(nline)):
        if "Standard orientation" in nline[i]:
            start = i + 4
            for j in range(start + 2, len(nline)):
                if "---" in nline[j]:
                    end = j
                    coord = getVector(nline, start, end)
                    all_coord.append(coord)
                    # atomic symbol is in the second column
                    all_atom.append([index2symbol(coord[a][1]) for a in range(len(coord))])
                    break

    # Take only the last three column (coordinate) from the inner most of 3D array
    # and then convert all element to float

    if len(all_coord) == 0:
        return [], [], [], []
    else:
        all_coord = np.asarray(all_coord)[:, :, 3:].astype(float)

        return all_atom, all_energy, all_coord


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Extract all optimized geometries from Gaussian output file and save as an XYZ file.'
    )
    parser.add_argument(
        'gau_dir',
        metavar='Gaussian-output-directory',
        type=str,
        help='path to directory storing Gaussian output files',
    )
    parser.add_argument(
        '-o',
        '--output',
        metavar='name-of-output',
        action='store',
        help='name of XYZ file',
    )

    args = parser.parse_args()
    gaudir = args.gau_dir
    
    # Try to find *.out files
    files = glob.glob(str(gaudir) + '/' + '*.out')
    # Try to find *.log files
    if len(files) == 0:
        files = glob.glob(str(gaudir) + '/' + '*.log')

    if not files:
        sys.stderr.write("Error: Cannot find any Gaussian output file (.out or .log) in the directory you entered!\n")
        exit()

    ######### Extract atomic symbols, coordinates and energies #########

    sys.stdout.write("===> Starting extracing cartesian coordinates from Gaussian output files\n")

    A, E, C = [], [], []
    for no, file in enumerate(files):
        a, e, c = getMol(file)
        A.append(a)
        E.append(e)
        C.append(c)
        sys.stdout.write("===> File {0} : {1}\n".format(no + 1, file))

    ######### Write an XYZ file #########

    sys.stdout.write("===> Writing output files\n")
    
    name_out = args.output
    if name_out is None:
        name_out = 'extract_geom_opted'
        sys.stdout.write("===> Output name not specified, default name '%s' will be used\n"
                         % name_out
                         )
    basename = name_out.split('.')[-1]
    extxyz_out = basename + '.xyz'
    energies_out = basename + '_energies.dat'
    
    extxyz_out = gaudir + '/' + extxyz_out
    energies_out = gaudir + '/' + energies_out

    o = open(extxyz_out, 'w+')
    oe = open(energies_out, 'w+')
    # loop over files
    for i in range(len(files)):
        # loop over conformers
        for j in range(len(E[i])):
            o.write('{0}\n'.format(len(A[i][j])))
            o.write('{0:16.12f}\n'.format(float(E[i][j])))
            oe.write('{0:16.12f}\n'.format(float(E[i][j])))
            # loop over geometries
            for k in range(len(C[i][j])):
                o.write('{0}\t{1:9.8f}\t{2:9.8f}\t{3:9.8f}\n'
                        .format(A[i][j][k],
                                C[i][j][k][0],
                                C[i][j][k][1],
                                C[i][j][k][2]
                                ))
    o.close()
    oe.close()
    sys.stdout.write("===> XYZ has been saved to '%s'\n" % extxyz_out)
    sys.stdout.write("===> Energies has been saved to '%s'\n" % energies_out)
    sys.stdout.write("===> All tasks completed!\n")
