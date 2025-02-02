#!/usr/bin/env python
"""
This script runs qvSZP for all elements in the periodic table,
exports the 2-el integrals from ORCA,
sorts and averages them,
and writes them into a file in the format required by GP3.
"""

# Python script for reading in an JSON file and inserting numbers into numpy arrays
import os
import sys
import argparse
import subprocess as sp
import numpy as np
from inthandler import jsonhandler, modtwoelints
from fortranarray import write_fortran_array, write_fortran_data
from strucIO import xyzwriter
from q_cn_import import read_q_cn

QVSZP_PATH = "qvSZP"
print("Binary used:")
sp.run(["which", "qvSZP"], check=True)

pesdict = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
}

# Parse command line arguments
parser = argparse.ArgumentParser(description="Run qvSZP for all elements in pesdict.")
parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    default=False,
    help="increase output verbosity",
)
parser.add_argument(
    "-ext", "--external_charges", action="store_true", help="use external charges"
)
parser.add_argument(
    "-dry", "--dry_run", action="store_true", help="only perform first input generation"
)
args = parser.parse_args()

if args.external_charges:
    if args.verbose:
        print("Using external charges")
    q_cn_dict = read_q_cn("q_cn.dat", args.verbose)

onecxcints = np.zeros((5, 87))

for i in range(1, 87):
    if 57 < i < 72:
        print(f"Skipping element {pesdict[i]}")
        continue
    print(f"Running for element {pesdict[i]}")

    # print current directory
    print("Current working directory:", os.getcwd())
    # check if a directory with the element name exists and if not create it
    if not os.path.exists(pesdict[i].lower()):
        try:
            os.mkdir(pesdict[i].lower())
        except OSError:
            print(f"Creation of the directory {pesdict[i].lower()} failed")
        else:
            print(f"Successfully created the directory {pesdict[i].lower()}")
    else:
        print(f"Directory {pesdict[i].lower()} already exists.")

    # change into the directory with the element name
    try:
        os.chdir(pesdict[i].lower())
    except OSError as err:
        print(f"Could not change directory to {pesdict[i].lower()}")
        raise SystemExit(1) from err
    print(f"Successfully changed directory to {pesdict[i].lower()}")

    # Copy "hf_q-vSZP.json.conf" to current directory
    try:
        process = sp.run(
            ["cp", "../hf_q-vSZP.json.conf", "."],
            capture_output=True,
            text=True,
            check=True,
        )
    except sp.CalledProcessError as err:
        print("Error: ", err.stderr)
        raise SystemExit(1) from err

    # Write the xyz file with the uppercase element symbols and the lower case file name
    strucfile = pesdict[i].lower() + ".xyz"
    if i % 2:
        with open(".UHF", "w", encoding="utf8") as f:
            f.write("1")
        f.close()
    xyzwriter(pesdict[i].upper(), strucfile)

    # if external charges are used, write the charges to the file "ext.charges"
    if args.external_charges:
        with open("ext.charges", "w", encoding="utf8") as f:
            f.write(str(q_cn_dict[str(i)]["q"]) + " " + str(q_cn_dict[str(i)]["CN"]))
        f.close()
        CHARGEMODEL = "ext"
    else:
        CHARGEMODEL = "ceh_external"

    try:
        process = sp.run(
            [
                QVSZP_PATH,
                "--struc",
                strucfile,
                "--bfile",
                "/Users/marcelmueller/source/qvSZP/q-vSZP_basis/basisq-3.0.0",
                "--efile",
                "/Users/marcelmueller/source/qvSZP/q-vSZP_basis/ecpq",
                "--mpi",
                "4",
                "--guess",
                "hcore",
                "--hfref",
                "--cm",
                CHARGEMODEL,
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    except sp.CalledProcessError as err:
        print("Error: ", err.stderr)
        raise SystemExit(1) from err
    else:
        if args.verbose:
            print("Output: ", process.stdout)
    if args.dry_run:
        sys.exit(0)

    with open("orca.out", "w", encoding="utf8") as f:
        try:
            process = sp.run(
                ["orca", "hf_q-vSZP.inp"],
                stdout=f,
                check=True,
                text=True,
            )
        except sp.CalledProcessError as err:
            print("Error: ", err.stderr)
            raise SystemExit(1) from err
    f.close()

    with open("orca_2json.out", "w", encoding="utf8") as f:
        try:
            process = sp.run(
                ["orca_2json", "hf_q-vSZP.gbw"],
                stdout=f,
                check=True,
                text=True,
            )
        except sp.CalledProcessError as err:
            print("Error: ", err.stderr)
            raise SystemExit(1) from err
    f.close()

    # Read in the json file
    twoelints = jsonhandler("hf_q-vSZP.json", pesdict[i].lower(), args.verbose)

    if args.verbose:
        # print the twoelints numpy array
        print(twoelints)

    msindo_xc_ints = modtwoelints(twoelints, i, args.verbose)
    # incorporate the msindo xc integrals into the onecenterxcints array for the current element
    # the whole vector of size 5 is copied into the array at the position of the current element
    onecxcints[:, i] = msindo_xc_ints

    # change back into the parent directory
    try:
        os.chdir("..")
    except OSError as err:
        print(f"Could not change directory to {os.getcwd()}.")
        raise SystemExit(1) from err

if args.verbose:
    # print the onecenterxcints array
    print(onecxcints)

# write the onecenterxcints array to Fortran code.

write_fortran_array(onecxcints, "onecxcints_array.f90")
write_fortran_data(onecxcints, "onecxcints_data.f90")
