#!/usr/bin/env python
"""
This script runs qvSZP for all elements in the periodic table,
exports the 2-el integrals from ORCA,
sorts and averages them,
and writes them into a file in the format required by GP3.
"""

# Python script for reading in an JSON file and inserting numbers into numpy arrays
import shutil
from pathlib import Path
import sys
import argparse
import subprocess as sp
import numpy as np
from inthandler import (
    jsonhandler_resorting_legacy,
    jsonhandler_no_resorting,
    modtwoelints_analytic_average_legacy,
    average_shell_exchange_integrals,
)
from plot import plot_onexc_ints
from fortranarray import write_fortran_array, write_fortran_data
from strucio import xyzwriter
from q_cn_import import read_q_cn

QVSZP_PATH = "qvSZP"
qvszp_binary = shutil.which(QVSZP_PATH)
if qvszp_binary is None:
    raise ImportError(f"Could not find {QVSZP_PATH} in $PATH.")
print("Binary used:")
print(qvszp_binary)
ORCA_PATH = "orca"
orca_binary = shutil.which(ORCA_PATH)
if orca_binary is None:
    raise ImportError(f"Could not find {ORCA_PATH} in $PATH.")
print("Binary used:")
print(orca_binary)

PSE: dict[int, str] = {
    0: "X",
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
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}
PSE_NUMBERS: dict[str, int] = {k.lower(): v for v, k in PSE.items()}
PSE_SYMBOLS: dict[int, str] = {v: k.lower() for v, k in PSE.items()}


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
parser.add_argument(
    "--legacy",
    action="store_true",
    default=False,
    help="use the legacy version of ORCA_2JSON",
)
parser.add_argument(
    "--read-only",
    "-r",
    action="store_true",
    default=False,
    help="Do not perform calculations. Just read in the JSON files.",
)
parser.add_argument(
    "--specific_element",
    "-se",
    type=str,
    default=None,
    help="Only run for a specific element",
)
args = parser.parse_args()

if args.external_charges:
    if args.verbose:
        print("Using external charges")
    q_cn_dict = read_q_cn(Path("q_cn.dat").resolve(), args.verbose)

if args.specific_element:
    if args.specific_element not in PSE_NUMBERS:
        raise ValueError(f"Element {args.specific_element} not in the periodic table.")

onecxcints = np.zeros((9, 104))
# if onexcints.npy is a file, load it and plot it
if Path("onecxcints.npy").is_file():
    onecxcints = np.load("onecxcints.npy")
    plot_onexc_ints(onecxcints)
    # write the onecenterxcints array to Fortran code.
    write_fortran_array(onecxcints, "onecxcints_array.f90")
    write_fortran_data(onecxcints, "onecxcints_data.f90")
    if args.read_only:
        sys.exit(0)

# print current directory via pathlib
print("Current working directory:", Path.cwd())
for i in range(1, 104):
    if args.specific_element:
        if i != PSE_NUMBERS[args.specific_element]:
            continue
    print(f"Running for element {PSE_SYMBOLS[i]}")

    # check if a directory with the element name exists and if not create it
    element_path = Path(PSE_SYMBOLS[i]).resolve()
    # do the steps in the if clause only if calculating integrals from scratch is desired.
    if not args.read_only:
        element_path.mkdir(exist_ok=True)
        print(f"Successfully created the directory {element_path}")

        # Copy "hf_q-vSZP.json.conf" to current directory
        shutil.copy("hf_q-vSZP.json.conf", element_path)

        # Write the xyz file with the uppercase element symbols and the lower case file name
        strucfile = element_path / (PSE_SYMBOLS[i] + ".xyz")
        if i % 2:
            with open(element_path / ".UHF", "w", encoding="utf8") as f:
                f.write("1")
            f.close()
        xyzwriter(PSE_SYMBOLS[i].upper(), strucfile)

        # if external charges are used, write the charges to the file "ext.charges"
        if args.external_charges:
            with open(element_path / "ext.charges", "w", encoding="utf8") as f:
                f.write(
                    str(q_cn_dict[str(i)]["q"]) + " " + str(q_cn_dict[str(i)]["CN"])
                )
            f.close()
            CHARGEMODEL = "ext"
        else:
            CHARGEMODEL = "ceh_external"

        try:
            process = sp.run(
                [
                    qvszp_binary,
                    "--struc",
                    strucfile.name,
                    "--bfile",
                    "/Users/marcelmueller/source/qvSZP/q-vSZP_basis/basisq-3.0.0",
                    "--efile",
                    "/Users/marcelmueller/source/qvSZP/q-vSZP_basis/ecpq",
                    "--mpi",
                    "4",
                    "--guess",
                    "hcore",
                    "--hfref",
                    "--scf-cycles",
                    "1000",
                    "--cm",
                    CHARGEMODEL,
                    "--notrahf",
                ],
                cwd=element_path,
                capture_output=True,
                text=True,
                check=True,
            )
        except sp.CalledProcessError as err:
            print(f"Error in qvSZP execution:\n{err.stderr}")
            raise SystemExit(1) from err
        if args.verbose:
            print("Output: ", process.stdout)
        if args.dry_run:
            sys.exit(0)

        with open(element_path / "orca.out", "w", encoding="utf8") as f:
            try:
                process = sp.run(
                    [orca_binary, "hf_q-vSZP.inp"],
                    stdout=f,
                    cwd=element_path,
                    check=True,
                    text=True,
                )
            except sp.CalledProcessError as err:
                print(f"Error in ORCA execution:\n{err.stderr}")
                raise SystemExit(1) from err
        f.close()

        with open(element_path / "orca_2json.out", "w", encoding="utf8") as f:
            try:
                process = sp.run(
                    ["orca_2json", "hf_q-vSZP.gbw"],
                    stdout=f,
                    cwd=element_path,
                    check=True,
                    text=True,
                )
            except sp.CalledProcessError as err:
                print(f"Error in orca_2json execution:\n{err.stderr}")
                raise SystemExit(1) from err
        f.close()

    # Read in the json file
    if args.legacy:
        twoelints = jsonhandler_resorting_legacy(
            element_path / "hf_q-vSZP.json", PSE_SYMBOLS[i], args.verbose
        )
    else:
        twoelints = jsonhandler_no_resorting(
            element_path / "hf_q-vSZP.json", PSE_SYMBOLS[i], args.verbose
        )

    if args.verbose and args.legacy:
        # print the twoelints numpy array
        print(twoelints)

    if args.legacy:
        msindo_xc_ints = modtwoelints_analytic_average_legacy(
            twoelints, i, args.verbose
        )
    else:
        msindo_xc_ints = average_shell_exchange_integrals(twoelints, i, args.verbose)
    # incorporate the msindo xc integrals into the onecenterxcints array for the current element
    # the whole vector of size 5 is copied into the array at the position of the current element
    onecxcints[:, i] = msindo_xc_ints

if args.verbose:
    # print the onecenterxcints array
    print("Final 1c-XC ints:")
    print(onecxcints)

# dump onexcints to a file for later use
np.save("onecxcints.npy", onecxcints)
# plot the onecenterxcints array
plot_onexc_ints(onecxcints)
# write the onecenterxcints array to Fortran code.
write_fortran_array(onecxcints, "onecxcints_array.f90")
write_fortran_data(onecxcints, "onecxcints_data.f90")
