"""
This script reads the 2-el integrals from a JSON file and writes them into a numpy array.
The indices are swapped to match the order of the integrals in the GP3 method.
Relevant integrals are separated and written into a file in the format required by GP3.
"""

import json
from pathlib import Path
import numpy as np


def jsonhandler_resorting_legacy(inpfile: Path, outprefix: str, verb: bool):
    """
    Read in the JSON file and write the integrals into a numpy array and a file.
    Indices are swapped to match the order of the integrals in the GP3 method.
    """
    # Open the JSON file
    with open(inpfile, encoding="utf8") as json_file:
        data = json.load(json_file)
    json_file.close()

    # Create numpy arrays for the data
    print("Creating numpy arrays...")

    # Get the 2elIntegrals list from the data
    integrals_list = data["Molecule"]["2elIntegrals"]["AO_PQRS"][0]

    integrals_array = np.array(integrals_list)
    # Swap the indices
    # exchange each occurence of the integer "3" with "11"
    integrals_swapped = integrals_array
    # reordering for p functions
    integrals_swapped = np.where(integrals_swapped == 1, 13, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 2, 11, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 3, 12, integrals_swapped)
    # reordering for d functions
    integrals_swapped = np.where(integrals_swapped == 4, 15, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 5, 17, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 6, 18, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 7, 14, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 8, 16, integrals_swapped)
    # exchange back
    integrals_swapped = np.where(integrals_swapped == 11, 1, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 12, 2, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 13, 3, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 14, 4, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 15, 5, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 16, 6, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 17, 7, integrals_swapped)
    integrals_swapped = np.where(integrals_swapped == 18, 8, integrals_swapped)

    if verb:
        print("Swapped indices:")
        print(integrals_swapped)
        print("Original indices:")
        print(integrals_array)

    outfile_orcaorder = inpfile.parent / (outprefix + "_integrals_ORCAorder.dat")
    outfile_gp3order = inpfile.parent / (outprefix + "_integrals_GP3order.dat")

    # convert the 2-dimensional numpy array into a 4-dimensional
    # numpy array with the first four fields as indices
    twoelints = np.zeros((9, 9, 9, 9))
    for row in integrals_swapped:
        twoelints[int(row[0]), int(row[1]), int(row[2]), int(row[3])] = row[4]
    # Write both numpy arrays in a file (within the array format)
    # but write the first four columns as integers and the last one as float
    print("Writing numpy arrays to file...")
    # Introduce a header line with the indices "p q r s" and "integral"
    np.savetxt(
        outfile_orcaorder,
        integrals_array,
        fmt="%3i %3i %3i %3i %8.5f",
        header="Array of 2-el integrals.\n\
0 -> s, 1 -> pz, 2 -> px, 3 -> py\np   q   r   s  <integral value>",
    )
    np.savetxt(
        outfile_gp3order,
        integrals_swapped,
        fmt="%3i %3i %3i %3i %8.5f",
        header="Array of 2-el integrals.\n\
0 -> s, 1 -> px, 2 -> py, 3 -> pz\np   q   r   s  <integral value>",
    )
    return twoelints


def jsonhandler_no_resorting(inpfile: Path, outprefix: str, verb: bool):
    """
    Read in the JSON file and write the integrals into a numpy array and a file.
    """
    # Open the JSON file
    with open(inpfile, encoding="utf8") as json_file:
        data = json.load(json_file)
    json_file.close()

    # Create numpy arrays for the data
    print("Creating numpy arrays...")

    # Get the 2elIntegrals list from the data
    # integrals_list: list of 2-el integral data (rows with 5 elements, 4 indices and 1 value)
    integrals_list = data["Molecule"]["2elIntegrals"]["AO_PQRS"][0]
    # integrals_array: two-dimensional numpy array of 2-el integrals data
    integrals_array = np.array(integrals_list)
    # outfile_orcaorder: output file name for the 2-el integrals in ORCA order
    outfile_orcaorder = inpfile.parent / (outprefix + "_integrals_ORCAorder.dat")
    # twoelints: four-dimensional numpy array of 2-el integrals data
    #            all indices are in the range 0-15 for the 16 basis functions
    #            s: 0
    #            pz: 1, px: 2, py: 3
    #            dz2: 4, dxz: 5, dyz: 6, dx2-y2: 7, dxy: 8
    #            fz3: 9, fxz2: 10, fyz2: 11, fzx2-y2: 12, fxyz: 13, fx(x2-3y2): 14, fy(3x2-y2): 15
    np.savetxt(
        outfile_orcaorder,
        integrals_array,
        fmt="%3i %3i %3i %3i %8.5f",
        header="Array of 2-el integrals.\n\
0 -> s, 1 -> pz, 2 -> px, 3 -> py\n4 -> dz2, 5 -> dxz, 6 -> dyz, 7 -> dx2-y2, 8 -> dxy\n\
9 -> fz3, 10 -> fxz2, 11 -> fyz2, 12 -> fzx2-y2, 13 -> fxyz, 14 -> fx(x2-3y2), 15 -> fy(3x2-y2)\n\
p   q   r   s  <integral value>",
    )
    twoelints = np.zeros((29, 29, 29, 29))
    for row in integrals_array:
        if verb:
            print(
                f"<p>: {int(row[0])}, <q>: {int(row[1])}, <r>: {int(row[2])}, <s>: {int(row[3])}, integral: {row[4]:.6f}"
            )
        twoelints[int(row[0]), int(row[1]), int(row[2]), int(row[3])] = row[4]
    return twoelints


def modtwoelints_analytic_average_legacy(twoelints: np.ndarray, ati: int, verb: bool):
    """
    Modify the two-electron integrals to match the MSINDO-XC method.
    """
    print("Modifying two-electron integrals...")
    msindo_xc_ints = np.zeros((5))

    # 1c-XC integral between s and p functions
    msindo_xc_ints[0] = twoelints[1, 0, 1, 0]
    if msindo_xc_ints[0] <= 1e-7:
        msindo_xc_ints[0] = twoelints[0, 1, 0, 1]
    print(f"s <-> p : {msindo_xc_ints[0]:.6f}")
    # 1c-XC integral between p and p' functions
    msindo_xc_ints[1] = twoelints[2, 1, 2, 1]
    if msindo_xc_ints[1] <= 1e-7:
        msindo_xc_ints[1] = twoelints[1, 2, 1, 2]
    print(f"p <-> p': {msindo_xc_ints[1]:.6f}")

    if ati > 2:
        # 1c-XC integral between s and d functions
        msindo_xc_ints[2] = twoelints[6, 0, 6, 0]
        if msindo_xc_ints[2] <= 1e-7:
            msindo_xc_ints[2] = twoelints[0, 6, 0, 6]
        print(f"s <-> d : {msindo_xc_ints[2]:.6f}")
        # 1c-XC integral between p and d functions
        pydxy = twoelints[2, 6, 2, 6]
        if pydxy <= 1e-7:
            pydxy = twoelints[6, 2, 6, 2]
        pydxz = twoelints[2, 7, 2, 7]
        if pydxz <= 1e-7:
            pydxz = twoelints[7, 2, 7, 2]
        pydz2 = twoelints[2, 5, 2, 5]
        if pydz2 <= 1e-7:
            pydz2 = twoelints[5, 2, 5, 2]
        pzdz2 = twoelints[3, 5, 3, 5]
        if pzdz2 <= 1e-7:
            pzdz2 = twoelints[5, 3, 5, 3]
        msindo_xc_ints[3] = (
            pydxy * 8.0 + pydxz * 4.0 + pydz2 * 2.0 + pzdz2 * 1.0
        ) / 15.0
        if verb:
            print("pydxy:    ", pydxy)
            print("pydxz:    ", pydxz)
            print("pydz2:    ", pydz2)
            print("pzdz2:    ", pzdz2)
        print(f"p <-> d : {msindo_xc_ints[3]:.6f}")

        # 1c-XC integral between d and d' functions
        dxydx2y2 = twoelints[6, 4, 6, 4]
        if dxydx2y2 <= 1e-7:
            dxydx2y2 = twoelints[4, 6, 4, 6]
        dxydyz = twoelints[6, 8, 6, 8]
        if dxydyz <= 1e-7:
            dxydyz = twoelints[8, 6, 8, 6]
        dxydz2 = twoelints[6, 5, 6, 5]
        if dxydz2 <= 1e-7:
            dxydz2 = twoelints[5, 6, 5, 6]
        dxzdz2 = twoelints[7, 5, 7, 5]
        if dxzdz2 <= 1e-7:
            dxzdz2 = twoelints[5, 7, 5, 7]
        msindo_xc_ints[4] = (
            dxydx2y2 * 1.0 + dxydyz * 5.0 + dxydz2 * 2.0 + dxzdz2 * 2.0
        ) / 10.0
        if verb:
            print("dxydx2y2: ", dxydx2y2)
            print("dxydyz:   ", dxydyz)
            print("dxydz2:   ", dxydz2)
            print("dxzdz2:   ", dxzdz2)
        print(f"d <-> d': {msindo_xc_ints[4]:.6f}")
    return msindo_xc_ints


def average_shell_exchange_integrals(
    ints: np.ndarray, ati: int, verb: bool
) -> np.ndarray:
    """
    Modify the two-electron integrals to match the MSINDO-XC method.
    NOTE: # Numbering of spherical AOs:

    General case:
    S -> 1
    P -> 4
    D -> 9
    F -> 16

    Caution: Numbering seems weird as f basis for Ln's/Ac's is not large-core ECP:
    S -> 1
    S -> 2
    S -> 3
    P -> 6
    P -> 9
    P -> 12
    D -> 17
    D -> 22
    F -> 29

    Caution: Neither for Fr, Ra:
    S -> 1
    S -> 2
    P -> 5
    P -> 8
    D -> 13
    F -> 20
    """
    s_valence_range = [0]
    if 57 < ati < 72 or ati > 88:
        s_valence_range = [2]
    elif 86 < ati < 89:
        s_valence_range = [1]
    p_valence_range = range(1, 4)
    if 57 < ati < 72 or ati > 88:
        p_valence_range = range(9, 12)
    elif 86 < ati < 89:
        p_valence_range = range(5, 8)
    d_valence_range = range(4, 9)
    if 57 < ati < 72 or ati > 88:
        d_valence_range = range(17, 22)
    elif 86 < ati < 89:
        d_valence_range = range(8, 13)
    # NOTE: Only conducted for Ln's/Ac's!
    f_valence_range = range(22, 29)
    print("Modifying two-electron integrals...")
    msindo_xc_ints = np.zeros((9))
    #################
    # s-p integrals:
    #################
    print("s-p integrals:")
    counter: int = 0
    sum_s_p: float = 0.0
    effective_counter: int = 0
    zero_integrals: int = 0
    for i in s_valence_range:
        for j in p_valence_range:
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_s_p += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_s_p = sum_s_p / counter
    print(f"Average s-p integral: {average_s_p:.6f}")
    print(f"# contributing s-p integrals: {effective_counter}")
    msindo_xc_ints[0] = average_s_p
    #################
    # p-p' integrals:
    #################
    print("p-p' integrals:")
    counter = 0
    sum_p_p: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in p_valence_range[:-1]:
        for j in range(i + 1, p_valence_range[-1]):
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_p_p += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_p_p = sum_p_p / counter
    print(f"Average p-p' integral: {average_p_p:.6f}")
    print(f"# contributing p-p' integrals: {effective_counter}")
    msindo_xc_ints[1] = average_p_p
    #################
    # s-d integrals:
    #################
    print("s-d integrals:")
    counter = 0
    sum_s_d: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in s_valence_range:
        for j in d_valence_range:
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_s_d += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_s_d = sum_s_d / counter
    print(f"Average s-d integral: {average_s_d:.6f}")
    print(f"# contributing s-d integrals: {effective_counter}")
    msindo_xc_ints[2] = average_s_d
    #################
    # p-d integrals:
    #################
    print("p-d integrals:")
    counter = 0
    sum_p_d: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in p_valence_range:
        for j in d_valence_range:
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_p_d += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_p_d = sum_p_d / counter
    print(f"Average p-d integral: {average_p_d:.6f}")
    print(f"# contributing p-d integrals: {effective_counter}")
    msindo_xc_ints[3] = average_p_d
    #################
    # d-d' integrals:
    #################
    print("d-d' integrals:")
    counter = 0
    sum_d_d: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in d_valence_range[:-1]:
        for j in range(i + 1, d_valence_range[-1]):
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_d_d += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_d_d = sum_d_d / counter
    print(f"Average d-d' integral: {average_d_d:.6f}")
    print(f"# contributing d-d' integrals: {effective_counter}")
    msindo_xc_ints[4] = average_d_d
    if not (57 < ati < 72 or ati > 88):
        return msindo_xc_ints
    #################
    # s-f integrals:
    #################
    print("s-f integrals:")
    counter = 0
    sum_s_f: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in s_valence_range:
        for j in f_valence_range:
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_s_f += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_s_f = sum_s_f / counter
    print(f"Average s-f integral: {average_s_f:.6f}")
    print(f"# contributing s-f integrals: {effective_counter}")
    msindo_xc_ints[5] = average_s_f
    #################
    # p-f integrals:
    #################
    print("p-f integrals:")
    counter = 0
    sum_p_f: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in p_valence_range:
        for j in f_valence_range:
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_p_f += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_p_f = sum_p_f / counter
    print(f"Average p-f integral: {average_p_f:.6f}")
    print(f"# contributing p-f integrals: {effective_counter}")
    msindo_xc_ints[6] = average_p_f
    #################
    # d-f integrals:
    #################
    print("d-f integrals:")
    counter = 0
    sum_d_f: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in d_valence_range:
        for j in f_valence_range:
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_d_f += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_d_f = sum_d_f / counter
    print(f"Average d-f integral: {average_d_f:.6f}")
    print(f"# contributing d-f integrals: {effective_counter}")
    msindo_xc_ints[7] = average_d_f
    #################
    # f-f integrals:
    #################
    print("f-f integrals:")
    counter = 0
    sum_f_f: float = 0.0
    effective_counter = 0
    zero_integrals = 0
    for i in f_valence_range[:-1]:
        for j in range(i + 1, f_valence_range[-1]):
            if verb:
                print(
                    f"<p>: {j}, <q>: {i}, <r>: {j}, <s>: {i}, integral: {ints[j, i, j, i]:.6f}"
                )
            counter += 1
            sum_f_f += ints[j, i, j, i]
            if ints[j, i, j, i] > 1e-7:
                effective_counter += 1
            else:
                zero_integrals += 1
    average_f_f = sum_f_f / counter
    print(f"Average f-f integral: {average_f_f:.6f}")
    print(f"# contributing f-f integrals: {effective_counter}")
    msindo_xc_ints[8] = average_f_f

    #################
    # RETURN AVERAGE INTEGRALS PER ATOM
    #################
    return msindo_xc_ints
