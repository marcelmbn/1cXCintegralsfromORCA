'''
This script reads the 2-el integrals from a JSON file and writes them into a numpy array.
The indices are swapped to match the order of the integrals in the GP3 method.
Relevant integrals are separated and written into a file in the format required by GP3.
'''

import json
import numpy as np

def jsonhandler(inpfile, outprefix, verb):
    '''
    Read in the JSON file and write the integrals into a numpy array and a file.
    Indices are swapped to match the order of the integrals in the GP3 method.
    '''
    # Open the JSON file
    with open(inpfile, encoding="utf8") as json_file:
        data = json.load(json_file)
    json_file.close()

    # Create numpy arrays for the data
    print("Creating numpy arrays...")

    # Get the 2elIntegrals list from the data
    integrals_list = data["Molecule"]["2elIntegrals"]

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

    outfile_orcaorder = outprefix + "_integrals_ORCAorder.dat"
    outfile_gp3order = outprefix + "_integrals_GP3order.dat"

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

def modtwoelints(twoelints,ati,verb):
    '''
    Modify the two-electron integrals to match the MSINDO-XC method.
    '''
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
    # 1c-XC integral between s and d functions
    msindo_xc_ints[2] = twoelints[6, 0, 6, 0]
    if msindo_xc_ints[2] <= 1e-7:
        msindo_xc_ints[2] = twoelints[0, 6, 0, 6]
    print(f"p <-> d : {msindo_xc_ints[2]:.6f}")

    if ati > 2:
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
        msindo_xc_ints[3] = ( pydxy * 8.0 + pydxz * 4.0 + pydz2 * 2.0 + pzdz2 * 1.0 ) / 15.0
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
        dxydyz   = twoelints[6, 8, 6, 8]
        if dxydyz <= 1e-7:
            dxydyz = twoelints[8, 6, 8, 6]
        dxydz2   = twoelints[6, 5, 6, 5]
        if dxydz2 <= 1e-7:
            dxydz2 = twoelints[5, 6, 5, 6]
        dxzdz2   = twoelints[7, 5, 7, 5]
        if dxzdz2 <= 1e-7:
            dxzdz2 = twoelints[5, 7, 5, 7]
        msindo_xc_ints[4] = (dxydx2y2 * 1.0 + dxydyz * 5.0 + dxydz2 * 2.0 + dxzdz2 * 2.0) / 10.0
        if verb:
            print("dxydx2y2: ", dxydx2y2)
            print("dxydyz:   ", dxydyz)
            print("dxydz2:   ", dxydz2)
            print("dxzdz2:   ", dxzdz2)
        print(f"d <-> d': {msindo_xc_ints[4]:.6f}")
    return msindo_xc_ints
