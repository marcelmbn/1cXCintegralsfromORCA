"""
Module reading element-specific individual charges and CNs from a file.
"""

from pathlib import Path


def read_q_cn(filename: Path, verb: bool) -> dict[str, dict[str, float]]:
    """
    Read the element-specific individual charges and CNs from the file.
    The file looks as follows:
    # for the elements 1-103 and all mol. in the GP3 fit set
    # average Hirshfeld charges (col 2) and q-vSZP CN (col 3)
    # CN/AC values averaged
    # Thu Apr 18 08:53:53 CEST 2024
           1  3.943834858121028E-002  0.392107130194450
           2  8.350069593223589E-002  8.105578447136906E-002
           3  0.224099080403601       0.991005665795985
           4  0.124357658093702       0.749951258771420
           5  5.444579768400955E-002   1.15437253100212
           6 -8.490959279665408E-003   1.66914161822950
           7 -0.127500860446933        1.42502747232808
           8 -0.248910794931116       0.871805073234203
           9 -0.217047059422607       0.633401380888131
          10  9.877375631009190E-002  8.767100282191988E-002
          11  0.310076137779227       0.874064634821176
          12  0.306909141188192       0.875476607806435

    Args:
        filename (str): The name of the file to read.
        verb (bool): Verbose output.

    Returns:
        dict: The dictionary with the element-specific individual charges
              and CNs: dict[element, dict["q","CN"]].

    """
    # Read the file
    try:
        with open(filename, "r", encoding="utf8") as file:
            lines = file.readlines()
        file.close()
    except Exception as e:
        print(f"Could not open file {filename}.")
        raise e

    # Initialize the dictionary
    q_cn_dict: dict[str, dict[str, float]] = {}

    # Read the lines
    for line in lines:
        if line[0] == "#":
            continue
        # Split the line
        values = line.split()
        # Get the element
        element = values[0]
        # Get the charge
        q = float(values[1])
        # Get the CN
        cn = float(values[2])
        # Add the element to the dictionary
        q_cn_dict[element] = {"q": q, "CN": cn}

    if verb:
        print(q_cn_dict)

    return q_cn_dict
