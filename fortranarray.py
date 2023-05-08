'''
This module contains a function that obtains a NumPy array
and writes a Fortran code for an equivalent array.
'''

import numpy as np

def write_fortran_array(array, outfile):
    # Get the shape of the array
    print("Shape of the array:")
    print(array.shape)
