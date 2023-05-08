#!/usr/bin/env python

# Python script for reading in an JSON file and inserting numbers into numpy arrays

import json
import numpy as np

# Open the JSON file
with open('hf.json') as json_file:
    data = json.load(json_file)

# Create numpy arrays for the data
print("Creating numpy arrays...")
print(type(data))

# Get the 2elIntegrals list from the data
integrals_list = data['Molecule']['2elIntegrals']

integrals_array = np.array(integrals_list)
# Swap the indices
# exchange each occurence of the integer "3" with "11"
integrals_swapped = integrals_array
integrals_swapped = np.where(integrals_swapped == 1, 13, integrals_swapped)
integrals_swapped = np.where(integrals_swapped == 2, 11, integrals_swapped)
integrals_swapped = np.where(integrals_swapped == 3, 12, integrals_swapped)
# exchange back
integrals_swapped = np.where(integrals_swapped == 11, 1, integrals_swapped)
integrals_swapped = np.where(integrals_swapped == 12, 2, integrals_swapped)
integrals_swapped = np.where(integrals_swapped == 13, 3, integrals_swapped)

print("Swapped indices:")
print(integrals_swapped)

print("Original indices:")
print(integrals_array)

# Write both numpy arrays in a file (within the array format) but write the first four columns as integers and the last one as float
print("Writing numpy arrays to file...")
# Introduce a header line with the indices "p q r s" and "integral"
np.savetxt("integrals_ORCAorder.txt", integrals_array, fmt='%3i %3i %3i %3i %8.5f', header="Array of 2-el integrals.\n0 -> s, 1 -> pz, 2 -> px, 3 -> py\np   q   r   s  <integral value>")
np.savetxt("integrals_GP3order.txt", integrals_swapped, fmt='%3i %3i %3i %3i %8.5f', header="Array of 2-el integrals.\n0 -> s, 1 -> px, 2 -> py, 3 -> pz\np   q   r   s  <integral value>")
