'''
This module contains a function that obtains a NumPy array
and writes a Fortran code for an equivalent array.
'''

def write_fortran_array(array, outfile):
    '''
    Write a NumPy array to a Fortran code so that it can be used in the
    Fortran code for the GP3 method.
    '''
    # open the output file
    with open(outfile, "w", encoding="utf8") as f:
        for i in range(array.shape[1]):
            if i == 0:
                continue
            for j in range(array.shape[0]):
                print(f"gmunu({j+1},{i}) = {array[j,i]:.8f}", file=f)
    f.close()
