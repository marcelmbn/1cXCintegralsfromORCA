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

def write_fortran_data(array, outfile: str) -> None:
    '''
    Write a NumPy array to a Fortran data statement so that it can be used in the
    Fortran code for the GP3 method.
    '''
    with open(outfile, "w", encoding="utf8") as f:
        print(f"real(wp), save, dimension(5, {array.shape[1]-1}) :: gmunu", file=f)
        for i in range(array.shape[1]):
            if i == 0:
                continue
            print(f"data gmunu(:,{i}) / &", file=f)
            for j in range(0,9):
                if j == 0:
                    print(f"& {array[j,i]:.7f}_wp, ", end="", file=f)
                elif j == 8:
                    print(f"0.0000000_wp /", file=f)
                elif j > 4:
                    print(f"0.0000000_wp, ", end="", file=f)
                else:
                    print(f"{array[j,i]:.7f}_wp, ", end="", file=f)
