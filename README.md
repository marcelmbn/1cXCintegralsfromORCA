# 1cXCintegralsfromORCA

Script for ...
- calculating the 1-center exchange integrals using ORCA with HF and the q-vSZP basis set
- averaging the integrals shell-wise
- plot shell-wise averaged integrals (s-p, p-p', ...)
- write the integrals as Fortran data statement to be used in development versions of the g-xTB code


## Usage


Copy the following files to your working directory (should be empty before):
- `hf_q-vSZP.json.conf` (`orca_2json` configuration file)
- `q_cn.dat` (average charges and coordination numbers for each element)

In the working directory, execute:

```bash
python ~/source/ORCA_2JSON_reader/main.py --ext
```