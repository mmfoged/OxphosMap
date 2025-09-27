# OxphosMap
Mapping logFCs onto one of 5 OXPHOS complexes using PyMOL. 

# Arguments

    -p: Path to the log2FC data. This should include a column called Log2FC and GeneID (With NCBI gene IDs)
    -s: Path to the structure information data file. Must contain columns Chain, pdb and GeneID. Chain must refer to the chain name in the structure, such as “chain P”.
    -l: the log2FC limit (absolute values). Fx, a limit of 2 means that values above 2 and below -2 will just be set to 2 or -2, respectively. 
    -c: The OXPHOS complex to visualize (CI, CII, CIII, CIV or CV)

## Example: 

    pymol OxphosMap.py -- -l 0.5 -p ./Data/data_example.csv -s ./Data/All_complexes.csv  -c CI

# Required:

    pandas
    matplotlib
    argparse
    pathlib
    PyMOL (open-source)

# instructions (Mac)

Install PyMOL via homebrew

Get the path of the python interpreter used for PyMOL. Open PyMOL and type:

    import sys
    print(sys.executable)

This gives a path like /opt/homebrew/Cellar/pymol/3.1.0_1/libexec/bin/python

Then install pandas and other libraries like this (replacing the path with the path returned by PyMOL): 

    /opt/homebrew/Cellar/pymol/3.1.0_1/libexec/bin/python -m pip install pandas matplotlib pathlib agparse 