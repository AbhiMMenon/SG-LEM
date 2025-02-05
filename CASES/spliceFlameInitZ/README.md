# 1D Premixed flamelet initialization

This solver was used to initialize the solution matrices for the MRB case. It
uses premixed "flamelets" using the 1D LEM domain (no triplet-maps), holds the
flame in the domain using splicing long enough for conditioning in *c*-space for
each *z* bin, i.e., each premixed domain is initialised with the appropriate
air-fuel ratio. MPI is used to advance multiple flamelets for multiple *Z*-bins
at once.

1. Compile using 
        
        cmake -B build
        cd build
        make
        cd ..

2. Run using 
    
        mpirun build/spliceFlame

    At the end of the run (all *Z*-bins done), the files `scalarZC.dat` and
    `scalarZC_flag.dat` are output. The file `initZC.dat` can be visulazed using `gnuplot`

    **Note** cti files must be converted to `.yaml` , e.g., 

        python -m cantera.cti2yaml sk17/chem.cti 


