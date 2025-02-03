# Installation

Requirements:
1. OpenFOAM-9
2. gcc (not lower than 5.x)
3. yaml-cpp (optional)
4. Cantera 3.0.0
5. Cmake (minimum version 3.5.0)

## Environment variables

1. Edit the sourceRC file 

        vim etc/sourceRC
        
    for just the variables `$CANTERAINCL` and  `$CANTERALIB`. These should
    contain the exact paths to `/include` and `/lib` folders of your Cantera
    installation.

        

2. Soure the RC file in trunk

        . trunk/sourceRC

## LEM library

Requires CMake (ver 3.5 minimum)

    cd trunk/lean-LEM-for-combustion
    cmake -B build
    cd build
    make 
    make install
    cd ..
    cd ..
    cd ..

## (optional) 1D splicing-stabilised flame (needs yaml-cpp)

1D flamelet subject to triplet-map stirring, stabilized in the domain by
splicing operations. Note that Cantera 3.0.0 requires chemistry files in the
`.yaml` format, the previous `.cti` formats need to be converted using the
`ctml2yaml` utility distributed with Cantera. This has been done for the 'sk17'
mechanism for this test case, but not the others listed in the `input.yaml`
file.

    cd cases/spliceFlameSingle/
    cmake -B build
    cd build
    make
    cd ..

Test using 

    ./build/spliceFlame
    

    cd ..
    cd ..




## ParMGridGen 1.0

The library used to generate super-grids though mesh agglomeration. First we need to compile the serial version of the library.

    cd trunk/ParMGridGen-1.0/
    make serial
    cd ..
    cd ..

## dbns library

`OpenFOAM` library that uses `ParMGridGen-1.0` (serial) for coarse-graining CFD meshes. Tested for meshes made with `blockMesh` and made with Ansys ICEM.

1. Ensure `OpenFOAM` environment variables are sourced, for `wmake`
2. Ensure `$(FOAM_USER_LIBBIN)` is set

        echo  $FOAM_USER_LIBBIN

    if not 

        mkdir -p $FOAM_USER_LIBBIN
3.      cd trunk/foam9/dbns
        wmake lib
        cd ..
        cd ..

4. `libdbns` will be installed in `$(FOAM_USER_LIBBIN)`

