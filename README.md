# Installation

Requirements:
1. OpenFOAM-9
2. gcc (not lower than 5.x)
3. yaml-cpp (optional)
4. Cantera 3.0.0
5. Cmake (minimum version 3.5.0)

## Environment variables

1. Set root variable to current working directory

       sed -i -e 's|/home/abhi/SG-LEM|'`pwd`'|g' etc/sourceRC

2. Edit the sourceRC file for just the variables `$CANTERAINCL` and  `$CANTERALIB`. These should
    contain the exact paths to `/include` and `/lib` folders of your Cantera
    installation.

        vim etc/sourceRC
        

3. Soure the RC file in trunk

        source etc/sourceRC

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

1. Compile solver       

        cd cases/spliceFlameSingle/ 
        cmake -B build 
        cd build 
        make 
        cd ..


2. Test using 

        mkdir  results 
        ./build/spliceFlame
    
3. Results can be visualized using gnuplot

         gnuplot plot.gnu
         evince flame.pdf

Exit directory 

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
2. Ensure `$(FOAM_USER_LIBBIN)` points to a folder

        ls  $FOAM_USER_LIBBIN

    if not 

        mkdir -p $FOAM_USER_LIBBIN
3.      cd trunk/foam9/dbns
        wmake lib
        cd ..
        cd ..
        cd ..

4. `libdbns` will be installed in `$(FOAM_USER_LIBBIN)`


## (optional) userField library function

Minor bug in RMS  calculated by OpenFOAM leads to negative values. Use this if using an earlier version of OpenFOAM-9 compiled from source before bug fix.

1. Ensure `OpenFOAM` environment variables are sourced, as before
2. Ensure `$(FOAM_USER_LIBBIN)` is valid, as before

        cd trunk/foam9/userField
        wmake lib
        cd ..
        cd ..
        cd ..

## Minor Cantera bugs

The following lines have to be corrected in Cantera header files distributed with Cantera 3.0.0 (and earlier releases) to avoid overload ambiguity, likely due to FOAM or boost using similarly named functions.



