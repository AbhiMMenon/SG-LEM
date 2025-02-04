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
        

3. Source the file

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

`OpenFOAM` library that uses `ParMGridGen-1.0` (serial) for coarse-graining CFD meshes. Tested for meshes made with `blockMesh` and also for those made with Ansys ICEM. *Should* work on any METIS format mesh.

1. Ensure `$FOAM` environment variables are sourced, for `wmake`, as well as `sourceRC`
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

Major bug in RMS  calculated by OpenFOAM leads to negative values. Use this if using an earlier version of OpenFOAM-9 compiled from source, i.e, before bug-fix/patch.

1. Ensure `$FOAM` and `sourceRC` environment variables are sourced, as before
2. Ensure `$(FOAM_USER_LIBBIN)` is valid, as before

        cd trunk/foam9/userField
        wmake lib
        cd ..
        cd ..
        cd ..

## Compile error -- "call of overloaded XX( ) is ambiguious"

The following lines have to be modified in Cantera header file `Func1.h`, for Cantera 3.0.0 and earlier releases, to avoid overload ambiguity, likely due to FOAM or boost using similarly named functions.
Might need `sudo` privilege, i.e, `sudo bash`

    sed -ie "s/return sin/return std::sin/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -ie "s/return cos/return std::cos/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -ie "s/return pow/return std::pow/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -ie "s/return exp/return std::exp/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -ie "s/return log/return std::log/g" $CANTERAINCL/cantera/numerics/Func1.h

### Overload ambiguious function in boost

For the same reason as the above, needed for `boost` version 1_83. Root privilege is required if `boost` headers are in `/usr/`

    sudo sed -ie "115s/sqrt/std::sqrt/" /usr/include/boost/math/special_functions/detail/bernoulli_details.hpp
    
    

    sed -ie "s/return sin/return std::sin/g" $CANTERAINCL/cantera/numerics/Func1.h

