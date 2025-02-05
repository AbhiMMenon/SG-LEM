# Notes/TODO

- 0D reactor initialisation for premixed solver no longer works, this is needed for the backward-step case, but not for the Volvo case. This is due to the change to Cantera-3.0.0
- OpenFOAM test cases to be added
- Solver description



# Installation

## Requirements

1. OpenFOAM-9
    - sourced environment variables, e.g., `$FOAM_LIBBIN` etc.
    - wmake
2. gcc (not lower than 5.x)
3. GNU make (tested 4.3 and lower)
4. cmake (minimum version 3.5.0)
5. A working installation of Cantera 3.0.0 or higher

    - Tested with installation from source, see https://github.com/Cantera/cantera/releases, 
    and instructions in included file `INSTALL.md` or https://cantera.org/stable/develop/index.html

    - Or from PPA (not tested), see https://cantera.org/stable/install/ubuntu.html#sec-install-ubuntu

    - Either way, the exact file paths to `/include` and `/lib` folders in the installation are required.

6. boost library (tested version 1.83)
7. OpenMPI
6. yaml-cpp (optional)


## 0. git clone

`git clone https://github.com/AbhiMMenon/SG-LEM.git` to a directory `{XYZ}` of your choice and `cd {XYZ}/SG-LEM` 

## 1. Environment variables

1. Set root variable to current working directory

       sed -i -e 's|/home/abhi/SG-LEM|'`pwd`'|g' etc/sourceRC

2. Edit the sourceRC file for just the variables `$CANTERAINCL` and  `$CANTERALIB`. These should
    contain the exact paths to `/include` and `/lib` folders of your Cantera
    installation. See **Requirements** above

        vim etc/sourceRC
        

3. Source the file

        source etc/sourceRC

4. Optional

        echo source $(pwd)/etc/sourceRC >> ~/.bashrc
        

## 2. LEM library

Requires CMake (ver 3.5 minimum)

    cd trunk/lean-LEM-for-combustion
    cmake -B build
    cd build
    make 
    make install
    cd ..
    cd ..
    cd ..

This installs both the shared (`libLEM.so`) and static (`libLEM.a`) versions of the LEM library.

## 3. (optional) 1D splicing-stabilised flame, needs yaml-cpp

1D flamelet subject to triplet-map stirring, stabilized in the domain by
splicing operations. Note that Cantera 3.0.0 requires chemistry files in the
`.yaml` format, the previous `.cti` formats need to be converted using the
`cti2yaml` utility distributed with Cantera. This has been done for the 'sk17'
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

**NOTE**: If using openFOAM was compiled with gcc > 10.x,  `functionObjects` can show an error 
    
    error in IOstream "OSHA1stream.sinkFile_"

during runtime. See commit https://github.com/OpenFOAM/OpenFOAM-9/commit/b0c15bebd37142f3902901ed5e9a60e33ed456eb 
 for the fix if it was compiled from an pre-patch source file.



## 4. ParMGridGen 1.0

The library used to generate super-grids though mesh agglomeration. We need to compile the serial version of the library.

    cd trunk/ParMGridGen-1.0/
    make serial
    cd ..
    cd ..

## 5. dbns library

`OpenFOAM` library that uses `ParMGridGen-1.0` (serial) for coarse-graining CFD
meshes. Tested for meshes made with `blockMesh` and also for those made with
Ansys ICEM. *Should* work on any METIS format mesh.

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


## 6. (optional) userField library function

Major bug in RMS  calculated by OpenFOAM leads to negative values. Use this if
using an earlier version of OpenFOAM-9 compiled from source, i.e, before
bug-fix/patch. Library functions can be switched on during runtime.

1. Ensure `$FOAM` and `sourceRC` environment variables are sourced, as before
2. Ensure `$(FOAM_USER_LIBBIN)` is valid, as before

        cd trunk/foam9/userField
        wmake lib
        cd ..
        cd ..
        cd ..


## 7. Compile error -- "call of overloaded XX( ) is ambiguous"

If we try to compile the SG-LEM solvers, the above error is encountered. The
following lines have to be modified in Cantera header file `Func1.h`, for
Cantera 3.0.0 and earlier releases, to avoid overload ambiguity. This is due to
either FOAM or the `boost` library using similarly named functions.  May
require `sudo` privileges, i.e, `sudo bash`

    sed -i "s/return sin/return std::sin/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -i "s/return cos/return std::cos/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -i "s/return pow/return std::pow/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -i "s/return exp/return std::exp/g" $CANTERAINCL/cantera/numerics/Func1.h
    sed -i "s/return log/return std::log/g" $CANTERAINCL/cantera/numerics/Func1.h

### 7.1 Overload ambiguous function in boost

For the same reason as the above, needed for `boost` version 1_83. Root privilege is required if `boost` headers are in `/usr/`

    sudo sed -i "115s/sqrt/std::sqrt/" /usr/include/boost/math/special_functions/detail/bernoulli_details.hpp


## 8. Premixed SG-LEM solver

1. Ensure `$FOAM` and `sourceRC` environment variables are sourced, as before
2. Ensure `$(FOAM_USER_APPBIN)` is valid,

        ls  $FOAM_USER_APPBIN

    else

        mkdir -p $FOAM_USER_APPBIN

3. ***Required:*** `ParMGridGen-1.0`, `dbns`, `libLEM.so` (shared) in the above steps.
        
4. Change dir

        cd trunk/foam9/SGLEMPFoam

4. Check if 

        echo $MPI_ARCH_PATH

    or (for clusters using cent-os or similar)
        
        echo $EBROOTMPI
    give valid paths. If not, check if MPI module has been loaded, or modifiy `Make/options` line 23 or 24 using the  the correct path for MPI given by:

        mpicc --showme:compile

    Do not forget the `\` at the end of the path.

3. Compile SG-LEM premixed solver

        wmake
        cd ..
        cd ..
        cd ..

The premixed SG-LEM solver should be  installed in `$FOAM_USER_APPBIN`.

## 9. Non-premixed SG-LEM solver

 The procedure for the non-premixed solver `SGLEMFoam` is identical. It can be found in `trunk/foam9/SGLEMPFoam/SGLEMFoam`

