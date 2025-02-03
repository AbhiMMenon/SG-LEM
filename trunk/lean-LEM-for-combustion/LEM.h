#ifndef LEM_H
#define LEM_H

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "cantera/core.h"
//#include "cantera/transport.h"
#include "cantera/zerodim.h"

/* Header files with a description of contents used */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
//#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <boost/multi_array.hpp>



// scalar positions for ODE solver
#define IT   0
#define IP   1
#define IRHO 2
#define IY0  3

typedef boost::multi_array<double, 2> matrix2D;
typedef boost::multi_array<double, 3> matrix3D;
typedef boost::multi_array<bool, 2> matrix2Dbool;

//using namespace std;
//using namespace Cantera;

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
class LEM
{
  public:

    LEM(std::string& infile,
        std::string id_="");
    ~LEM()
    {
    //N_VDestroy(y);
    //N_VDestroy(abstol);
    //CVodeFree(&cvode_mem);
      delete [] Ms;
      delete [] Rs;
      delete [] tmp_s1;
      delete [] tmp_s2; //delete tr;
    //  free(tr);
    }

    friend int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

    std::shared_ptr<Cantera::ThermoPhase> gas;
    std::shared_ptr<Cantera::Transport> tr;
    std::shared_ptr<Cantera::Kinetics> kin;

 // Cantera::Transport*   tr;
    
    int ns, neq;             // number of species, equations
    int nel;                 // number of elements        -A.Menon

    realtype reltol;
    N_Vector y, abstol;
//  void     *cvode_mem;


    double *Ms;      // species molecular weight
    double *Rs;      // species gas constant
    double *tmp_s1;  // temporary arrays with size ns
    double *tmp_s2;

    std::vector<double> tmp_s3;  // temp arrays of size nel
    std::vector<double> tmp_s4;

    //-----------------------------------------------------------------------------------
    // --- Fuel and oxidiser stream porperties.                     -A.Menon

    std::string fuelStream;
    std::string   oxStream;

    // For CFD, default is all species in chemistry file
    std::vector<std::string> specieList;

    double fuelStreamT;
    double oxStreamT;
    double fuelRho;
    double oxRho;
    double pressure;
    double F;       // thickening factor

    std::vector<double> fuelY;
    std::vector<double> oxY;
    std::vector<int> C_select;        // list of species selected for calculating C


    // coefficients for elemental C,H,O and N in Bilger's formula.
    // eg. for  Methane in air gammaCoeff = {2.0, 0.5, -1, 0} *default*
    std::vector<double> gammasCoeff;
    std::vector<double> gammas;
    double beta1, beta2;       // Coeff's for Bilgers' formula


    std::vector<double> hFormation;

    //-----------------------------------------------------------------------------------
    //  Zero D reactor April 8th 2021; trial for enthalpy source term

//  Cantera::Reactor combustor;
    Cantera::IdealGasConstPressureReactor combustor;
    Cantera::ReactorNet sim;


    // -- mem for non-premixed data
    std::vector<double> TEquil_z;
    std::vector<double> TUnbZ;
    matrix2D  YEquilZ;
    matrix2D  YUnbZ;
    matrix3D  BetaShapes;

    bool Zonly;
};

#endif
