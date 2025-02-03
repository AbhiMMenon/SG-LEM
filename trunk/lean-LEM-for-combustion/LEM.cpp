#include "LEM.h"
#include "CellData.h"

#define Ith(v,i)    NV_Ith_S(v,i)         /* Ith numbers components 1..NEQ */

static int check_flag(void *returnvalue, const char *funcname, int opt);
int chemistrySource(realtype _t, N_Vector _y, N_Vector _ydot, void *user_data);
/**-----------------------------------------------------------------------------
 * @brief       constructor for LEM class. Initializes cantera IdealGasMix
 *              object usig the file string for cti file. Initializes CVODE and
 *              performs memory checks.
 *
 * @param       \input  infile  file String for cti file
 * @author      M.Oevermann
 * --------------------------------------------------------------------------*/
LEM::LEM ( std::string& infile, std::string id_)
    {

    /* default values for CH4 -A.Menon*/
    gammasCoeff   = {2.0, 0.5, -1.0, 0};
    gammas        = {0.0, 0.0, 0.0, 0.0};

    Zonly = false;
    F = 1.0;        // change for complex mechanisms for flame stability


    auto sol = Cantera::newSolution(infile, id_);
    gas = sol->thermo();
    tr = sol->transport();
    kin = sol->kinetics();

    
    ns = gas->nSpecies();
    nel=gas->nElements();
    specieList    = gas->speciesNames();
    CellData::ns  = ns;
    neq =ns+3;
    
    std::cerr << "Made soln object, neq = " << neq << std::endl;
    

        /* memory for fuelY and oxY */
    fuelY.resize(ns);
    oxY.resize(ns);

  /* 

    reltol        = 1.e-5;
    y = abstol = NULL;

    y = N_VNew_Serial(neq);
    check_flag((void *)y,  "N_VNew_Serial", 0);

    abstol    = N_VNew_Serial(neq);
    if (check_flag((void *)abstol, "N_VNew_Serial", 0)) std::cerr << "FAILED\n";

    // Set the vector absolute tolerance 
    Ith(abstol,0)   = RCONST(0);
    Ith(abstol,IT)   = RCONST(1.e-4);
    Ith(abstol,IP)   = RCONST(1.e-3);
    Ith(abstol,IRHO) = RCONST(1.e-4);

    for (int s = 0; s < ns; ++s)
        Ith(abstol, IY0+s) = 1.e-11;

    // Set dummy initial values 
    ydata = NV_DATA_S(y);
    for (int s = 0; s < neq; ++s)
        ydata[s] = 1.0;

    int flag;
    double t0 = 0.0;

 // cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    cvode_mem = CVodeCreate(CV_BDF);
    flag = CVodeInit(cvode_mem, chemistrySource, t0, y);

    check_flag(&flag,  "CVodeInit", 1);
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    check_flag(&flag,  "CVodeSVtolerances", 1);
    
  //flag = CVDense(cvode_mem, neq);
  //check_flag(&flag, (char*) "CVDense", 1);

    flag = CVodeSetMaxOrd(cvode_mem, 2);
    check_flag(&flag,  "CVodeSetMaxOrd", 1);
    flag = CVodeSetMaxNumSteps(cvode_mem, 1500);
    check_flag(&flag,  "CVodeSetMaxNumSteps", 1);
    void *f_data = (void *) this;
    flag = CVodeSetUserData(cvode_mem, f_data);
    check_flag(&flag,  "CVodeSetUserData", 1);
    */

 ///tr = Cantera::newTransportMgr("Mix", &gas);
    Ms = new double[ns];
    Rs = new double[ns];
    tmp_s1 = new double[ns];
    tmp_s2 = new double[ns];

    tmp_s3.resize(nel);
    tmp_s4.resize(nel);



    gas->getMolecularWeights(Ms);
    for (int s = 0; s < ns; ++s){
        Rs[s] = Cantera::GasConstant / Ms[s]; // J/kmol/k * kmol/kg = J/kg/K
    }


    // -- calculate formation enthalpies -A.Menon

    double* H_RT = tmp_s2;

    hFormation.resize(ns,0.0);
    gas->setState_TP(298.0, 101325.0);
    gas->getEnthalpy_RT_ref(H_RT);
    for (size_t s = 0; s < ns; s++){
        hFormation[s] = H_RT[s] * Rs[s] * 298.0; //  1 * J/kg/K * K = J/kg
      //hFormation[s] = H_RT[s] * GasConstant * 298.0;
    }


    // -- Zero D reactor
    combustor.insert(sol);
    combustor.setInitialVolume(1.0);
    sim.addReactor(combustor);

}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
//static int check_flag(void *flagvalue, char *funcname, int opt)
//{
//  int *errflag;
//
//  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
//  if (opt == 0 && flagvalue == NULL) {
//      std::cerr<< "\nSUNDIALS_ERROR: (" << funcname<< ") failed - returned NULL pointer\n\n"; 
//      return(1); 
//  }
//
//  /* Check if flag < 0 */
//  else if (opt == 1) {
//    errflag = (int *) flagvalue;
//    if (*errflag < 0) {
//        std::cerr<<  "\nSUNDIALS_ERROR: ("<< funcname<< ") failed with flag = "<< *errflag<<"\n\n";
//      return(1); }
//  }
//
//  /* Check if function returned NULL pointer - no memory allocated */
//  else if (opt == 2 && flagvalue == NULL) {
//      std::cerr<< "\nMEMORY_ERROR: ("<< funcname<< ") failed - returned NULL pointer\n\n";
//    return(1); 
//  }
//
//  return(0);
//}


static int check_flag(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
