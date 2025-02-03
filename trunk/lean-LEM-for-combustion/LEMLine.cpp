#include "LEM.h"
#include "LEMLine.h"
#include "CellData.h"
//#include <cantera/equilibrium.h>

/**
 * @brief       Constructor for class LEMLINE. No physically meaningfull data,
 * needs initializeDataOnLine() or some variant.
 * @param       _nlem          \input          number of LEM wafers
 * @param       _length        \input          length (m)
 * @param       _ut            \input          integral length scale (m)
 * @param       Re_t           \input          Turbulent Reynolds Number
 * @param       clambda        \input          Model constant
 * @param       Neta           \input          Model constant
 * @param       Neta           \input          Model constant
 * @param       tfile          \input          File pointer for output file
 *
 * @author      M. Oevermann
 */
LEMLINE::LEMLINE(
        const int       _nlem,
        const double    _length,
        const double    _lt,

        std::ofstream& tfile,
        LEM     &_lem,

        int     _ZBinSize,
        int     _CBinSize
    )
    :
    nlem        (_nlem),
    length      (_length),
    lt          (_lt),
    lem         (_lem),
    tfile_lem   (tfile),

    ZBinSize    (_ZBinSize),
    CBinSize    (_CBinSize),

    wallPresent (false),
    time        (0.0)
{

    const double dx = length / nlem;
    const int    ns = lem.ns;

    CellData c(ns, dx);

    for (int i = 0; i < nlem; ++i){
        cells.push_back(c);
    }


    // -- LEM line stuff
    counter     = 0;
    dtLEM       = 1E-7;
    Energy      = 0.0;
    noofEddies  = 5;
    timeNextTripletMap   = 1e5;
    Total_EddiesLength   = 0.1;
    Average_EddiesLength = 0.05;

    // -- Conditioning and stream data
    nSpecies    = lem.ns;
    Cstep       = 1./(CBinSize-1);
    Zstep       = 1./(ZBinSize-1);
    bFac        = 1.0;
    ZMax        = 1.0;
    N_ETA       =NETA_DEF;
    C_LAMBDA    =CLAMBDA_DEF;
    //setStreams();
}



/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::initializeDataOnLine()
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;

   Energy = 0.0;

   double T0=700;
   gas->setState_TPX(T0, Cantera::OneAtm, "C7H16:1.9090909090, O2:21.0, N2:79.0");
   double* Y_ = lem.tmp_s1;


   gas->getMassFractions(Y_);
   for (it = cells.begin(); it != cells.end(); ++it)
   {

     (*it).p   = gas->pressure();
     (*it).T   = gas->temperature();
     (*it).rho = gas->density();
     gas->getMassFractions((*it).Y);

     (*it).m = gas->density() * (*it).dx;
     double E_cell = (gas->enthalpy_mass() - (*it).p / (*it).rho) * (*it).m;
     (*it).E = E_cell;
     Energy += E_cell;
     T0 += 2;		//Why do this?
     gas->setState_TPY(T0, Cantera::OneAtm, Y_);

   }
}

/*-----------------------------------------------------------------------------
 * O3/air on bottom and NO/air on top at 298 K, 1 atm.
 * --------------------------------------------------------------------------*/
void LEMLINE::initializeDataOnLineForMixingLayerTest()
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;

   long unsigned int nlem_div2 = cells.size()/2;
   long unsigned int counter=0;

   double *Y = lem.tmp_s1;

   for (it = cells.begin(); it != cells.end(); ++it)
   {

     if (counter<nlem_div2)
     {
       gas->setState_TPY(298, Cantera::OneAtm, "O3:4E-6, O2:0.23, N2:0.769996");
       gas->getMassFractions(Y);
       gas->setState_TPY(298, Cantera::OneAtm, Y); // redundant??
     }
     else
     {
       gas->setState_TPY(298, Cantera::OneAtm, "NO:4E-6, O2:0.23, N2:0.769996");
       gas->getMassFractions(Y);
       gas->setState_TPY(298, Cantera::OneAtm, Y);
     }

     (*it).p   = gas->pressure();
     (*it).T   = gas->temperature();
     (*it).rho = gas->density();
     gas->getMassFractions((*it).Y);
     (*it).m = gas->density() * (*it).dx;
     (*it).E = (gas->enthalpy_mass() - (*it).p / (*it).rho) * (*it).m;

     counter++;

   }

}

/*-----------------------------------------------------------------------------
 * Hot products on bottom and 300 K reactants on top for C3H8/air phi=0.65 1 atm.
 * --------------------------------------------------------------------------*/
void LEMLINE::initializeDataOnLineForPremixedFlameTest(long unsigned int
        nlem_div)
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;

   long unsigned int counter=0;

   double *Y = lem.tmp_s1;

   for (it = cells.begin(); it != cells.end(); ++it)
   {

     if (counter<nlem_div) // products
     {
       gas->setState_TPY(1800, Cantera::OneAtm,
               "O2:7.826E-2, N2:7.3638E-1, H2O:6.546E-2, CO2:1.199E-1");
     }
     else // reactants
     {
       gas->setState_TPY(300, Cantera::OneAtm, "C3H8:4E-2, O2:2.236E-1, N2:7.364E-1");
     }

     (*it).p   = gas->pressure();
     (*it).T   = gas->temperature();
     (*it).rho = gas->density();
     gas->getMassFractions((*it).Y);
     (*it).m = gas->density() * (*it).dx;
     (*it).E = (gas->enthalpy_mass() - (*it).p / (*it).rho) * (*it).m;

     counter++;

   }

}



/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::initializeDataOnLineFromLES(double T_in, double P_in, std::vector<double> Y_in)
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;

   double T0=T_in;
   double P0=P_in;
   std::vector<double> Ycop = Y_in;

   gas->setState_TPY(T0,P0, Ycop.data());
   for (it = cells.begin(); it != cells.end(); ++it)
   {

     (*it).p   = gas->pressure();
     (*it).T   = gas->temperature();
     (*it).rho = gas->density();
     gas->getMassFractions((*it).Y);

     (*it).m = gas->density() * (*it).dx;
   }
   setEnergy();

}

void LEMLINE::setLEMT(double T_in)
{
   std::list<CellData>::iterator it;
   for (it = cells.begin(); it != cells.end(); ++it)
   {

     (*it).T   = T_in;
   }

}


/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::eddyEvent()
{

  int i0  = eddyLocation();

  //   cout << "finding eddy size" << std::endl;
  int nEddy = eddySize();

  //   cout << "doing triplet map" << std::endl;

  int nocell=cells.size();
  if ((i0+nEddy) > nocell-1)
      // nlem-1 // same is the condition for shifting elements back defined in
      // function tripletMap(i0, nEddy)
  {
  }
  else
  {
    tripletMap(i0, nEddy);
  }

}

/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
int LEMLINE::eddyLocation(void)
{
  nlem = cells.size();
  int i =  rand() % nlem + 1;
  return i;
}

/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
int LEMLINE::eddySize(void)
{

  const double F = (double) rand() / RAND_MAX;

//   cout << "random number is" << F << std::endl;
//   cout << "kalmogorov scale is" << eta << std::endl; // 'eta' is changing, while
//   'lt' remains the same as it is always the cube root of volume of LES cell

  const double a = -5./3.;
  const double eta_a = pow(eta, a);
  const double tmp = F * (pow(lt,a) - eta_a) + eta_a;
  const double lEddy = pow(tmp,-3./5.);
  nlem = cells.size();

  int nEddy;
  if(lEddy <= lt)
  {
    const double dx = lengthOfLine() / nlem;
    // change to this line because neddy is getting more than actual nlem on
    // line
    //   cout << "dx for eddy size is " << dx << std::endl;
    nEddy = (int)(lEddy / (3. * dx)) * 3;
    //   cout << "eddy size in lem cells is" << nEddy << std::endl;

    if (nEddy < 6)
    nEddy = 6;
  }
  else
  {
    nEddy =10000; // a very high value
    //   cout << "Info: eddy size is bigger than integral length scale" << std::endl;
  }

  return nEddy;
}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::tripletMap(
    int unsigned iStart,
    const unsigned int nEddy)
{

  if(cells.size() < 6) return;
  noofEddies++;
  double EddyLengthnow = 0.0;
  double dxnow = lengthOfLine() / nlem;
//  EddyLengthnow = (nEddy * (3. * dxnow)) / 3;
  EddyLengthnow = nEddy * dxnow;
  Total_EddiesLength += EddyLengthnow;
//   cout << "inputs to triplet map are: eddy start location:" << iStart << "
//   ,no of eddies:"  << nEddy << std::endl;
  std::list<CellData>::iterator it, it1, it2, it3, itStart, itEnd;
  const unsigned int nTriplet3  = nEddy / 3;
  const unsigned int nTriplet   = nTriplet3 * 3;

  // shift elements from front to back in case the eddy extends
  // the right end of the domain
  bool shift = false;
  int nShift = 0;
  if (iStart + nTriplet > cells.size() - 1)
  {
    shift = true;
    nShift = iStart + nTriplet - cells.size() - 1;
    std::list<CellData>::iterator i0 = cells.begin();
    std::list<CellData>::iterator i1 = i0;
    for (int i = 0; i < nShift; i++) ++i1;
    cells.splice(cells.end(), cells, i0, i1);
  }
  iStart -= nShift;
//  cout << "shifting ends" << std::endl;
  itStart = cells.begin();
  for(unsigned int i = 1; i < iStart; i++) ++itStart;
  itEnd = itStart;
  for(unsigned int i = 1; i <= nTriplet; i++) ++itEnd;


  it1 = itStart; ++it1;
  it2 = it1;

  for (unsigned int i = 1; i < nTriplet3; i++)
  {
    it2++, it2++;
    it3 = it2; ++it3;
    cells.splice(it1, cells, it2);
    it2 = it3;
  }
// cout << "1st for loop end" << std::endl;

  it1 = itEnd; --it1;
  it2 = it1; --it2;
  for (unsigned int i = 1; i < nTriplet3; i++)
  {
    --it2;
    it3 = it2; --it3;
    cells.splice(it1, cells, it2);
    // ++it3;
    it2 = it3;
    --it1;
  }
// cout << "2nd for loop end" << std::endl;

  it3 = it1; it1 = it2; it2 = it3;
  for (unsigned int i = 1; i <= nTriplet3 / 2; i++)
  {
    it3 = it1; it3++;
    cells.splice(it2, cells, it1);
    --it1; --it2;
    cells.splice(it3, cells, it1);
    ++it1;
  }
// cout << "3rd for loop end" << std::endl;
   // shifting back elements
//    cout << "shifting back elements" << std::endl;
   if (shift)
   {
      std::cout << "Warning: Shifting back elements in triplet map" << std::endl;
      std::list<CellData>::iterator i0 = cells.end();
      std::list<CellData>::iterator i1 = i0;
      for (int i = 0; i < nShift; i++) --i0;
      cells.splice(cells.begin(), cells, i0, i1);
   }
}

/*-----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/
int LEMLINE::getTotalNoofEddies()
{
  return noofEddies;
}

/*-----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/
double LEMLINE::getAverageEddiesLength()
{
  if (noofEddies > 0)
  {
    Average_EddiesLength = Total_EddiesLength / noofEddies;
  }
  else
  {
        Average_EddiesLength = 3.0387e-05; // ?? what is this constant?
  }
    return Average_EddiesLength;
}


/*-----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/
double LEMLINE::getLEMDomainLength()
{

  double LEMDomainLength = lengthOfLine();

  return LEMDomainLength;

}

/*-----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/
int LEMLINE::getLEMDomainCells()
{

  int LEMDomainCells = cells.size();

  return LEMDomainCells;

}

/**-----------------------------------------------------------------------------
 * @brief       advances the temperature and Fickian diffusion. Fickian
 * diffusivity coefficients are based on unity Lewis number assumption.
 *
 * @param       \input  dt      Diffusion time step
 * @author      M. Oevermann
 *----------------------------------------------------------------------------*/
void LEMLINE::advanceDiffusion(const double &dt)
{
  std::list<CellData>::iterator it, it1, it2;
  auto gas = lem.gas;

  if(cells.size()<3){ // TDMA solver needs at least 3 cells
    for (it = cells.begin(); it != cells.end(); ++it)
    {
      gas->setState_TPY((*it).T, (*it).p, (*it).Y);
      (*it).rho = gas->density();
      (*it).dx  = (*it).m / (*it).rho;
    }

  return;
  }

  for (it = cells.begin(); it != cells.end(); ++it){

    gas->setState_TPY((*it).T, (*it).p, (*it).Y);
    (*it).lambda = lem.tr->thermalConductivity();
    (*it).cp     = gas->cp_mass();

    (*it).a = (*it).b = (*it).c = (*it).d = 0.0;
  }

  // temperature equation
  it1 = it2 = cells.begin(); ++it2;

  while (it2 != cells.end() )
  {
    double dx     = 0.5 * ((*it1).dx     + (*it2).dx);
    double lambda = 0.5 * ((*it1).lambda + (*it2).lambda);
    double dtODx2_1 = dt / (dx * (*it1).dx);
    double dtODx2_2 = dt / (dx * (*it2).dx);

    double tmp1 = lambda / ((*it1).rho * (*it1).cp) * dtODx2_1;
    double tmp2 = lambda / ((*it2).rho * (*it2).cp) * dtODx2_2;

    (*it1).b += tmp1;
    (*it1).c -= tmp1;
    (*it2).b += tmp2;
    (*it2).a -= tmp2;

    ++it1; ++it2;
  }


  for (it = cells.begin(); it != cells.end(); ++it)
  {
    (*it).b += 1.0;
    (*it).d  = (*it).T;
  }


  solveWithThomasAlgorithm();
  for (it = cells.begin(); it != cells.end(); ++it)
  {
    (*it).T = (*it).x;
  }

  // -- Wall diffusion, Menon 8/6/2022
  if(wallPresent){
      it  = cells.begin();
      (*it).T = wallT;
  }


  // species diffusion, Le_s = 1 for now only
  for (int s = 0; s < lem.ns; s++)
  {
    for (it = cells.begin(); it != cells.end(); ++it)
    {
      (*it).a = (*it).b = (*it).c = (*it).d = 0.0;
    }
    it1 = it2 = cells.begin(); ++it2;
    while (it2 != cells.end() )
    {
      double dx   = 0.5 * (  (*it1).dx     + (*it2).dx);
      double rhoD = 0.5 * (  (*it1).lambda / (*it1).cp
                           + (*it2).lambda / (*it2).cp );
      double dtODx2_1 = dt / (dx * (*it1).dx);
      double dtODx2_2 = dt / (dx * (*it2).dx);

      double tmp1 = rhoD / ((*it1).rho) * dtODx2_1;
      double tmp2 = rhoD / ((*it2).rho) * dtODx2_2;

      (*it1).b += tmp1;
      (*it1).c -= tmp1;
      (*it2).b += tmp2;
      (*it2).a -= tmp2;

      ++it1; ++it2;
    }

    for (it = cells.begin(); it != cells.end(); ++it)
    {
      (*it).b += 1.0;
      (*it).d  = (*it).Y[s];
    }

    solveWithThomasAlgorithm();

    for (it = cells.begin(); it != cells.end(); ++it)
    {
      (*it).Y[s] = (*it).x;
    }
    for (it = cells.begin(); it != cells.end(); ++it)
    {
      gas->setState_TPY((*it).T, (*it).p, (*it).Y);
      (*it).rho = gas->density();
      (*it).dx  = (*it).m / (*it).rho;
    }
  }

}

/**-----------------------------------------------------------------------------
 * @brief       advances the temperature and Fickian diffusion. Fickian
 * diffusivity coefficients are based on unity Lewis number assumption.
 *
 * @param       \input  dt      Diffusion time step
 * @author      M. Oevermann
 *----------------------------------------------------------------------------*/
void LEMLINE::advanceDiffDiffusion(const double &dt)
{
  if(cells.size()<3) return;
  std::list<CellData>::iterator it, it1, it2;
  auto gas = lem.gas;

  for (it = cells.begin(); it != cells.end(); ++it){

    gas->setState_TPY((*it).T, (*it).p, (*it).Y);
    (*it).lambda = lem.tr->thermalConductivity();
    (*it).cp     = gas->cp_mass();
    lem.tr->getMixDiffCoeffsMass((*it).DD);


    (*it).a = (*it).b = (*it).c = (*it).d = 0.0;
  }

  // temperature equation
  it1 = it2 = cells.begin(); ++it2;

  while (it2 != cells.end())
  {
    double dx     = 0.5 * ((*it1).dx     + (*it2).dx);
    double lambda = 0.5 * ((*it1).lambda + (*it2).lambda);
    double dtODx2_1 = dt / (dx * (*it1).dx);
    double dtODx2_2 = dt / (dx * (*it2).dx);

    double tmp1 = lambda / ((*it1).rho * (*it1).cp) * dtODx2_1;
    double tmp2 = lambda / ((*it2).rho * (*it2).cp) * dtODx2_2;

    (*it1).b += tmp1;
    (*it1).c -= tmp1;
    (*it2).b += tmp2;
    (*it2).a -= tmp2;

    ++it1; ++it2;
  }

  for (it = cells.begin(); it != cells.end(); ++it)
  {
    (*it).b += 1.0;
    (*it).d  = (*it).T;
  }


  solveWithThomasAlgorithm();

  for (it = cells.begin(); it != cells.end(); ++it)
  {
    (*it).T = (*it).x;
  }

  // -- Wall diffusion, Menon 8/6/2022
  if(wallPresent){
      it  = cells.begin();
      (*it).T = wallT;
  }

  // species diffusion, differential diffusion enabled
  for (int s = 0; s < lem.ns; s++)
  {
    for (it = cells.begin(); it != cells.end(); ++it)
    {
      (*it).a = (*it).b = (*it).c = (*it).d = 0.0;
    }
    it1 = it2 = cells.begin(); ++it2;
    while (it2 != cells.end())
    {

      double dx   = 0.5 * (  (*it1).dx     + (*it2).dx);
      double rhoD = 0.5 * (  (*it1).DD[s] + (*it2).DD[s] );

      double dtODx2_1 = dt / (dx * (*it1).dx);
      double dtODx2_2 = dt / (dx * (*it2).dx);

    //double tmp1 = rhoD / ((*it1).rho) * dtODx2_1;
    //double tmp2 = rhoD / ((*it2).rho) * dtODx2_2;
      double tmp1 = rhoD /  dtODx2_1;
      double tmp2 = rhoD /  dtODx2_2;

      (*it1).b += tmp1;
      (*it1).c -= tmp1;
      (*it2).b += tmp2;
      (*it2).a -= tmp2;

      ++it1; ++it2;
    }

    for (it = cells.begin(); it != cells.end(); ++it)
    {
      (*it).b += 1.0;
      (*it).d  = (*it).Y[s];
    }

    solveWithThomasAlgorithm();

    for (it = cells.begin(); it != cells.end(); ++it)
    {
      (*it).Y[s] = (*it).x;
    }
    for (it = cells.begin(); it != cells.end(); ++it)
    {
      gas->setState_TPY((*it).T, (*it).p, (*it).Y);
      (*it).rho = gas->density();
      (*it).dx  = (*it).m / (*it).rho;
    }
  }

}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::solveWithThomasAlgorithm()
{
  std::list<CellData>::iterator it, itPrev, itNext;
  double tmp;

  it = cells.begin();
  (*it).c /= (*it).b;
  (*it).d /= (*it).b;
  itPrev = it;
  ++it;
  itNext = it;
  ++itNext;
  for (; itNext != cells.end(); ++itNext, ++it, ++itPrev)
  {
    tmp = (*it).b - (*it).a * (*itPrev).c;
    (*it).c /= tmp;
    (*it).d = ((*it).d - (*it).a * (*itPrev).d ) / tmp;
  }
  (*it).d = ((*it).d - (*it).a * (*itPrev).d )
          / ((*it).b - (*it).a * (*itPrev).c);

  it = itNext = cells.end();
  --it;
  (*it).x = (*it).d;

  do
  {
    --it;
    --itNext;
    (*it).x = (*it).d - (*it).c * (*itNext).x;
  }
  while (it != cells.begin());

}

/**
 * @brief       advances chemistry using CVODE
 *
 * @param       dt      \input          time step
 *
 * @author      M.Oevermann
 */
void LEMLINE::advanceSourceTerm(const double &dt)
{
  realtype tret;
  std::list<CellData>::iterator it;
  auto gas = lem.gas;
  Cantera::IdealGasConstPressureReactor& comb =  lem.combustor;
  Cantera::ReactorNet& rNet =  lem.sim;

    for (it = cells.begin(); it != cells.end(); ++it){
        //-- setup reactor
        gas->setState_TPY((*it).T,lem.pressure,(*it).Y);
        ySol[0] = gas->density();
     // ySol[1] = gas->enthalpy_mass()*gas->density();
        ySol[1] = gas->temperature();
        gas->getMassFractions(ySol.data()+2);
        rNet.setInitialTime(0.0);

        comb.setInitialVolume(1.0);
        comb.updateState(ySol.data());

        // -- react
        rNet.setInitialTime(0.0);
        rNet.advance(dt);

        // -- copy to domain

        (*it).T   = comb.temperature();
        (*it).p   = comb.pressure();
        (*it).rho = comb.density();
        gas->getMassFractions((*it).Y);

  //realtype *y = NV_DATA_S(lem.y);
  //y[IT]       = (*it).T;
  //y[IP]       = (*it).p;
  //y[IRHO]     = (*it).rho;

  //for (int s = 0; s < lem.ns; ++s) {
  //  y[IY0+s] = ((*it).Y)[s];
  //}

  //CVodeReInit(lem.cvode_mem, time, lem.y);
  //CVode(lem.cvode_mem, time+dt, lem.y, &tret, CV_NORMAL);

  //(*it).T   = y[IT];
  //(*it).p   = y[IP];
  //(*it).rho = y[IRHO];

  //for (int s = 0; s < lem.ns; ++s){

  //    ((*it).Y)[s] = y[IY0+s];
  //}

    (*it).dx = (*it).m / (*it).rho;
  }
}




/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
int chemistrySource(realtype _t, N_Vector _y, N_Vector _ydot, void *user_data)
{
  LEM *lem = (LEM *) user_data;
  auto gas = lem->gas;
  const int ns = lem->ns;

  realtype *y  = NV_DATA_S(_y);
  realtype *yp = NV_DATA_S(_ydot);

  const double T   = y[IT];
  const double p   = y[IP];
  const double rho = y[IRHO];
  const double *Y  = &(y[IY0]);

  gas->setState_TPY(T, p, Y);

  const double invRho = 1. / rho;
  const double M      = gas->meanMolecularWeight();
  const double R      = Cantera::GasConstant / M;
  const double cp     = gas->cp_mass();

  double *RR   = lem->tmp_s1;
  double *H_RT = lem->tmp_s2;

  lem->kin->getNetProductionRates(RR);
  gas->getEnthalpy_RT(H_RT);

  yp[IT]   = 0.0;
  yp[IP]   = 0.0;
  yp[IRHO] = 0.0;
  for(int s = 0; s < ns; ++s)
  {
 // double hs = H_RT[s] * T * Cantera::GasConstant;
    double hs = H_RT[s] * T * lem->Rs[s];
    double dYdt = lem->Ms[s] * RR[s] * invRho;
    yp[IY0+s] = dYdt;
    yp[IT] -= hs * dYdt;
    yp[IRHO] -= lem->Rs[s] * dYdt;
  }
  yp[IRHO] *= rho / R;
  yp[IT] /= cp;
  yp[IT] += yp[IP] / (rho * cp);

  yp[IRHO] += yp[IP] / (R * T) - rho / T * yp[IT];

  return(0);
}

// --- Means, medians and mins S.Arshad stuff {{{
/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getSpeciesMeanMixtureFractions()
{
auto gas = lem.gas;
std::list<CellData>::iterator it;
double rho_sum=0.0;

std::vector<double> Y_sum;
Y_sum.resize(lem.ns,0.0);

std::vector<double> Y_mean;
Y_mean.resize(lem.ns,0.0);

for (it = cells.begin(); it != cells.end(); ++it)
{
  rho_sum += (*it).rho * (*it).dx;
  for (unsigned int s=0; s< lem.ns; s++)
  {
    Y_sum[s]+= (*it).Y[s] * (*it).rho * (*it).dx;
  }
}

for(unsigned int s=0; s < lem.ns; s++)
{
  Y_mean[s]=Y_sum[s]/rho_sum;
}

return Y_mean;
}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanTemperature()
{
std::list<CellData>::iterator it;

double rho_sum = 0.0;
double T_sum   = 0.0;
double T_mean  = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  rho_sum += (*it).rho * (*it).dx;
  T_sum   += (*it).T   * (*it).rho * (*it).dx;
}

T_mean = T_sum/rho_sum;

return T_mean;
}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMinTemperature()
{
std::list<CellData>::iterator it;

double T_min  = 10000.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  T_min =std::min(T_min, (*it).T);
}

return T_min;
}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMaxTemperature()
{
std::list<CellData>::iterator it;
double T_max  = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  T_max = std::max(T_max, (*it).T);
}

return T_max;
}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanTemperatureCpWeighted()
{
auto gas = lem.gas;
std::list<CellData>::iterator it;

double rho_sum = 0.0;
double cp      = 0.0;
double T_sum   = 0.0;
double T_mean  = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  gas->setState_TPY((*it).T, (*it).p,(*it).Y);
  cp = gas->cp_mass();
  rho_sum += (*it).rho * cp * (*it).dx;
  T_sum += (*it).T * (*it).rho * cp * (*it).dx;
}

T_mean = T_sum/rho_sum;

return T_mean;
}


/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanKinematicViscosity()
// getMeanKinematicViscosity is only to be used by setdtStirr and cannot
// be used in the same way as, e.g., getMeanDensity or getMeanDynamicViscosity.
{
std::list<CellData>::iterator it;

double rho_sum  = 0.0;
double rho_mean = 0.0;

double nu_sum   = 0.0;
double nu_mean  = 0.0;

double dx_sum = 0.0;

double MeanKinematicViscosity = 0.0;
auto gas = lem.gas;

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum  += (*it).dx;
  rho_sum += (*it).rho * (*it).dx;
  gas->setState_TPY((*it).T, (*it).p, (*it).Y);

  double nu_l = 0.0;
  nu_l = lem.tr->viscosity();
  double nu_Lem = 0.0;
  nu_Lem = nu_l * (*it).rho * (*it).dx;
  nu_sum += nu_Lem;
}

rho_mean = rho_sum/dx_sum;

nu_mean = nu_sum/rho_sum;

MeanKinematicViscosity = nu_mean/rho_mean;

//cout << "mu1= " << nu_mean << std::endl;
//cout << "nu= " << MeanKinematicViscosity << std::endl;

return MeanKinematicViscosity;

}


/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanDynamicViscosity()
{
auto gas = lem.gas;
//Transport* trmix = newTransportMgr("Mix", &gas);
auto tr = lem.tr;
std::list<CellData>::iterator it;

double rho_sum  = 0.0;

double mu_sum   = 0.0;
double mu_mean  = 0.0;

double dx_sum = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum  += (*it).dx;
  rho_sum += (*it).rho * (*it).dx;

  //double mu_l = 0.0;
  //mu_l = lem.tr->viscosity();
  //double mu_Lem = 0.0;
  //mu_Lem = mu_l * (*it).rho * (*it).dx;
  //mu_sum += mu_Lem;

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);
  double mu_cell = tr->viscosity() * (*it).rho * (*it).dx;
  mu_sum += mu_cell;

}

mu_mean = mu_sum/rho_sum;

//cout << "mu2= " << mu_mean << std::endl;

return mu_mean;

}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanGasConstant()
{
auto gas = lem.gas;
//Transport* tr = lem.tr;
std::list<CellData>::iterator it;

double rho_sum  = 0.0;

double R_sum   = 0.0;
double R_mean  = 0.0;

double dx_sum = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum  += (*it).dx;
  rho_sum += (*it).rho * (*it).dx;

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);
  double R = Cantera::GasConstant / gas->meanMolecularWeight();
  double R_cell = R * (*it).rho * (*it).dx;
  R_sum += R_cell;

}

R_mean = R_sum/rho_sum;

//cout << "R_mean= " << R_mean << std::endl;

return R_mean;

}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanEnthalpy()
{
auto gas = lem.gas;
//Transport* tr = lem.tr;
std::list<CellData>::iterator it;

double rho_sum  = 0.0;

double h_sum   = 0.0;
double h_mean  = 0.0;

double dx_sum = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum  += (*it).dx;
  rho_sum += (*it).rho * (*it).dx;

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);
  double h = gas->enthalpy_mass();
  double h_cell = h * (*it).rho * (*it).dx;
  h_sum += h_cell;
}

h_mean = h_sum/rho_sum;

//cout << "h_mean= " << h_mean << std::endl;

return h_mean;

}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanFormationEnthalpy()
{
auto gas = lem.gas;
//Transport* tr = lem.tr;
std::list<CellData>::iterator it;

double rho_sum  = 0.0;

double h_sum   = 0.0;
double h_mean  = 0.0;

double dx_sum = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum  += (*it).dx;
  rho_sum += (*it).rho * (*it).dx;

  gas->setState_TPY(298.15, 101325, (*it).Y);
  double h = gas->enthalpy_mass();
  double h_cell = h * (*it).rho * (*it).dx;
  h_sum += h_cell;
}

h_mean = h_sum/rho_sum;

//cout << "h_mean= " << h_mean << std::endl;

return h_mean;

}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanDensity()
{
std::list<CellData>::iterator it;

double rho_sum  = 0.0;
double rho_mean = 0.0;

double dx_sum = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum  += (*it).dx;
  rho_sum += (*it).rho * (*it).dx;
}

rho_mean = rho_sum/dx_sum;

return rho_mean;

}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanPressure()
{
std::list<CellData>::iterator it;

double p_sum  = 0.0;
double p_mean = 0.0;

double dx_sum = 0.0;

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum  += (*it).dx;
  p_sum   += (*it).p * (*it).dx;
}

p_mean = p_sum/dx_sum;

return p_mean;

}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
//double LEMLINE::getMeanDynamicViscosity()
//{
//
//  double mu_mean  = 0.0;
//
//  double nu = getMeanKinematicViscosity();
//
//
//  double rho = getMeanDensity();
//
//  mu_mean = rho * nu;
//
//  return mu_mean;
//
//}


/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getSpeciesMeanNetProductionRate()
{

auto gas = lem.gas;
std::list<CellData>::iterator it;

double* Ms = lem.Ms;
double* RR = lem.tmp_s1;

double dx_sum=0.0;

std::vector<double> YSource_sum;
YSource_sum.resize(lem.ns,0.0);

std::vector<double> YSource_mean;
YSource_mean.resize(lem.ns,0.0);

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum += (*it).dx;

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);

  lem.kin->getNetProductionRates(RR);

  for (unsigned int s=0; s< lem.ns; s++)
  {
    YSource_sum[s]+= Ms[s] * RR[s] * (*it).dx;
  }

}

for(unsigned int s=0; s < lem.ns; s++)
{
  YSource_mean[s]=YSource_sum[s]/dx_sum;
}

return YSource_mean;
}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMeanHeatProductionRate()
{

auto gas = lem.gas;
std::list<CellData>::iterator it;

double* Ms = lem.Ms;
double* Rs = lem.Rs;
double* RR = lem.tmp_s1;
double* H_RT = lem.tmp_s2;

double dx_sum=0.0;
double hsSource_sum=0;
double hsSource_mean=0;

std::vector<double> hformation;
hformation.resize(lem.ns,0.0);

gas->setState_TP(298.0, 101325.0);
gas->getEnthalpy_RT(H_RT);
for (unsigned int s=0; s< lem.ns; s++)
{
  hformation[s] = H_RT[s] * Rs[s] * 298.0;
}

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum += (*it).dx;

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);

  lem.kin->getNetProductionRates(RR);
  //gas->getEnthalpy_RT(H_RT);

  for (unsigned int s=0; s< lem.ns; s++)
  {
    double speciesProd = Ms[s] * RR[s];
    hsSource_sum+= - hformation[s] * speciesProd * (*it).dx;
  }

}

hsSource_mean=hsSource_sum/dx_sum;

return hsSource_mean;
}

/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::getMedianHeatProductionRate()
{

auto gas = lem.gas;
std::list<CellData>::iterator it;

double* Ms = lem.Ms;
double* Rs = lem.Rs;
double* RR = lem.tmp_s1;
double* H_RT = lem.tmp_s2;

std::vector<double> hsSource;
hsSource.resize(nlem,0.0);

double hsSourcei;
double hsSource_median;

std::vector<double> hformation;
hformation.resize(lem.ns,0.0);

gas->setState_TP(298.0, 101325.0);
gas->getEnthalpy_RT(H_RT);
for (unsigned int s=0; s< lem.ns; s++)
{
  hformation[s] = H_RT[s] * Rs[s] * 298.0;
}

int ilem=0;
for (it = cells.begin(); it != cells.end(); ++it)
{
  gas->setState_TPY((*it).T, (*it).p, (*it).Y);
  lem.kin->getNetProductionRates(RR);
  //gas->getEnthalpy_RT(H_RT);

  hsSourcei=0.0;

  for (unsigned int s=0; s< lem.ns; s++)
  {
    double speciesProd = Ms[s] * RR[s];
    hsSourcei += - hformation[s] * speciesProd;
    //cout << "Specie " << s << " hformation= " << hformation[s] << " Rs= " << Rs[s] << std::endl;
  }

  hsSource[ilem] = hsSourcei;

  ilem++;

}

hsSource_median=calculateMedian(hsSource);

return hsSource_median;
}


/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getSpeciesAndHeatMeanNetProductionRate()
{

auto gas = lem.gas;
std::list<CellData>::iterator it;

double* Ms = lem.Ms;
double* Rs = lem.Rs;
double* RR = lem.tmp_s1;
double* H_RT = lem.tmp_s2;

double dx_sum=0.0;
double hsSource_sum=0;

std::vector<double> YSource_sum;
YSource_sum.resize(lem.ns+1,0.0);
std::vector<double> YSource_mean;
YSource_mean.resize(lem.ns+1,0.0);

std::vector<double> hformation;
hformation.resize(lem.ns,0.0);
gas->setState_TP(298.0, 101325.0);
gas->getEnthalpy_RT(H_RT);
for (unsigned int s=0; s< lem.ns; s++)          // why recalculate this every single time?
{
  hformation[s] = H_RT[s] * Rs[s] * 298.0;
}

for (it = cells.begin(); it != cells.end(); ++it)
{
  dx_sum += (*it).dx;

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);

  lem.kin->getNetProductionRates(RR);

  for (unsigned int s=0; s< lem.ns; s++)
  {
    YSource_sum[s]+= Ms[s] * RR[s] * (*it).dx; // species

    hsSource_sum  += -hformation[s] * Ms[s] * RR[s] * (*it).dx;
    //cout << "Specie " << s << " hformation= " << hformation[s] << " Rs= " << Rs[s] << std::endl;
  }
}

//cout << "hsSource_mean 1 =" << hsSource_sum/dx_sum << std::endl;

YSource_sum[lem.ns]=hsSource_sum;

for(unsigned int s=0; s < lem.ns+1; s++)
{
  YSource_mean[s]=YSource_sum[s]/dx_sum;
}

//cout << "hsSource_mean 2 =" << YSource_mean[lem.ns] << std::endl;

return YSource_mean;

}


/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getSpeciesAndHeatMedianNetProductionRate()
{

auto gas = lem.gas;
std::list<CellData>::iterator it;

double* Ms = lem.Ms;
double* Rs = lem.Rs;
double* RR = lem.tmp_s1;
double* H_RT = lem.tmp_s2;

int ns=lem.ns;

double YSourcei;
double hSourcei;

std::vector< std::vector<double> > YSource(ns+1, std::vector<double>(nlem));

std::vector<double> YSource_median;
YSource_median.resize(ns+1,0.0);

std::vector<double> hformation;
hformation.resize(lem.ns,0.0);
gas->setState_TP(298.0, 101325.0);
gas->getEnthalpy_RT(H_RT);
for (unsigned int s=0; s< lem.ns; s++)
{
  hformation[s] = H_RT[s] * Rs[s] * 298.0;
}

int ilem=0;
for (it = cells.begin(); it != cells.end(); ++it)
{
  hSourcei=0;

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);
  lem.kin->getNetProductionRates(RR);

  for (unsigned int s=0; s< lem.ns; s++)
  {
    YSourcei = Ms[s] * RR[s]; // species

    YSource[s][ilem] = YSourcei;

    hSourcei += -hformation[s] * YSourcei;
  }

  YSource[ns][ilem] = hSourcei;

  ilem++;
}

for (unsigned int s=0; s<= lem.ns; s++)
{
/*      ilem=0;
  for (it = cells.begin(); it != cells.end(); ++it)
  {
    cout << "ilem=" << ilem << " s=" << s << " YSource=" << YSource[s][ilem] << std::endl;
  ilem++;
  }
*/
  YSource_median[s] = calculateMedian(YSource[s]);

//      cout << "      s=" << s << " YSource_median=" << YSource_median[s] << std::endl;

}

return YSource_median;

}


/*-----------------------------------------------------------------------------
*
* --------------------------------------------------------------------------*/
double LEMLINE::calculateMedian(std::vector<double> V)
{

double median;

std::size_t size = V.size();

std::sort(V.begin(), V.end());

if (size % 2 == 0)
{
 median = (V[size / 2 - 1] + V[size / 2]) / 2;
}
else
{
 median = V[size / 2];
}

return median;

}

// ---  }}}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getLEMcellsTemperature()
{

  auto gas = lem.gas;
  std::list<CellData>::iterator it;

  std::vector<double> tlemcellall;
  tlemcellall.resize(nlem, 285.0);

  int ilem = 0;
  for (it = cells.begin(); it != cells.end(); ++it)
  {
    tlemcellall[ilem] = (*it).T;

    ilem++;
  }

  return tlemcellall;
}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
double LEMLINE::getMaxYfuel(int fuelindex)
{
  std::list<CellData>::iterator it;

  double Yfuel_max  = 1e-80;

  for (it = cells.begin(); it != cells.end(); ++it)
  {
      Yfuel_max =std::max(Yfuel_max, (*it).Y[fuelindex]);
  }

  return Yfuel_max;
}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
double LEMLINE::internalEnergyOnLine()
{
  auto gas = lem.gas;
  std::list<CellData>::iterator it;
  double E = 0.0;
  for (it = cells.begin(); it != cells.end(); ++it)
  {
    gas->setState_TPY((*it).T, (*it).p, (*it).Y);
    double E_cell = (gas->enthalpy_mass() - (*it).p / (*it).rho) * (*it).m;
    E += E_cell;
    (*it).E = E_cell;
  }
  return E;
}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
double LEMLINE::lengthOfLine()
{
  std::list<CellData>::iterator it;
  double l = 0.0;
  for (it = cells.begin(); it != cells.end(); ++it)
  {
         (*it).dx = (*it).m/(*it).rho;
         l += (*it).dx;
  }

  return l;
}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
double LEMLINE::massOnLine()
{
  std::list<CellData>::iterator it;
  double mass = 0.0;
  for (it = cells.begin(); it != cells.end(); ++it)
  {
    mass += (*it).m;
  }
  return mass;
}


/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::correctExpansion()
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;
   double l = 0.0;
   double cp = 0.0;
   double cv = 0.0;

   // find the current state on LEM line, length and mass averaged cp and cv
   for (it = cells.begin(); it != cells.end(); it++)
   {
      const double m = (*it).m;
      gas->setState_TPY((*it).T, (*it).p,(*it).Y);
      l  += m / (*it).rho;
      cp += gas->cp_mass() * m;
      cv += gas->cv_mass() * m;
   }

   // mean ratio of specific heats
   double gamma = cp / cv; // expansion factor of the LEM domain
   double alpha = l / length;  // ratio of current length to original

   // pressure after isentropic compression
   // double p = gas->pressure() * pow(alpha, gamma);

   double e = 0.0;
   double dedp = 0.0;
   for (it = cells.begin(); it != cells.end(); it++)
   {
     // temperature in cell after isentropic compression
     double T = (*it).T * pow(alpha, (gamma-1.));
     gas->setState_TRY(T, (*it).rho * alpha, (*it).Y);
     (*it).p = gas->pressure();
     (*it).T = gas->temperature();
     (*it).rho = gas->density();
     e += (gas->enthalpy_mass() - (*it).p / (*it).rho) * (*it).m;
     dedp += (*it).m / ((*it).rho * (gas->cp_mass() / gas->cv_mass() - 1.));
     (*it).dx = (*it).m / (*it).rho;
   }

   // pressure iteration to match the current internal energy (e) to the
   // initial (E) of LEM line
   double delta_e = Energy - e;
   double delta_p = delta_e / dedp;
   while (fabs(delta_e / Energy) > 1e-8)
   {
     e = 0.0;
     dedp = 0.0;
     for (it = cells.begin(); it != cells.end(); it++)
     {
       (*it).p += delta_p;
       gas->setMassFractions((*it).Y);
       const double M = gas->meanMolecularWeight();
       const double R = Cantera::GasConstant / M;
       (*it).T = (*it).p / ((*it).rho * R);
       gas->setState_TRY((*it).T, (*it).rho, (*it).Y);
       double E_cell = (gas->enthalpy_mass() - (*it).p / (*it).rho) * (*it).m;
       e += E_cell;
       (*it).E = E_cell;
       dedp += (*it).m / ((*it).rho * (gas->cp_mass() / gas->cv_mass() - 1.));
     }
     delta_e = Energy - e;
     delta_p = delta_e / dedp;
   }

   it = cells.begin();
}

/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
void LEMLINE::CorrectionLocalPressure( double Pnew)
{
   std::list<CellData>::iterator it;
   for (it = cells.begin(); it != cells.end(); ++it)
   {
   //   cout << "on lem cell number" << std::distance(cells.begin(), it) << "Pnew is" << Pnew << "and pressure is " << (*it).p << std::endl;
     (*it).p = Pnew;
   }

}

/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
void LEMLINE::CorrectionGlobalPressure( double Pnew)
{
   std::list<CellData>::iterator it;
   for (it = cells.begin(); it != cells.end(); ++it)
   {
     double FactorCorr = 0.0;
     FactorCorr =  Pnew / (*it).p;
     (*it).p = Pnew;
     (*it).T *= FactorCorr;
   }

}

/*------------------------------------------------------------------------------
 * @author      S.Arshad
 *
 *----------------------------------------------------------------------------*/

void LEMLINE::CorrectionGlobalPressureAssumingAdiabaticProcess( double Pnew)
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;

   for (it = cells.begin(); it != cells.end(); ++it)
   {
     gas->setState_TPY((*it).T, (*it).p,(*it).Y);
     double cp = gas->cp_mass();
     double cv = gas->cv_mass();
     double gamma = cp/cv;
     double exponent = (gamma-1)/gamma;
     double FactorCorr = 0.0;
     FactorCorr =  Pnew / (*it).p;
     (*it).p = Pnew;
     (*it).T *= pow(FactorCorr,exponent);
     (*it).rho = gas->density();
   }

}


/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
double LEMLINE::CorrectLEMTemperatures(double Tles, double Tlem, double Tmin, double Tmax)
{
   std::list<CellData>::iterator it;

//   double Tmin=273;
//   double Tmax=1714;

   double Terror = Tles - Tlem;

   double rho_sum=0.0;
   double T_sum=0.0;
   double T_mean=0.0;

   for (it = cells.begin(); it != cells.end(); ++it)
   {
//     cout << "Told=" << (*it).T << std::endl;
//     (*it).T += Terror;
//     (*it).T = max((*it).T, Tmin);
//     (*it).T = min((*it).T, Tmax);
     (*it).T = std::min( std::max( (*it).T + Terror , Tmin)  , Tmax);
//     cout << "Tnew=" << (*it).T << std::endl;

     rho_sum += (*it).rho * (*it).dx;
     T_sum += (*it).T * (*it).rho * (*it).dx;
   }

    T_mean = T_sum / rho_sum;

    //double Terror_new = Tles - T_mean;

    return T_mean;

}


/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
double LEMLINE::CorrectLEMTemperaturesByFactor(double Tles, double Tlem, double Tmin, double Tmax)
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;

   double rho_sum=0.0;
   double cp= 0.0;
   double T_sum=0.0;
   double T_mean=0.0;

   for (it = cells.begin(); it != cells.end(); ++it)
   {
     //cout << "Told=" << (*it).T << std::endl;
     (*it).T *= Tles/Tlem;
     (*it).T = std::max((*it).T, Tmin);
     (*it).T = std::min((*it).T, Tmax);
     //(*it).T = std::min( std::max( (*it).T, Tmin)  , Tmax);
     //cout << "Tnew=" << (*it).T << std::endl;

     gas->setState_TPY((*it).T, (*it).p,(*it).Y);
     cp = gas->cp_mass();
     rho_sum += (*it).rho * cp * (*it).dx;
     T_sum += (*it).T * (*it).rho * cp * (*it).dx;
   }

    T_mean = T_sum / rho_sum;

    //double Terror_new = Tles - T_mean;

    return T_mean;

}

/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
double LEMLINE::CorrectLEMSingleSpecies(int s, double Yles, double Ylem)
{
   std::list<CellData>::iterator it;

   double Yerror = Yles - Ylem;

   double rho_sum=0.0;
   double Y_sum=0.0;
   double Y_mean=0.0;

   for (it = cells.begin(); it != cells.end(); ++it)
   {
//      cout << "specie= " << s << ": Yold=" << (*it).Y[s] << std::endl;
      (*it).Y[s] += Yerror;
      (*it).Y[s] = std::max((*it).Y[s], 0.0);
      (*it).Y[s] = std::min((*it).Y[s], 1.0);
//      cout << "        " << s << ": Ynew=" << (*it).Y[s] << std::endl;

      rho_sum += (*it).rho * (*it).dx;
      Y_sum += (*it).Y[s] * (*it).rho * (*it).dx;
    }

    Y_mean = Y_sum / rho_sum;

    //double Yerror_new = Yles - Y_mean;

    return Y_mean;
}



/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
void LEMLINE::CorrectLEMDensities()
{
   auto gas = lem.gas;
   std::list<CellData>::iterator it;

   for (it = cells.begin(); it != cells.end(); ++it)
   {
     gas->setState_TPY((*it).T, (*it).p,(*it).Y);
     (*it).rho = gas->density();
   }

}

/*-----------------------------------------------------------------------------
 *            		Binary RESTART LEM lines
 * --------------------------------------------------------------------------*/
void LEMLINE::binary_write(std::ofstream& o)
{
  std::list<CellData>::iterator it;

  o.write((char*)&length,sizeof(length));
  o.write((char*)&Energy, sizeof(Energy));
  o.write((char*)&time, sizeof(time));

  for (it = cells.begin(); it != cells.end(); it++)
  {
    o.write((char*)&(*it).dx, sizeof((*it).dx));
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    o.write((char*)&(*it).T, sizeof((*it).T));
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    o.write((char*)&(*it).m, sizeof((*it).m));
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    o.write((char*)&(*it).E, sizeof((*it).E));
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    o.write((char*)&(*it).rho, sizeof((*it).rho));
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    o.write((char*)&(*it).p, sizeof((*it).p));
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    for (int s = 0; s < (*it).ns; s++)
    {
          o.write((char*)&(*it).Y[s], sizeof((*it).Y[s]));
    }
  }
}



/*-----------------------------------------------------------------------------
 *            		Output LEM data for debugging/LEM performance
 * --------------------------------------------------------------------------*/
void LEMLINE::LEMoutput_write(std::ofstream& o, double timeNow, int noLES)
{
  std::list<CellData>::iterator it;

  double x= 0.0;

  for (it = cells.begin(); it != cells.end(); it++)
  {
    x += (*it).dx;
    o << timeNow << "    " << noLES << "    " << x << "    " << (*it).T << std::endl ;
  }

}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::binary_read(std::ifstream& i)
{

   std::list<CellData>::iterator it;

   double length_res = 0.0;
   double Energy_res = 0.0;
   double time_res = 0.0;
   std::vector<double> dx_res;
   std::vector<double> T_res;
   std::vector<double> m_res;
   std::vector<double> E_res;
   std::vector<double> rho_res;
   std::vector<double> p_res;
   std::vector<std::vector<double> > Y_res;

   dx_res.resize( cells.size(), 0.0);
   T_res.resize( cells.size(), 0.0);
   m_res.resize( cells.size(), 0.0);
   E_res.resize( cells.size(), 0.0);
   rho_res.resize( cells.size(), 0.0);
   p_res.resize( cells.size(), 0.0);
   Y_res.resize( cells.size());
   for (it = cells.begin(); it != cells.end(); it++)
   {
    Y_res[std::distance(cells.begin(), it)].resize( (*it).ns, 0.0);
   }

   i.read((char*)&length_res,sizeof(length_res));
   i.read((char*)&Energy_res, sizeof(Energy_res));
   i.read((char*)&time_res, sizeof(time_res));

  for (it = cells.begin(); it != cells.end(); it++)
  {
   i.read((char*)&dx_res[std::distance(cells.begin(), it)], sizeof(dx_res[std::distance(cells.begin(), it)]));
//    cout << "the read dx of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << dx_res[std::distance(cells.begin(), it)] << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
   i.read((char*)&T_res[std::distance(cells.begin(), it)], sizeof(T_res[std::distance(cells.begin(), it)]));
//    cout << "the read temperature of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << T_res[std::distance(cells.begin(), it)] << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
   i.read((char*)&m_res[std::distance(cells.begin(), it)], sizeof(m_res[std::distance(cells.begin(), it)]));
//    cout << "the read mass of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << m_res[std::distance(cells.begin(), it)] << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
   i.read((char*)&E_res[std::distance(cells.begin(), it)], sizeof(E_res[std::distance(cells.begin(), it)]));
//    cout << "the read energy of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << E_res[std::distance(cells.begin(), it)] << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
   i.read((char*)&rho_res[std::distance(cells.begin(), it)], sizeof(rho_res[std::distance(cells.begin(), it)]));
//    cout << "the read density of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << rho_res[std::distance(cells.begin(), it)] << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
   i.read((char*)&p_res[std::distance(cells.begin(), it)], sizeof(p_res[std::distance(cells.begin(), it)]));
//    cout << "the read pressure of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << p_res[std::distance(cells.begin(), it)] << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    for (int s = 0; s < (*it).ns; s++)
    {
      i.read((char*)&Y_res[std::distance(cells.begin(), it)][s], sizeof(Y_res[std::distance(cells.begin(), it)][s]));
//       cout << "the read mass fraction of specie " << s << " of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << Y_res[std::distance(cells.begin(), it)][s] << std::endl;
    }
  }

//   std::clock_t startini;
//   startini = std::clock();
  initializeRestartData(length_res, Energy_res, time_res, dx_res, T_res, m_res, E_res, rho_res, p_res, Y_res);
//   cout << "Time for function initializeRestartData() is: " << (std::clock() - startini) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::binary_read_wo_ini(std::ifstream& i)
{
  //   std::clock_t startrb;
  //   startrb = std::clock();

   std::list<CellData>::iterator it;

   i.read((char*)&length,sizeof(length));
   i.read((char*)&Energy, sizeof(Energy));
   i.read((char*)&time, sizeof(time));
//   cout << "the read length from binary file is " << length << std::endl;

   for (it = cells.begin(); it != cells.end(); it++)
   {
      i.read((char*)&(*it).dx, sizeof((*it).dx));
//    cout << "the read dx of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << (*it).dx << std::endl;
   }

   for (it = cells.begin(); it != cells.end(); it++)
   {
      i.read((char*)&(*it).T, sizeof((*it).T));
//    cout << "the read temperature of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << (*it).T << std::endl;
   }

   for (it = cells.begin(); it != cells.end(); it++)
   {
      i.read((char*)&(*it).m, sizeof((*it).m));
//    cout << "the read mass of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << (*it).m << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    i.read((char*)&(*it).E, sizeof((*it).E));
//    cout << "the read energy of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << (*it).E << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    i.read((char*)&(*it).rho, sizeof((*it).rho));
//    cout << "the read density of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << (*it).rho << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
   i.read((char*)&(*it).p, sizeof((*it).p));
//    cout << "the read pressure of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << (*it).p << std::endl;
  }

  for (it = cells.begin(); it != cells.end(); it++)
  {
    for (int s = 0; s < (*it).ns; s++)
    {
      i.read((char*)&(*it).Y[s], sizeof((*it).Y[s]));
//       cout << "the read mass fraction of specie " << s << " of lem cell number " <<  std::distance(cells.begin(), it) <<  " from binary file is " << (*it).Y[s] << std::endl;
    }
  }
//     cout << "Time for reading binary data is: " << (std::clock() - startrb) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
}


/*-----------------------------------------------------------------------------
 *
 *---------------------------------------------------------------------------*/
bool LEMLINE::noLEMcells()
{
  std::list<CellData>::iterator it;
  int count =0;

  for (it = cells.begin(); it != cells.end(); ++it)
  {
    count++;
  }

  if(count!=nlem)
  {
      nlem = count;
      return true;

  }
  else
  {
        return false;
  }
}

/*-----------------------------------------------------------------------------
 *			Assignment Operator Overloading
 * --------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::operator=(const LEMLINE &lemline)
{

    nlem = lemline.nlem;
    length = lemline.length;
    Energy = lemline.Energy;
    counter = lemline.counter;
    time = lemline.time;
    timeNextTripletMap = lemline.timeNextTripletMap;
    dtLEM = lemline.dtLEM;
    dtStirr = lemline.dtStirr;
    ut = lemline.ut;
    lt = lemline.lt;
    eta = lemline.eta;
    N_ETA = lemline.N_ETA;
    C_LAMBDA = lemline.C_LAMBDA;

    //cout << "entering Celldata = operator from LEMLINE" << std::endl;
    cells.assign(lemline.cells.begin(),lemline.cells.end());
}

/*-----------------------------------------------------------------------------
 *			Stream RESTART LEM lines
 * --------------------------------------------------------------------------*/
//{{{

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const LEMLINE &lem)
{

  std::list<CellData> c = lem.cells; // making copy of original lem.cells
  std::list<CellData>::iterator it;

  os << lem.length << std::endl;
  os << lem.Energy << std::endl;
  os << lem.time << std::endl;

  for (it = c.begin(); it != c.end(); it++)
  {
    os << (*it).dx << " ";
  }
  os << std::endl;

  for (it = c.begin(); it != c.end(); it++)
  {
    os << (*it).T << " ";
  }
  os << std::endl;

  for (it = c.begin(); it != c.end(); it++)
  {
    os << (*it).m << " ";
  }
  os << std::endl;

  for (it = c.begin(); it != c.end(); it++)
  {
    os << (*it).E << " ";
  }
  os << std::endl;

  for (it = c.begin(); it != c.end(); it++)
  {
    os << (*it).rho << " ";
  }
  os << std::endl;

  for (it = c.begin(); it != c.end(); it++)
  {
    os << (*it).p << " ";
  }
  os << std::endl;

  for (it = c.begin(); it != c.end(); it++)
  {
    for (int s = 0; s < (*it).ns; s++)
    {
        os << (*it).Y[s] << " ";
    }
  }
  os << std::endl;

  return os;
}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/

std::istream& operator>> (std::istream& is, LEMLINE &lem)
{

   std::list<CellData> c = lem.cells;
   std::list<CellData>::iterator it;

   double length_res = 0.0;
   double Energy_res = 0.0;
   double time_res = 0.0;

   std::vector<double> dx_res;
   std::vector<double> T_res;
   std::vector<double> m_res;
   std::vector<double> E_res;
   std::vector<double> rho_res;
   std::vector<double> p_res;
   std::vector<std::vector<double> > Y_res;

   dx_res.resize( lem.cells.size(), 0.0);
   T_res.resize( lem.cells.size(), 0.0);
   m_res.resize( lem.cells.size(), 0.0);
   E_res.resize( lem.cells.size(), 0.0);
   rho_res.resize( lem.cells.size(), 0.0);
   p_res.resize( lem.cells.size(), 0.0);
   Y_res.resize( lem.cells.size());
   for (it = c.begin(); it != c.end(); it++)
   {
    Y_res[std::distance(c.begin(), it)].resize( (*it).ns, 0.0);
   }

   is >> length_res;
   is >> Energy_res;
   is >> time_res;

  for (it = c.begin(); it != c.end(); it++)
  {
    is >> dx_res[std::distance(c.begin(), it)] ;
  }

  for (it = c.begin(); it != c.end(); it++)
  {
    is >> T_res[std::distance(c.begin(), it)] ;
  }

  for (it = c.begin(); it != c.end(); it++)
  {
    is >> m_res[std::distance(c.begin(), it)] ;
  }

  for (it = c.begin(); it != c.end(); it++)
  {
    is >> E_res[std::distance(c.begin(), it)] ;
  }

  for (it = c.begin(); it != c.end(); it++)
  {
    is >> rho_res[std::distance(c.begin(), it)] ;
  }

  for (it = c.begin(); it != c.end(); it++)
  {
    is >> p_res[std::distance(c.begin(), it)] ;
  }

  for (it = c.begin(); it != c.end(); it++)
  {
    for (int s = 0; s < (*it).ns; s++)
    {
      is >> Y_res[std::distance(c.begin(), it)][s];
    }
  }

  lem.initializeRestartData(length_res, Energy_res, time_res, dx_res, T_res,
          m_res, E_res, rho_res, p_res, Y_res);

  return is;

}

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
void LEMLINE::initializeRestartData
        (
                double              length_restart      ,
                double              Energy_restart      ,
                double              time_restart        ,
                std::vector<double> dx_restart          ,
                std::vector<double> T_restart           ,
                std::vector<double> m_restart           ,
                std::vector<double> E_restart           ,
                std::vector<double> rho_restart         ,
                std::vector<double> p_restart           ,
                std::vector<std::vector<double>> Y_restart
        )
{
    std::list<CellData>::iterator it;

    length = length_restart;
    Energy = Energy_restart;
    time = time_restart;
//  cout << "length is " << length << " , energy is " << Energy
//       << " while internalEnergyOnLine() is "
//       << internalEnergyOnLine()  << std::endl;

    for (it = cells.begin(); it != cells.end(); ++it)
    {
      (*it).dx  = dx_restart[std::distance(cells.begin(), it)];
      (*it).T   = T_restart[std::distance(cells.begin(), it)];
      (*it).m   = m_restart[std::distance(cells.begin(), it)];
      (*it).E   = E_restart[std::distance(cells.begin(), it)];
      (*it).rho = rho_restart[std::distance(cells.begin(), it)];
      (*it).p   = p_restart[std::distance(cells.begin(), it)];

      for (int s = 0; s < (*it).ns; s++)
      {
	(*it).Y[s] = Y_restart[std::distance(cells.begin(), it)][s];
      }
    }

  //for (it = cells.begin(); it != cells.end(); ++it)
  //{
  //  if (abs((*it).E - E_restart[std::distance(cells.begin(), it)]) > 1.0)
  //  cout << "Warning: In restart, Energy on lem cell number"
  //      << std::distance(cells.begin(), it)
  //      <<  "is not same as before" << std::endl;
  //}
}
// }}}
