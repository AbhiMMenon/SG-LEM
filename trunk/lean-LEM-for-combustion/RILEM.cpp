#include "LEM.h"
#include "LEMLine.h"
#include "CellData.h"
#include "cantera/equil/ChemEquil.h"

#define MEANBINS 200
#define VARBINS  200



/**
 * @brief       Returns the dy/dx of vector y using a central differencing for
 * interior elements, forward and backward differencing for the ends.
 *
 * @param       \input  y
 * @param       \input  dx      spacing vector
 *
 * @author      A.Menon
 */

std::vector<double> gradient(std::vector<double> y, std::vector<double> dx){

    int nel = y.size();
    std::vector<double> dydx(nel, 0.0);

    // central differencing for interior
    for(size_t ix = 1; ix < nel-1; ix++){
        double dxTot = 0.5*(dx[ix-1] + dx[ix+1]) + dx[ix];
        dydx[ix] = (y[ix+1] - y[ix -1])/(dxTot);
    }

    // forward differencing for near end
    dydx[0] = (y[1] - y[0])/(dx[0]);
    // backward differencing for far end
    dydx[nel -1] = (y[nel -2] - y[nel -1])/(dx[nel-1]);

    return dydx;
}


/**
 * Heavyside fuction
 *
 * EDIT: Now is Alan-step
 *
 * @author      A.Menon
 */


std::vector<double> Heavyside(double mean,size_t size){

    std::vector<double> fx(size, 0.0);

    mean = mean == 0.0? 1e-7: mean;
    mean = mean == 1.0? 1-1e-7: mean;

    double ss = 1./(size-1);
    // -- generate the C PDF
    for(size_t ix = 0; ix < size; ix++){
        double C;
		C = ix* ss;

        // integral areas like the incomplete beta function
        if(C < mean)  fx[ix] = (1.-mean)/mean*C;
        else fx[ix] = (1-mean)+mean/(1.-mean)*(C-mean);
    }
    return fx;
}

/**
 * Top-hat filter based on J.Floyd 2008, computes FDF and the CDF for return
 *
 * @param \input mean
 * @param \input var
 * @param \input number of elements
 *
 * @author      A.Menon
 */

std::vector<double> TopHat(double mean,double var, double ss, size_t size){

    var = var < 5e-6? 5e-6:var; // necessary condition
    std::vector<double> fx(size, 1.0);
    double l  = std::sqrt(12.*var);
    double fa = std::max(0.0, (mean - l/2.));
    double fb = std::min(1.0, (mean + l/2.));

    for(size_t ix = 0; ix < size; ix++){
        double C;
        C = ix* ss;
        if(C < fa) fx[ix] = 0.0;
        else if(C > fb) fx[ix] = 1.0;
        else fx[ix] = (C-fa)/(fb-fa);
    }

    return fx;
}

/**
 * Top-hat filter based on J.Floyd 2008, for bigger variances
 *
 * @param \input mean
 * @param \input var
 * @param \input number of elements
 *
 * @author      A.Menon
 */

std::vector<double> TopHat2(double mean,double var, double ss, size_t size){

    var = std::min(mean-mean*mean, var);
    var = std::max(var,1e-9);
    std::vector<double> fx(size+2, 0.0);
    double l  = std::sqrt(12.*var);
    double fa = mean - l/2;
    double fb = mean + l/2;
    double w0=0, w1=0;

    double t;
    double varMin  = 1/3.*std::min(mean*mean, (1.-mean)*(1. - mean));

    // condition 0 -- Delta peak, handled at solver level

    // condition 1 -- Small variances
    if(0. < var && var <= varMin){
        t = 1./l;
    }
    // condition 2 -- Larger variance, weighted delta at 0
    else if (1/3.*mean*mean < var && var <= mean*(2./3-mean)){
        w0 = 1. - 4./3*mean*mean/(var + mean*mean);
        l = 2*mean/(1. - w0);
        t = (1.-w0)/l;
        fa=0; fb = l;
    }
    // condition 3 -- Larger variance, weighted delta at 1
    else if (1/3.*(1.-mean)*(1.-mean) < var && var <= (1-mean)*(2./3-(1.-mean))){
        w1 = 1. - 4./3*(1.-mean)*(1.-mean)/(var + (1.-mean)*(1.-mean));
        l = 2.*(1.-mean)/(1. - w1);
        t = (1.-w1)/l;
        fa=1.-l; fb = 1.;
    }
    // condition 4 -- very large variance, Bray-Moss-Libby bi-modal PDF

    else{
        fa = 0.0; fb =1.0;
        t = 6.0*(mean- mean*mean-var);
        w1 = var - 1./3.*t + mean*mean;
        w0 = 1.- t -w1;
    }

    for(size_t ix = 0; ix < size; ix++){
        double C;
        C = (ix)* ss;
        if(C < fa) fx[ix+2] = 0.0 + w0;
        else if(C > fb) fx[ix+2] = 1.0 -w1;
        else fx[ix+2] = (C-fa)*t + w0;
    }

    fx[0] = w0;  // needs to be added separately
    fx[1] = w1;


    return fx;
}


/**---------------------------------------------------------------------------
 * @brief               Allocates memory for the matrices and vectors, Z space
 * then C space and ZC space in the last
 *
 * @author              A.Menon 14th Jan 2021
 *---------------------------------------------------------------------------*/
void LEMLINE::initZCMemory()
{


    initZMemory();

    // --- allocate for scalarZC and scalarZC_flag and initialize flags
    scalarZC_flag.resize(boost::extents[ZBinSize][CBinSize]);
    dCdt_flag.resize(boost::extents[ZBinSize][CBinSize]);
    scalarZC.resize(boost::extents[ZBinSize][CBinSize][nSpecies + NSCALAR]);

    for(size_t ix = 0; ix < ZBinSize; ix++){
        for(size_t jx = 0; jx < CBinSize; jx ++){
            dCdt_flag[ix][jx]=false;
            scalarZC_flag[ix][jx]=false;
        }
    }

    ySol.resize(lem.combustor.neq());
}

/**---------------------------------------------------------------------------
 * @author              A.Menon 1/1/2024
 *---------------------------------------------------------------------------*/
void LEMLINE::initBetaShapes(size_t bins,double mMax)
{
    lem.BetaShapes.resize(boost::extents[MEANBINS][VARBINS][bins]);

    for(size_t ix = 1; ix < MEANBINS-1; ix++){
        double m = ix*(1./(MEANBINS-1.))*mMax;
        double vMax = m*(1.-m);
        for(size_t jx = 0; jx < VARBINS; jx++){
            double v = jx*(1./(VARBINS-1.))*vMax*0.999999;
            v = v < 1e-5? 1e-5:v;
        //  std::cout << "m , v= " <<m <<" " <<v<<"\n";

            double gamma = m * (1. - m)/v - 1.;
            double alpha = gamma * m;
            double beta  = gamma * (1. - m);
            for(size_t kx = 0; kx < bins; kx++){
                double x = kx/(bins-1.);
                lem.BetaShapes[ix][jx][kx] = boost::math::ibeta(alpha, beta, x);
            }
        }
    }

    double m = 0.0001;
    double vMax = m*(1.-m);
    for(size_t jx = 0; jx < VARBINS; jx++){
        double v = jx*(1./(VARBINS-1.))*vMax*0.999999;
        v = v < 1e-5? 1e-5:v;
        double gamma = m * (1. - m)/v - 1.;
        double alpha = gamma * m;
        double beta  = gamma * (1. - m);
        for(size_t kx = 0; kx < bins; kx++){
            double x = kx/(bins-1.);
            lem.BetaShapes[0][jx][kx] = boost::math::ibeta(alpha, beta, x);
        }
    }
    m = 0.9999*mMax;
    vMax = m*(1.-m);
    for(size_t jx = 0; jx < VARBINS; jx++){
        double v = jx*(1./(VARBINS-1.))*vMax*0.999999;
        v = v < 1e-5? 1e-5:v;
        double gamma = m * (1. - m)/v - 1.;
        double alpha = gamma * m;
        double beta  = gamma * (1. - m);
        for(size_t kx = 0; kx < bins; kx++){
            double x = kx/(bins-1.);
            lem.BetaShapes[MEANBINS-1][jx][kx] = boost::math::ibeta(alpha, beta, x);
        }
    }

}

/**---------------------------------------------------------------------------
 * @author              A.Menon 1/1/2024
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::interpBeta(double mean, double var,double mMax){
    size_t bins = lem.BetaShapes[0][0].size();
    std::vector<double> fx(bins, 0.0);

    double MEANstep = 1./(MEANBINS-1)*mMax;

    int Mx = std::floor(mean/MEANstep);
    Mx = std::max(Mx, 0);
    Mx = std::min(Mx, MEANBINS);

    double vMax = mean*(1.-mean);
    double VARstep = 1./(MEANBINS-1)*vMax;
    int Vy = var/VARstep;
    Vy = std::max(Vy, 0);
    Vy = std::min(Vy, VARBINS-1);

    int MxInc = std::min(Mx+1,MEANBINS-1);
    int VyInc = std::min(Vy+1,VARBINS-1);


    double diffM = mean - Mx*MEANstep;
    double diffV = var  - Vy*VARstep;

    for(size_t ix = 0; ix < bins;ix++){
        // interpolate along M for Vy
        double m = (lem.BetaShapes[MxInc][Vy][ix] - lem.BetaShapes[Mx][Vy][ix])/MEANstep;
        double sum1 = lem.BetaShapes[Mx][Vy][ix] + m*diffM;

        // interpolate along M for VyInc
               m = (lem.BetaShapes[MxInc][VyInc][ix] - lem.BetaShapes[Mx][VyInc][ix])/MEANstep;
        double sum2 = lem.BetaShapes[Mx][VyInc][ix] + m*diffM;

        // interpolate along Var
        m = (sum2-sum1)/VARstep;
        fx[ix] = sum1 + m*diffV;
    }

//  std::memcpy(fx.data(), lem.BetaShapes[Mx][Vy].origin(),bins*sizeof(double));
    return fx;

}


/**---------------------------------------------------------------------------
 * @brief               Allocates memory for the matrices and vectors for Z
 *                      space only
 *
 *                      26-nov-2021
 *---------------------------------------------------------------------------*/
void LEMLINE::initZMemory(){

    if(lem.Zonly){
        scalarZ.resize(boost::extents[ZBinSize][nSpecies+NSCALAR]);
    }


    if(lem.TUnbZ.size() ==0){
        lem.YEquilZ.resize(boost::extents[ZBinSize][nSpecies]);
        lem.YUnbZ.resize(boost::extents[ZBinSize][nSpecies]);
        lem.TUnbZ.resize(ZBinSize);
        lem.TEquil_z.resize(ZBinSize);
    }



}

/**---------------------------------------------------------------------------
 * @brief               Allocates memory for the matrices and vectors for C
 *                      space only
 *
 *                      26-nov-2021
 *---------------------------------------------------------------------------*/
void LEMLINE::initCMemory(){

    scalarC.resize(boost::extents[CBinSize][nSpecies + NSCALAR]);
    tempC.resize(boost::extents[CBinSize][2*nSpecies + NSCALAR]);

    scalarC_flag.resize(CBinSize,0.0);
    speciesProdC.resize(boost::extents[CBinSize][nSpecies]);

    ySol.resize(lem.combustor.neq());
    Yu_p.resize(lem.ns,0.0);
    Yb_p.resize(lem.ns,0.0);

}





/**---------------------------------------------------------------------------
 * @brief    Destructor, (NO MORE) memory deletions for several matrices
 *
 * @author   A.Menon
 *---------------------------------------------------------------------------*/
  LEMLINE::~LEMLINE(){

  }


/**----------------------------------------------------------------------------
 *  @brief          Initializes LEM line with fuel in the centre, air in the
 *  	  	    rest.
 *
 *  @author         A.Menon
 *---------------------------------------------------------------------------*/
void LEMLINE::initializeJetTest()
{
    auto gas = lem.gas;
    std::list<CellData>::iterator it;

    float di = 0,
          th = 0.10 * length;

    // --- central 'jet' with fuel, rest with oxidiser
    double Zj = 0.044;
    std::vector<double> Yu(nSpecies, 0);
    std::vector<double> Yb(nSpecies, 0);

    for(size_t jx = 0; jx < nSpecies; jx++) {
        Yu[jx] = Zj*lem.fuelY[jx] + (1.-Zj)*lem.oxY[jx];
    }
    std::vector<double> Yu2 = Yu;

    Yu2[gas->speciesIndex("H")]= 1;
    gas->setState_TPY(lem.fuelStreamT, lem.pressure, Yu.data());
    gas->equilibrate("HP");
    gas->getMassFractions(Yb.data());
    double Tb = 2200;

    for (it = cells.begin(); it != cells.end(); ++it){
        di += (*it).dx;
        if((di >= length/2 - th) && (di <= length/2 + th))
            gas->setState_TPY(lem.fuelStreamT, lem.pressure, Yu.data());

        else if((di >= length/2 - th*1.1) && (di <= length/2 + th*1.1)){
         // gas->setState_TPY(300, lem.pressure,Yu2.data() );
            gas->setState_TPY(Tb, lem.pressure, Yb.data());
        }
        else
            gas->setState_TPY(lem.oxStreamT, lem.pressure, lem.oxY.data());

     (*it).p   = lem.pressure;
     (*it).T   = gas->temperature();
     (*it).rho = gas->density();
     (*it).m   =  (*it).dx*(*it).rho;

 //  (*it).m   =  1e-5;
 //  (*it).dx = (*it).m/(*it).rho;

     gas->getMassFractions((*it).Y);
    }

    setEnergy();

}
/**----------------------------------------------------------------------------
 *  @brief          Initializes LEM line with fuel in the centre, air in the
 *  	  	    rest.
 *
 *  @author         A.Menon
 *---------------------------------------------------------------------------*/
void LEMLINE::initCellByCell(double T_in, double P_in, std::vector<double> Y_in, int incr)
{
    auto gas = lem.gas;
    std::list<CellData>::iterator it;

    double T0=T_in;
    double P0=P_in;

    gas->setState_TPY(T0,P0, Y_in.data());

    it = cells.begin();
    std::advance(it,incr);

    (*it).p   = gas->pressure();
    (*it).T   = gas->temperature();
    (*it).rho = gas->density();
    gas->getMassFractions((*it).Y);

    (*it).m = gas->density() * (*it).dx;
    (*it).E= gas->intEnergy_mass() * (*it).m;
}

/**----------------------------------------------------------------------------
 *  @brief      Initializes mixing layer case, 1/3d fuel, 1/3 gradient, 1/3
 *  oxidiser
 *
 *  @author     A.Menon, 24th Aug
 *---------------------------------------------------------------------------*/
void LEMLINE::initializeMixTest()
{
    auto gas = lem.gas;
    std::list<CellData>::iterator it;


    int nCells = cells.size();
    int mPosStart = nCells/3;
    int mPosStop  = nCells*2/3;
    int counter = 0;
    double Y_fuel[nSpecies];
    double Y_ox[nSpecies];
    double Zstep_loc = 1./(mPosStart - 1);

    // fuel props
    gas->setState_TPX(lem.fuelStreamT, lem.pressure, lem.fuelStream);
    gas->getMassFractions(Y_fuel);

    gas->setState_TPX(lem.oxStreamT, lem.pressure, lem.oxStream);
    gas->getMassFractions(Y_ox);

    for (it = cells.begin(); it != cells.end(); ++it){
        if(counter < mPosStart)
            gas->setState_TPX(lem.fuelStreamT, lem.pressure, lem.fuelStream);
        else if (counter > mPosStop)
            gas->setState_TPX(lem.oxStreamT, lem.pressure, lem.oxStream);
        else{
            double Z = 1. - (counter-mPosStart)*Zstep_loc;
            double T = Z* lem.fuelStreamT +
                   (1.-Z)* lem.oxStreamT;
            double Ys[nSpecies];
            for(size_t ix = 0; ix < nSpecies; ix ++)
                Ys[ix] = Z*Y_fuel[ix] + (1.-Z)*Y_ox[ix];

            gas->setState_TPY(T, lem.pressure, &Ys[0]);

            // ingition test
            if(fabs(Z-0.02) < 2e-2){
             // equilibrate(gas,"HP");
                gas->equilibrate("HP");
                T = gas->temperature();      // set to equilibrium temperature
                gas->setState_TPY(T, lem.pressure, &Ys[0]);
            }
        }
        counter++;


        (*it).p   = gas->pressure();
        (*it).T   = gas->temperature();
        (*it).rho = gas->density();

        gas->getMassFractions((*it).Y);

        (*it).m   =  gas->density() * (*it).dx;
        double E   = (gas->enthalpy_mass() - (*it).p /
                 (*it).rho) * (*it).m;
         // added 16th Aug
         (*it).E   = E;
         Energy   += E;
    }

}

/**----------------------------------------------------------------------------
 *  @brief      Initializes line with only oxidier
 *
 *  @author     A.Menon, 4th Sep 2020
 *---------------------------------------------------------------------------*/
void LEMLINE::initializeOx()
{
    auto gas = lem.gas;
    std::list<CellData>::iterator it;

    gas->setState_TPY(lem.oxStreamT, lem.pressure, "O2:1,N2:3.76");
    for (it = cells.begin(); it != cells.end(); ++it){
        (*it).p   = gas->pressure();
        (*it).T   = gas->temperature();
        (*it).rho = gas->density();
        gas->getMassFractions((*it).Y);

        (*it).m    =  gas->density() * (*it).dx;
        double E   = (gas->enthalpy_mass() - (*it).p /
                     (*it).rho) * (*it).m;
    }

    it = cells.begin();
    setEnergy();

}

/**----------------------------------------------------------------------------
 *  @brief      Initializes line with Z
 *
 *  @author     A.Menon, 10th Oct 2020
 *---------------------------------------------------------------------------*/
void LEMLINE::initializeWithZ(double Z) {
    auto gas = lem.gas;
    std::list<CellData>::iterator it; 
    double Y0[nSpecies];
    double T = Z * lem.fuelStreamT + (1. - Z)* lem.oxStreamT;

    for(size_t nx = 0; nx < nSpecies; nx ++)
        Y0[nx] = Z * lem.fuelY[nx] + (1. - Z)*lem.oxY[nx];

    gas->setState_TPY(T, lem.pressure, Y0);

    for (it = cells.begin(); it != cells.end(); ++it){
        (*it).p   = gas->pressure();
        (*it).T   = gas->temperature();
        (*it).rho = gas->density();
        gas->getMassFractions((*it).Y);

        (*it).m    =  gas->density() * (*it).dx;
        (*it).Z    =  Z;
    }
    setEnergy();
}

/**----------------------------------------------------------------------------
 *  @brief      Initializes line with premixed gas of equivalence ratio phi
 *
 *  @author     A.Menon, 29th Nov 2021
 *---------------------------------------------------------------------------*/
void LEMLINE::initializeWithPhi(double phi) {
    setPremixedData(phi);

    auto gas = lem.gas;
    std::list<CellData>::iterator it;

    gas->setState_TPY(TUnb_p, lem.pressure, Yu_p.data());

    for (it = cells.begin(); it != cells.end(); ++it){
        (*it).p   = gas->pressure();
        (*it).T   = gas->temperature();
        (*it).rho = gas->density();
        gas->getMassFractions((*it).Y);

        (*it).m    =  gas->density() * (*it).dx;
     // (*it).dx   =  (*it).m/ gas->density();
    }
    setEnergy();
}

/**----------------------------------------------------------------------------
 *  @brief      Initializes line with non-premixed (alternate fuel and oxidiser)
 *              so as to maintain the average Z
 *
 *  @author     A.Menon, 10th Oct 2020
 *
 *  EDIT : this is completely un-necessary if the per-patch inlet line splicing
 *  in being used.
 *---------------------------------------------------------------------------*/
void LEMLINE::initializeWithZ_staggered(double Z) {
    auto gas = lem.gas;
    std::list<CellData>::iterator it;

    if(Z < 1e-3 ||Z > 1-1e-9){
        initializeWithZ(Z);
        return;
    }

    double dx2 = length/cells.size()*2;
    double rhoAv = Z*lem.fuelRho + (1.-Z)*lem.oxRho;
    double mF = Z*lem.fuelRho/rhoAv;
    double mO = (1-Z)*lem.oxRho/rhoAv;
    gas->setState_TPY(lem.fuelStreamT, lem.pressure, lem.fuelY.data());

    // -- set fuel cells
    for (it = cells.begin(); it != cells.end(); ++it,++it){
        (*it).p   = gas->pressure();
        (*it).T   = gas->temperature();
        (*it).rho = gas->density();
        gas->getMassFractions((*it).Y);

        // -- set the cell size
        (*it).dx  =  mF*dx2;
        (*it).m    =  gas->density() * (*it).dx;
        double E   =  gas->intEnergy_mass()*(*it).m;
        (*it).E   = E;
        Energy   += E;
    }

    // -- set Ox cells
    it = cells.begin(); it++;
    gas->setState_TPY(lem.oxStreamT, lem.pressure, lem.oxY.data());
    for (; it != cells.end(); ++it ,++it){
        (*it).p   = gas->pressure();
        (*it).T   = gas->temperature();
        (*it).rho = gas->density();
        gas->getMassFractions((*it).Y);

        // -- set the cell size
        (*it).dx  =  mO*dx2;
        (*it).m    =  gas->density() * (*it).dx;
        double E   =  gas->intEnergy_mass()*(*it).m;
        (*it).E   = E;
        Energy   += E;
    }
}

/**----------------------------------------------------------------------------
 * @brief  Ignites half the line, premixed only, intended for pilot injection
 *
 * EDIT -- The whole line
 *---------------------------------------------------------------------------*/
void LEMLINE::sparkPremixed(){
    std::list<CellData>::iterator it = cells.begin();
    std::list<CellData>::iterator itMid = std::next(cells.begin(), cells.size()/3*2.);

    auto gas = lem.gas;
    gas->setState_TPY(Tb_p, lem.pressure, Yb_p.data());
    for (it=itMid; it != cells.end(); ++it){
        (*it).T = Tb_p;
        gas->getMassFractions((*it).Y);
        (*it).p   = gas->pressure();
        (*it).rho = gas->density();
        (*it).dx =(*it).m/gas->density();
    }
    setEnergy();

}

/**----------------------------------------------------------------------------
 * @brief  Ignites half the line
 *---------------------------------------------------------------------------*/
void LEMLINE::sparkLine(double startFrac){
    std::list<CellData>::iterator it = cells.begin();
    auto gas = lem.gas;
    std::list<CellData>::iterator itMid = std::next(cells.begin(), cells.size()*startFrac);
    for (it=itMid; it != cells.end(); ++it){
        gas->setState_TPY((*it).T, lem.pressure, (*it).Y);
        gas->equilibrate("HP");
        (*it).T = gas->temperature();
        gas->getMassFractions((*it).Y);
        (*it).p   = gas->pressure();
        (*it).rho = gas->density();
        (*it).dx =(*it).m/gas->density(); //induce dilatation
     // (*it).m =(*it).dx*gas->density();
    }
    setEnergy();
}

/**----------------------------------------------------------------------------
 * @brief  boost combustion for initial time steps
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::injectRadicals(){
    auto gas = lem.gas;
    std::list<CellData>::iterator it = cells.end();
    std::advance(it,-10);

    int indH = gas->speciesIndex("H");

    for (it; it != cells.end(); ++it){
        (*it).Y[indH] = 0.01;
        gas->setState_TPY((*it).T, (*it).p, (*it).Y);
        gas->getMassFractions((*it).Y);// normalize
        (*it).rho = gas->density();
        (*it).dx = (*it).m/(*it).rho;
    }

    setEnergy();

}

/**----------------------------------------------------------------------------
 * @brief  Ignites whole line, intended for pilot injection
 *---------------------------------------------------------------------------*/
void LEMLINE::sparkNonPremixed(){
    std::list<CellData>::iterator it = cells.begin();

    for (it = cells.begin(); it != cells.end(); ++it){
        (*it).rho = (*it).rho*3000/(*it).T;
        (*it).T = 3000;
        (*it).dx =(*it).m/(*it).rho;
    }
    setEnergy();

}


/**----------------------------------------------------------------------------
 *  @brief      Re-Initializes mixing layer boundary, 1/3d fuel,1/3
 *
 *  @author     A.Menon, 1st Sep 2020
 *---------------------------------------------------------------------------*/
void LEMLINE::boundaryMixTest()
{
    auto gas = lem.gas;
    std::list<CellData>::iterator it;


    // fuel props

    double Y_fuel[nSpecies];
    double Y_ox[nSpecies];
    gas->setState_TPX(lem.fuelStreamT, lem.pressure, lem.fuelStream);
    gas->getMassFractions(Y_fuel);

    gas->setState_TPX(lem.oxStreamT, lem.pressure, lem.oxStream);
    gas->getMassFractions(Y_ox);

    int Zx = 0;
    double gradZ = 1./cells.size();
    for (it = cells.begin(); it != cells.end(); ++it){
        double Z = Zx*gradZ;
        Zx++;
        double Ys[nSpecies];
        for(size_t ix = 0; ix < nSpecies; ix ++){
            Ys[ix] = Z*Y_fuel[ix] + (1.-Z)*Y_ox[ix];
        }

        double T = Z* lem.fuelStreamT +
               (1.-Z)* lem.oxStreamT;
        gas->setState_TPY(T, lem.pressure, &Ys[0]);
        (*it).p   = gas->pressure();
        (*it).T   = gas->temperature();
        (*it).rho = gas->density();

        gas->getMassFractions((*it).Y);
        (*it).m   =  gas->density() * (*it).dx;
    }

    setEnergy();

}

/**----------------------------------------------------------------------------
 *  @brief          Computes fuel and oxidiser stream properties, i.e
 *                  gammas (from gammaCoeff), beta1 and beta 2 used in
 *                  Bilger's formula.
 *
 *  @author         A.Menon
 ----------------------------------------------------------------------------*/
void LEMLINE::setStreams(){
    auto gas = lem.gas;
    double Y[nSpecies];

    std::vector<double> elemMassFrac;  // elemental mass fractions

    // set gammas using gammaCoeff
    lem.gammas[0] = lem.gammasCoeff[0]/gas->atomicWeight(gas->elementIndex("C"));
    lem.gammas[1] = lem.gammasCoeff[1]/gas->atomicWeight(gas->elementIndex("H"));
    lem.gammas[2] = lem.gammasCoeff[2]/gas->atomicWeight(gas->elementIndex("O"));
    lem.gammas[3] = lem.gammasCoeff[3]/gas->atomicWeight(gas->elementIndex("N"));

    // --- Set the fuel stream mass fractions

    gas->setState_TPX(300, Cantera::OneAtm,lem.fuelStream);
    gas->getMassFractions(Y);
    gas->getMassFractions(lem.fuelY.data());
    lem.fuelRho=gas->density();
    elemMassFrac = getElementMassFrac(Y);

    // set coefficients for Bilger's formula
    // --- Fuel, beta 1

    lem.beta1 = lem.gammas[0] * elemMassFrac[gas->elementIndex("C")] +
                lem.gammas[1] * elemMassFrac[gas->elementIndex("H")] +
                lem.gammas[2] * elemMassFrac[gas->elementIndex("O")] +
                lem.gammas[3] * elemMassFrac[gas->elementIndex("N")];

    // --- Set the oxidiser stream mass fractions

    gas->setState_TPX(300, Cantera::OneAtm,lem.oxStream);
    gas->getMassFractions(Y);
    gas->getMassFractions(lem.oxY.data());
    lem.oxRho=gas->density();
    elemMassFrac = getElementMassFrac(Y);

    // set coefficients for Bilger's formula
    // --- oxidiser, beta2

    lem.beta2 = lem.gammas[0] * elemMassFrac[gas->elementIndex("C")] +
                lem.gammas[1] * elemMassFrac[gas->elementIndex("H")] +
                lem.gammas[2] * elemMassFrac[gas->elementIndex("O")] +
                lem.gammas[3] * elemMassFrac[gas->elementIndex("N")];

}


/* -----------------------------------------------------------------------------
 * sets the data for TUnb_p, Tb_p, Yu_p, Yb_p  for a premixed cases based on
 * equivalence ratio, phi
 *
 * m and n give the hydrocarbon molecule based on which the stoichiometric O2
 * is calculated
 *
 * CmHn + (m+n/4)O2 -> some CO2 and H2O
 *
 *
 * 29-Nov-2021
 *------------------------------------------------------------------------------*/

void LEMLINE::setPremixedData(double phi){

    auto gas = lem.gas;

    std::string fuel = lem.fuelStream;
    size_t sepPos = fuel.find(':');
    fuel.erase(fuel.begin()+sepPos, fuel.end());

    // -- find m and n
    gas->setState_TPX(300, 101325, lem.fuelStream);
    int indC = gas->elementIndex("C");
    int indH = gas->elementIndex("H");
    int indF  = gas->speciesIndex(fuel);
    int indO2 = gas->speciesIndex("O2");
    int indN2 = gas->speciesIndex("N2");

    int m = gas->nAtoms(indF, indC);
    int n = gas->nAtoms(indF, indH);

    if(m+n==0){
        std::cout << "ERROR setPremixedData || check fuel composition"<<
         ", only valid for single component hydrocarbon to burn in air, exiting";
        exit(0);

    }

    // -- find oxidizer composition

    gas->setState_TPX(300, 101325, lem.oxStream);

    double X[lem.ns]{};
    gas->getMoleFractions(X); // oxidiser mole fractions

    double Nrat = X[indN2]/X[indO2]; // nitrogen ratio


    // -- set  mole fractions for phi
    X[indF] = 1.0;
    X[indO2] = (m+n/4.) * (1./phi) ;
    X[indN2] = (m+n/4.) * (1./phi) * Nrat ;

    // -- unburnt premixed mass fractions
    gas->setState_TPX(lem.fuelStreamT, lem.pressure, X);
    gas->getMassFractions(Yu_p.data());
    TUnb_p = lem.oxStreamT;


    // -- burnt mass fractions and temperature
    gas->equilibrate("HP");
    gas->getMassFractions(Yb_p.data());
    Tb_p = gas->temperature();

}


/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advanceent
 * @param               diffFlag       /input  diffusion flag
 * @param               combFlag       /input  comustion flag
 *
 * @brief               advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source.
 *
 * @author              M Oevermann; modifications by S.Arshad & A.Menon
 ----------------------------------------------------------------------*/
//void LEMLINE::advanceLEM(const double &dt, double threshold)
//{
//  Total_EddiesLength    = 0.0;
//  Average_EddiesLength  = 0.0;
//
//  std::list<CellData>::iterator it(cells.begin());
//
//  dtLEM                 = dt / 10.;
//  const int nDt         = ceil(dt / dtStirr * 2.);
//  const double Dt       = dt / nDt;
//
//  for (int n = 0; n < nDt; n++){
//      if (time + Dt > timeNextTripletMap){
//          double dt1        = timeNextTripletMap - time;
//          double dt2        = time + Dt - timeNextTripletMap;
//
//          if (dt1 > threshold){
//
//              if(diffFlag)
//                advanceDiffusion(0.5 * dt1);
//              if(combFlag)
//                advanceSourceTerm(dt1);
//              if(diffFlag)
//                advanceDiffusion(0.5 * dt1);
//
//              setBilgerMixFrac();
//              setProgressVariable_T();
//              setEnergy();
//              calcScalarZC();
//          }
//          eddyEvent();
//       // timeNextTripletMap += dtStirr;
//          sampleEddyTime();
//
//          if (dt2 > threshold){
//              if(diffFlag)
//                advanceDiffusion(0.5 * dt2);
//              if(combFlag)
//                advanceSourceTerm(dt2);
//              if(diffFlag)
//                advanceDiffusion(0.5 * dt2);
//
//              setBilgerMixFrac();
//              setProgressVariable_T();
//              setEnergy();
//              calcScalarZC();
//          }
//      }
//
//      else{
//          if(diffFlag)
//            advanceDiffusion(0.5 * Dt);
//          if(combFlag)
//            advanceSourceTerm(Dt);
//          if(diffFlag)
//            advanceDiffusion(0.5 * Dt);
//
//          setBilgerMixFrac();
//          setProgressVariable_T();
//          setEnergy();
//          calcScalarZC();
//    }
//    time += Dt;
//  }
//
//}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               diffFlag       /input  diffusion flag
 * @param               combFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source. New method that should be able to
 *                      combine the advantages of brut clustering, time-based
 *                      clustering (with eddy sampling)
 *
 * @author              A.Menon 19th Jul 2021, modified 3Aug
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEM2(const double &dt, double threshold, bool fastOn, bool chemFlag)
{

    // Check for small time steps and high turbulence, if yes cluster an
    // average number of eddies (dt/dtStirr), else perform the usual sampling
    // with a time based threshold for the eddy clustering

    if(dt <= 1e-5 && (dtStirr < dt) && fastOn){
        size_t nEddy = dt/dtStirr;

        // Just perform the average number of eddies at the beginning
        for(size_t ix = 0; ix < nEddy; ix++) eddyEvent();

        // advance diffusion + reaction for the remainder
        advanceDiffusion(dt/2);
        if(chemFlag)advanceSourceTerm(dt);
        advanceDiffusion(dt/2);

        // diagnose and condition
        setBilgerMixFrac();
        setProgressVariable_T();
        setEnergy();
        calcScalarZC();
        return;
    }

    // Else..

    double elapsedTime = 0;
    double dt1 = sampleEddyTime(); //inter-eddy gap
    double accumTime = 0;

    while(elapsedTime < dt){
        // advance to earlier of eddyTime or end of timeStep
        double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);
        accumTime += diffStep;

        // if the accumulated eddy gaps satisfy the time threshold, advance the
        // diffusion-reaction for the accumulated time. This means that eddies
        // are implemented sequentially (diff-chem skip) if the eddy gaps are
        // too small.

        if(accumTime >= threshold){
            advanceDiffusion(accumTime/2);
            if(chemFlag){advanceSourceTerm(accumTime);}
            advanceDiffusion(accumTime/2);

            // now all the skipped diffusions so far have been accounted for
            accumTime = 0;

            // diagnose and condition
            setBilgerMixFrac();
            setProgressVariable_T();
            setEnergy();
            calcScalarZC();

        }

        elapsedTime  += dt1;
        eddyEvent();
        dt1 = sampleEddyTime(); // sample next eddy
    }
}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source term. Testing only fast advance
 *
 * @author              A.Menon 17th Aug 2021, modified
 * EDIT:                added bFac 25th Aug
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEM3(const double &dt, bool chemFlag)
{
     double elapsedTime = 0;
     double dt1 = sampleEddyTime(); //inter-eddy gap


     while(elapsedTime < dt){
         double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);
         elapsedTime  += dt1;
         if(dt1 < dt){
             eddyEvent();
             dt1 = sampleEddyTime(); // sample next eddy
         }

     }

    // advance diffusion + reaction for the remainder
    advanceDiffusion(dt/2);
    if(chemFlag) advanceSourceTerm(dt);
    advanceDiffusion(dt/2);

    // diagnose and condition
    setEnergy();
    setBilgerMixFrac();
    setProgVble_nonpremixSp();
    if(scalarZC.size()!=0)
        calcScalarZC();
    if(lem.Zonly) calcScalarZ();

}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source term. Testing only fast advance
 *
 * @author              A.Menon 17th Aug 2021, modified
 * EDIT:                added bFac 25th Aug
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEM3DD(const double &dt, bool chemFlag)
{

     double elapsedTime = 0;
     double dt1 = sampleEddyTime(); //inter-eddy gap

     while(elapsedTime < dt){
         double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);
         elapsedTime  += dt1;
         if(dt1 < dt){
             eddyEvent();
             dt1 = sampleEddyTime(); // sample next eddy
         }

     }

    // advance diffusion + reaction for the remainder
    advanceDiffDiffusion(dt/2);
    if(chemFlag) advanceSourceTerm(dt);
    advanceDiffDiffusion(dt/2);

    // diagnose and condition
    setEnergy();
    setBilgerMixFrac();
    setProgVble_nonpremixSp();
    calcScalarZC();

}



/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source term. Testing only
 *                      diffuse-eddy-diffuse
 *
 * @author              A.Menon 23th Aug 2021, modified
 *
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEM4(const double &dt, bool chemFlag)
{

    double elapsedTime = 0;
    double dt1 = sampleEddyTime(); //inter-eddy gap


    while(elapsedTime < dt){
        double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);

        // diffuse-eddy-diffuse
        advanceDiffusion(diffStep/2);
        if(chemFlag && (true)){advanceSourceTerm(diffStep);}
        advanceDiffusion(diffStep/2);

        // diagnose and condition
        setEnergy();
        setHeatRelease();
        setBilgerMixFrac();
        setProgVble_nonpremixSp();
        if(lem.Zonly) calcScalarZ();
        calcScalarZC();


        // --
        elapsedTime  += dt1;
        if(dt1 < dt){
            eddyEvent();
            dt1 = sampleEddyTime(); // sample next eddy
        }

    }

}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source term.
 *
 *                      Operator splitting: eddy-diffuse-eddy for LEM (no
 *                      chem). Advance chemistry at the end. For use with
 *                      heavier mechanisms.
 *
 * @author              A.Menon 23rd Nov 2023
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEM5(const double &dt, bool chemFlag)
{
    double elapsedTime = 0;
    double dt1 = sampleEddyTime(); //inter-eddy gap


    while(elapsedTime < dt){
        double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);
        // diffuse-eddy-diffuse
        advanceDiffusion(diffStep);
        // --
        elapsedTime  += dt1;
        if(dt1 < dt){
            eddyEvent();
            dt1 = sampleEddyTime(); // sample next eddy
        }
    }
    if(chemFlag ){advanceSourceTerm(dt);}

    // diagnose and condition
    setEnergy();
    setHeatRelease();
    setBilgerMixFrac();
    setProgVble_nonpremixSp();
    if(lem.Zonly) calcScalarZ();
    calcScalarZC();

    if(lem.Zonly) calcScalarZ();

}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source term. Testing only
 *                      diffuse-eddy-diffuse
 *
 * @author              A.Menon 23th Aug 2021, modified
 *
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEM1(const double &dt, bool chemFlag)
{
    int nEddy = std::round(dt/dtStirr);

    if (nEddy ==0){
        advanceDiffusion(dt/2);
        advanceSourceTerm(dt);
        advanceDiffusion(dt/2);
        return;

    }

    double diffStep = dt/(nEddy);

    for(size_t ix = 0; ix < nEddy; ix++){
        advanceDiffusion(diffStep/2);
        advanceSourceTerm(diffStep);
        advanceDiffusion(diffStep/2);
        eddyEvent();
        advanceDiffusion(diffStep/2);
        advanceSourceTerm(diffStep);
        advanceDiffusion(diffStep/2);
    }

}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffuse, no source term.
 *                      Uniform spacing of eddy-diffuse-eddy.
 *
 * @author              A.Menon 23th Aug 2021, modified
 *
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEM0(const double &dt)
{

  //int nEddy = floor(dt/dtStirr);
  //if (nEddy ==0){
  //    advanceDiffusion(dt);
  //    return;

  //}

  //double diffStep = dt/(nEddy);

  //for(size_t ix = 0; ix < nEddy; ix++){
  //    advanceDiffusion(diffStep/2);
  //    advanceSourceTerm(diffStep);
  //    advanceDiffusion(diffStep/2);
  //    eddyEvent();
  //    advanceDiffusion(diffStep/2);
  //    advanceSourceTerm(diffStep);
  //    advanceDiffusion(diffStep/2);
  //}

    double elapsedTime = 0;
    double dt1 = sampleEddyTime(); //inter-eddy gap


    while(elapsedTime < dt){
        double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);
        // diffuse-eddy-diffuse
        advanceDiffusion(diffStep);
        // --
        elapsedTime  += dt1;
        if(dt1 < dt){
            eddyEvent();
            dt1 = sampleEddyTime(); // sample next eddy
        }
    }

}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source term. Only diffuse-eddy-diffuse for
 *                      pemixed flames
 *
 *                      1. Sample TE1
 *                      2. Advance to min(TCFD, TE)
 *                      3. if (t== TE) perfom TripletMap
 *                      4. sample TE2
 *                              if TE2 < TE1
 *                                 TE1 = TE2
 *                      5. Repeat
 *
 * @author              A.Menon 29th Aug, 2021
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEMPremixed(const double &dt, bool chemFlag)
{

     double elapsedTime = 0;
     double dt1 = sampleEddyTime(); //inter-eddy gap

     while(elapsedTime < dt){
         double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);

         advanceDiffusion(diffStep/2);
         if(chemFlag){advanceSourceTerm(diffStep);}
         advanceDiffusion(diffStep/2);

         // diagnose and condition
         setEnergy();
         setHeatRelease();
         setProgVble_premixSp();
         calcScalarC();
         calcSpeciesC();

         // --
         elapsedTime  += dt1;
         if(dt1 < dt){
             eddyEvent();
             dt1 = sampleEddyTime(); // sample next eddy
         }

     }

     setEnergy();
     setHeatRelease();
     setProgVble_premixSp();
     calcScalarC();
     calcSpeciesC();

}


/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @brief               Advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source term. Only eddy culstering for  pemixed flames
 *
 * @author              A.Menon 29th Aug, 2021
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEMFastPremixed(const double &dt, bool chemFlag)
{

     double elapsedTime = 0;
     double dt1 = sampleEddyTime(); //inter-eddy gap
     while(elapsedTime < dt){
         double diffStep = (elapsedTime + dt1)<dt? dt1: (dt-elapsedTime);
         elapsedTime  += dt1;
         if(dt1 < dt){
             eddyEvent();
             dt1 = sampleEddyTime(); // sample next eddy
         }

     }


    advanceDiffusion(dt/2);
    if(chemFlag ){advanceSourceTerm(dt);}
    advanceDiffusion (dt/2);

    // diagnose and condition
    setEnergy();
    setHeatRelease(); // needed?
    setProgVble_premixSp();
    calcSpeciesC();
    calcScalarC();


}

/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advancement
 * @param               chemFlag       /input  comustion flag
 *
 * @author              A.Menon 19th Aug, 2023
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEMFast2(const double &dt, bool chemFlag)
{

    size_t nEddy = ceil(dt/dtStirr);
    // -- perform all the eddies
    double dT = dt;

    for(size_t ix = 0; ix < nEddy; ix++) {eddyEvent();}

    // -- diffuse+react once per substep

    advanceDiffDiffusion(dT/2);
    if(chemFlag ){advanceSourceTerm(dT);}
    advanceDiffDiffusion(dT/2);

    // diagnose and condition
    setEnergy();
    setHeatRelease();
    setBilgerMixFrac();
    setProgVble_nonpremixSp();
    calcScalarZ();
    calcScalarZC();


}


/**--------------------------------------------------------------------
 * @param               dt             /input  time step for advanceent
 * @param               combFlag       /input  comustion flag
 *
 * @brief               advances the time on the LEM line, performs triplet
 *                      mapping, thermal and Fickian diffusion, advances
 *                      chemistry source. Maps are done in bunches which fall
 *                      within a CFL diffusion criteria. Variation of operator
 *                      splitting from Sankaran 2003 (Phd. Thesis)
 *
 *                      UPDATED: 17th May
 *
 * @author              A.Menon, 22nd Dec 2020
 ----------------------------------------------------------------------*/
void LEMLINE::advanceLEMFast(const double &dt, double &threshold){

    double elapsedTime  = 0.0;
    double tMap         = 0.0;
    double tDiff        = 0.0;


/*
    while(elapsedTime < dt){

        //Eddy clustering with threshold check
        for (int n = 0; n < nDt; n++){
            // perform several eddies in sequence
            eddyEvent();
        }
        advanceDiffusion(dtDiffStep/2.);
  //    advanceSourceTerm(dtDiffStep);
        advanceDiffusion(dtDiffStep/2.);

        elapsedTime += dtDiffStep;

        setBilgerMixFrac();
        setProgressVariable_T();
        calcScalarZC();
    }

    */

    while(elapsedTime < dt){

        if(tDiff < tMap){
            advanceDiffusion(threshold/2.);
            advanceSourceTerm(threshold);
            advanceDiffusion(threshold/2.);

            elapsedTime += threshold;
            tDiff       += threshold;

            setBilgerMixFrac();
            setProgressVariable_T();
            setEnergy();
            calcScalarZC();
        }
        else{
            tMap += sampleEddyTime(); //returns eddyGap
            eddyEvent();
        }

    }
}


/**----------------------------------------------------------------------------
 * @brief       Returns the elemental mass fractions for a given mass fraction
 *              vector
 *              Z_j = sigma(1,k, (a_ij W_j / W_i * Y_i))
 *              from 'Combustion theory and application'- Prof. Heinz Pitsch
 *
 * @param       \input          Y               mass fraction vector
 * @return      \output         elemMassFrac    std::vector<double>
 * @author      A.Menon
 *
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getElementMassFrac(double *Y)
{
    auto gas = lem.gas;
    std::vector<double> &elemMassFrac= lem.tmp_s3;
    std::vector<double> &atomArray   = lem.tmp_s4;
    std::fill(elemMassFrac.begin(), elemMassFrac.end(), 0.0);
    double sum;

    gas->setMassFractions(Y);

    for(size_t ix = 0; ix < lem.nel; ix++){
       double atWt = gas->atomicWeight(ix);
        for(size_t jx = 0; jx < nSpecies; jx++){

            gas->getAtoms(jx, &atomArray[0]); // atoms for species jx

            elemMassFrac[ix] += atomArray[ix] * atWt /
                                   lem.Ms[jx] * Y[jx];
        }
    }
    return elemMassFrac;
}


/**----------------------------------------------------------------------------
 * @brief      returns mixture fraction using Bilger's Formula.
 *
 * @param       \input  y       species mass fractions (double*)
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/
double LEMLINE::getBilgerMixFrac(double *y){

    if(fabs(lem.beta1 - lem.beta2) < 1E-9)
        return 0.0;

    std::vector<double> &el1 = lem.tmp_s4;
    el1 =   getElementMassFrac(y);
    auto gas =  lem.gas;
    double beta;

    beta = lem.gammas[0] * el1[gas->elementIndex("C")] +
           lem.gammas[1] * el1[gas->elementIndex("H")] +
           lem.gammas[2] * el1[gas->elementIndex("O")] +
           lem.gammas[3] * el1[gas->elementIndex("N")];


    double mixF = (beta - lem.beta2)/(lem.beta1 - lem.beta2);
    if(mixF < 0.0) mixF = 0.0;  // Just to be sure
    if(mixF > 1.0) mixF = 1.0;

    return mixF;
}

/**----------------------------------------------------------------------------
 * @brief      returns mean mixture fraction on Line
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/
double LEMLINE::getMeanMixFrac(){

    double mixF = 0;
    double vol = 0;
    std::list<CellData>::iterator it = cells.begin();
    for(;it != cells.end(); it++){
        double mCell = (*it).dx;
        vol += mCell;
        mixF += (*it).Z*mCell;
     }

    return mixF/vol;
}

/**----------------------------------------------------------------------------
 * @brief      returns max mixture fraction on Line
 *
 * @author      A.Menon 14th Jan
 * EDIT         gets min and max NOT IN USE
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getMinMaxMixFrac()
{
    std::list<CellData>::iterator it;
    std::vector<double> Z_minMax(2,0);

    for (it = cells.begin(); it != cells.end(); ++it){
        Z_minMax[0] = std::min(Z_minMax[0], (*it).Z);
        Z_minMax[1] = std::max(Z_minMax[1], (*it).Z);
    }

    return Z_minMax;
}

/**----------------------------------------------------------------------------
 * @brief      returns full-ness of C space in %age
 *
 * @author      A.Menon 25Mar 2022
 * EDIT         gets min and max
 *---------------------------------------------------------------------------*/
double LEMLINE::getCFullness()
{
    double full = 0;

    for (size_t ix = 0; ix < CBinSize; ix++){
        full += scalarC_flag[ix];
    }

    return full/CBinSize * 100;
}

/**----------------------------------------------------------------------------
 * @brief       Sets mixture fraction for each cell using the Bilger Formula
 *
 *
 * @author      A.Menon
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::setBilgerMixFrac()
{
    std::list<CellData>::iterator it;
    for(it = cells.begin(); it != cells.end(); it++){
        (*it).Z = getBilgerMixFrac((*it).Y);
    }
}

/**----------------------------------------------------------------------------
 * @brief       Sets mixture fraction for each cell using inert specie mass
 * fraction
 *
 * @author      A.Menon 30th Nov
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::setTracerMixFrac(std::string specie)
{
    std::list<CellData>::iterator it;
    auto gas =  lem.gas;

    // check using elemental mass fractions
    size_t pos = gas->elementIndex(specie);
    std::vector<double> elemMassFracFuel = getElementMassFrac(lem.fuelY.data());
    double YC1 = elemMassFracFuel[pos];

    if(specie == "C" || specie == "H" || specie == "Ar"){
        for(it = cells.begin(); it != cells.end(); it++){
            std::vector<double> elemMassFrac = getElementMassFrac((*it).Y);
            double YC = elemMassFrac[pos];
            (*it).Z = YC/YC1 > 1.0 ? 1.0: YC/YC1;
        }
    return;
    }

    // check using molecular mass fractions
    pos = gas->speciesIndex(specie);
    double Y1 = lem.fuelY[pos];

    for(it = cells.begin(); it != cells.end(); it++){
        double  Yinert = (*it).Y[pos];
        double Z = Yinert / Y1;
        Z = Z >1.? 1:Z; // check
        (*it).Z = Z;

    }
}


/**----------------------------------------------------------------------------
 * @brief       Sets the equilibrium values for mass fractions for Zs based on
 *              ZBinSize (in constructor)
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/
void LEMLINE::setEquilTY_Z()
{
    auto gas = lem.gas;
    std::vector<double> Ynow(lem.ns, -1.0);

    for(size_t ix=0; ix<ZBinSize; ix++)
    {
        double Zm = ix*Zstep; // mean value of Z bin

        // --- set mass fractions
        for(size_t jx = 0; jx < lem.ns; jx++){
            Ynow[jx]= Zm*lem.fuelY[jx] + (1-Zm)*lem.oxY[jx];
            lem.YUnbZ[ix][jx] = Ynow[jx];
        }

        // --- equilibrate
        double T = Zm*lem.fuelStreamT + (1.-Zm)*lem.oxStreamT;
        lem.TUnbZ[ix] = T;

        gas->setState_TPY(T,lem.pressure,Ynow.data());
        gas->equilibrate("HP");
        gas->getMassFractions(Ynow.data());

        for(size_t jx=0; jx < lem.ns; jx++){
            lem.YEquilZ[ix][jx] = Ynow[jx];
        }


        lem.TEquil_z[ix] = gas->temperature();

    }

    // find the max temp in the distribution
   // Tstoich  = *(std::max_element(std::begin(TEquil_z), std::end(TEquil_z)));

}

/**----------------------------------------------------------------------------
 * @brief
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/
void LEMLINE::copyMats(LEMLINE &ref)
{

    size_t RSize = nSpecies + NSCALAR;
    for(size_t ix = 0; ix < scalarC.size(); ix++){

        TUnb_p = ref.TUnb_p;
        Tb_p   = ref.Tb_p  ;
        for(size_t kx = 0; kx < RSize;kx++){
            scalarC[ix][kx]  = ref.scalarC[ix][kx];
            scalarC_flag[ix] = ref.scalarC_flag[ix];
        }

        for(size_t kx = 0; kx < lem.ns;kx++){
            speciesProdC[ix][kx]  = ref.speciesProdC[ix][kx];
            Yu_p[kx] = ref.Yu_p[kx];
            Yb_p[kx] = ref.Yb_p[kx];
        }
    }

    Tstoich = ref.Tstoich;

}

/**----------------------------------------------------------------------------
 * @brief
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/

void LEMLINE::copyZCMats(LEMLINE &ref)
{

    size_t RSize = nSpecies + NSCALAR;

    for(size_t ix = 0; ix < ZBinSize; ix++){
      //  TEquil_z[ix] = ref.TEquil_z[ix];

        if(lem.Zonly)
        for(size_t kx = 0; kx < RSize;kx++)
            scalarZ[ix][kx]  = ref.scalarZ[ix][kx];

        for(size_t kx = 0; kx < nSpecies; kx++){
          //YEquilZ[ix][kx] = ref.YEquilZ[ix][kx];
          //  lem.YUnbZ[ix][kx] = ref.YUnbZ  [ix][kx];
        }

        for(size_t jx = 0; jx < CBinSize; jx++){
			scalarZC_flag[ix][jx]  = ref.scalarZC_flag[ix][jx];
			dCdt_flag[ix][jx]      = ref.dCdt_flag[ix][jx];
            for(size_t kx = 0; kx < RSize;kx++)
                scalarZC[ix][jx][kx]  = ref.scalarZC[ix][jx][kx];
        }
    }

}

/**----------------------------------------------------------------------------
 * @brief       Sets the equilibrium values for mass fractions for Zs based on
 *              ZBinSize (in constructor)
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/
void LEMLINE::clearCmats()
{

    size_t RSize = nSpecies + NSCALAR;

 // for(size_t ix = 0; ix < ZBinSize; ix++){

 //     for(size_t kx = 0; kx < RSize;kx++)
 //         scalarZ[ix][kx]  = 0;

 //     for(size_t jx = 0; jx < CBinSize; jx++){
 //         for(size_t kx = 0; kx < RSize;kx++)
 //             scalarZC[ix][jx][kx]  = 0;
 //     }
 // }

    for(size_t ix = 0; ix < CBinSize; ix++){

        for(size_t kx = 0; kx < RSize;kx++){
            scalarC[ix][kx]  = 0;
            scalarC_flag[ix] = 0;
        }

        for(size_t kx = 0; kx < lem.ns;kx++)
            speciesProdC[ix][kx]  = 0;
    }

}


/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarZ[Zbin][scalar] with unbunrt values
 *
 * @author      A.Menon, 18th Nov
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarZ()
{

    auto gas = lem.gas;
    //  ---
    for(size_t ix = 0; ix < ZBinSize; ix++){
        double Z = Zstep*(ix);

        // temperature
        double T = Z*lem.fuelStreamT + (1.-Z)*lem.oxStreamT;

        // Ys
        for(size_t jx = 0; jx < nSpecies; jx++) {
            scalarZ[ix][IY+jx] = Z*lem.fuelY[jx] + (1.-Z)*lem.oxY[jx];
         // scalarZ[ix][IY+jx] = YEquilZ[ix][jx]; // test
        }
        // -- formation enthalpy
        gas->setState_TPY(298.15,101325,scalarZ[ix].origin()+IY);
        double hForm = gas->enthalpy_mass();
        // -- sensible enthalpy
        gas->setState_TPY(T,lem.pressure,scalarZ[ix].origin()+IY);
        double hs = gas->enthalpy_mass() -hForm;

        scalarZ[ix][IT]   = T;
     // scalarZ[ix][IT]   = TEquil_z[ix];
        scalarZ[ix][IRho] = gas->density();
        scalarZ[ix][Ih]   = hs;
        scalarZ[ix][ITheta] = 0;
        scalarZ[ix][IC]     = 0;
    }
}

/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarC[Cbin][scalar] with unbunrt values (C=0)
 *              OBS: only for beta PDF integration and after setPremixedData()
 *
 * @author      A.Menon, 31Jan 2022
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarC()
{

    auto gas = lem.gas;
    //  ---
    for(size_t ix = 0; ix < CBinSize; ix++){

        double T = TUnb_p;
        // Ys

        for(size_t jx = 0; jx < nSpecies; jx++) {
           scalarC[ix][IY+jx] = Yu_p[jx];
        }
        // -- formation enthalpy
        gas->setState_TPY(298.15,101325,scalarC[ix].origin()+IY);
        double hForm = gas->enthalpy_mass();
        // -- sensible enthalpy
        gas->setState_TPY(T,lem.pressure,scalarC[ix].origin()+IY);
        double hs = gas->enthalpy_mass() -hForm;

        scalarC[ix][IT]   = T;
        scalarC[ix][IRho] = gas->density();
        scalarC[ix][Ih]   = hs;
        scalarC[ix][ITheta] = 0;
        scalarC[ix][IC]     = 0;
        scalarC[ix][IC2]    = 0;

        scalarC_flag[ix] = 1;
    }
}

/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarC[Cbin][scalar] with unbunrt values (C=0)
 *              OBS: only for beta PDF integration and after setPremixedData()
 *
 * @author      A.Menon, 31Jan 2022
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarC_linear()
{

    auto gas = lem.gas;
    //  ---
    for(size_t ix = 0; ix < CBinSize; ix++){


        double C = ix * Cstep;
        double T = C*Tb_p + (1. -C) * TUnb_p;

        for(size_t jx = 0; jx < nSpecies; jx++) {
           scalarC[ix][IY+jx] = C*Yb_p[jx] + (1.-C)*Yu_p[jx];
        }
        // -- formation enthalpy
        gas->setState_TPY(298.15,101325,scalarC[ix].origin()+IY);
        double hForm = gas->enthalpy_mass();
        // -- sensible enthalpy
        gas->setState_TPY(T,lem.pressure,scalarC[ix].origin()+IY);
        double hs = gas->enthalpy_mass() -hForm;

        scalarC[ix][IT]   = T;
        scalarC[ix][IRho] = gas->density();
        scalarC[ix][Ih]   = hs;
        scalarC[ix][ITheta] = 0;
        scalarC[ix][IC]     = 0;
        scalarC[ix][IC2]    = 0;

        scalarC_flag[ix] = 1;
    }
}


/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarZ[Zbin][Cbin][scalar] with equilibrium
 * values, simple linear in C for T and Yi; rest are 0
 *
 *
 * @author      A.Menon, 19th Nov
 *
 * UPDATE  used to initialized to zeros
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarZC_zero()
{
    for(size_t ix = 0; ix < ZBinSize; ix++){

        double Z  = Zstep * (ix);
        double Tu = Z * lem.fuelStreamT + (1.- Z)*lem.oxStreamT;
        scalarZC[ix][0][IT] = Tu;

        for(size_t kx = 0; kx < nSpecies + NSCALAR; kx++) {
            scalarZC[ix][0][kx+IY] = lem.YUnbZ[ix][kx];
        }
        scalarZC_flag[ix][0] = true;
        dCdt_flag[ix][0]     =  true;


        for(size_t jx = 1; jx < CBinSize; jx++){
            for(size_t kx = 0; kx < nSpecies + NSCALAR; kx++) {
                // incorrect, C is based on Tu and Tb, mass fractions
                // may be non-linear in C space
                scalarZC[ix][jx][kx] = 0.0;
            }
        }
    }
}
/**----------------------------------------------------------------------------
 *
 * @brief       Re-initializes the scalarZ[Zbin][Cbin][scalar] with (linear combination) equilibrium values, if a flag bin[][][] is false
 *
 * @author      A.Menon, 13th Jan 2021
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarZC_linear()
{
    auto gas = lem.gas;

    double *Y = lem.tmp_s1;

    for(size_t ix = 0; ix < ZBinSize; ix++){
        double Z  = Zstep * (ix);
        double Tu = Z * lem.fuelStreamT + (1.- Z)*lem.oxStreamT;

        // burnt Ts and Ys
        double *Yb = lem.YEquilZ[ix].origin();
        double *Yu = lem.YUnbZ[ix].origin();
        double Tb = lem.TEquil_z[ix];

        for(size_t jx = 0; jx < CBinSize; jx++){
                double C  = Cstep * jx;
                double T                 = C*Tb + (1.-C)*Tu; // correct
                scalarZC[ix][jx][IT]     = T;
              //scalarZC[ix][jx][ITheta] = 0; // also fine since we're not using it
              //scalarZC[ix][jx][IC]     = 0; // fine, Z alone conditioning


                for(size_t kx = 0; kx < lem.ns; kx++) {
                    scalarZC[ix][jx][IY + kx] = C*Yb[kx] + (1.-C)*Yu[kx];
                }
                // --just to be sure it's a realizable state
                gas->setState_TPY(T, lem.pressure, scalarZC[ix][jx].origin() + IY);
                gas->getMassFractions(Y);
                std::memcpy(scalarZC[ix][jx].origin() + IY, Y, nSpecies*sizeof(double));
                double enth= gas->enthalpy_mass();

                gas->setState_TPY(298.00, 101325, scalarZC[ix][jx].origin() + IY);
                double hForm = gas->enthalpy_mass();

                scalarZC[ix][jx][IRho] = gas->density();
                scalarZC[ix][jx][Ih]   = enth -hForm;
                scalarZC_flag[ix][jx] = true;
        }
    }
}

void LEMLINE::initScalarZC_unburnt()
{
    auto gas = lem.gas;

    double *Y = lem.tmp_s1;

    for(size_t ix = 0; ix < ZBinSize; ix++){
        double Z  = Zstep * (ix);
        double Tu = Z * lem.fuelStreamT + (1.- Z)*lem.oxStreamT;

        // burnt Ts and Ys
        double *Yb = lem.YEquilZ[ix].origin();
        double *Yu = lem.YUnbZ[ix].origin();
        double Tb = lem.TEquil_z[ix];

        for(size_t jx = 0; jx < CBinSize; jx++){

                scalarZC[ix][jx][IT]     = Tu;
                scalarZC[ix][jx][ITheta] = 0; // also fine since we're not using it
                scalarZC[ix][jx][IC]     = 0; // fine, Z alone conditioning

                for(size_t kx = 0; kx < nSpecies; kx++) {
                    scalarZC[ix][jx][IY + kx] = Yu[kx];
                }
                // --just to be sure it's a realizable state
                gas->setState_TPY(Tu, lem.pressure, scalarZC[ix][jx].origin() + IY);
                gas->getMassFractions(Y);
                std::memcpy(scalarZC[ix][jx].origin() + IY, Y, nSpecies*sizeof(double));
                double enth= gas->enthalpy_mass();

                gas->setState_TPY(298.00, 101325, scalarZC[ix][jx].origin() + IY);
                double hForm = gas->enthalpy_mass();

                scalarZC[ix][jx][IRho] = gas->density();
                scalarZC[ix][jx][Ih]   = enth -hForm;
         // }
        }
    }
}

/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarZ[Zbin][Cbin][scalar] with equilibrium
 * values, simple linear in C for T ; read 0D reactor solution for Yi;rest are
 * 0
 *
 * @author      A.Menon, 20th Nov
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarZC_read()
{

    auto gas = lem.gas;
    int RSize = scalarZC[0][0].size();

    //  --- read file, check header

    std::ifstream file;
    file.open("scalarZC.dat" );
    if(!file.is_open()){
        std::cout << "\nERROR!! file scalarZC.dat not found";
        file.close();
        exit(0);
    }

    int fZbins, fCbins,RS;
    file >> fZbins >> fCbins >>RS;
    if( (fZbins != ZBinSize) || (fCbins != CBinSize) || (RS != RSize))
    {
        std::cerr << "\nERROR scalarZC.dat has incorrect sizes:" << fZbins <<' ' <<fCbins <<' ' << RS<<'\n';
        file.close();
        exit(0);
    }
 // else std::cout<< "\nSizes match! "<< fZbins <<' ' <<fCbins <<' ' << RS<<'\n';

    // -- read scalarZC file
    std::vector<double> readData;
    double num;
    while (file >> num) {
        readData.push_back(num);
    }
    file.close();

    // -- read flag file

    std::vector<bool> readDataFlag;
    file.open("scalarZC_flag.dat" );
    if(!file.is_open()){
        std::cout << "\nERROR!! file scalarZC_flag.dat not found";
        file.close();
    }
    else
    while (file >> num) {
        readDataFlag.push_back(num);
    }
    file.close();

    size_t index = 0;
    size_t indexFlag = 0;
    for(int ix = 0;ix < ZBinSize; ix++){
        for(int jx = 0; jx < CBinSize; jx++){
            scalarZC_flag[ix][jx] = readDataFlag[indexFlag];
            dCdt_flag[ix][jx] = readDataFlag[indexFlag];
            indexFlag++;
            for(size_t kx = 0; kx < RSize; kx++){
                scalarZC[ix][jx][kx]= readData[index];
                index++;
            }
        }
    }


}

/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarC[Cbin][scalar]
 *
 * @author      A.Menon 2/1/2024
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarC_read()
{

    auto gas = lem.gas;
    int RSize = scalarC[0].size();

    //  --- read file, check header

    std::ifstream file;
    file.open("scalarC.dat" );
    if(!file.is_open()){
        std::cout << "\nERROR!! file scalarC.dat not found";
        file.close();
        exit(0);
    }

    int fCbins,RS;
    file >> fCbins >>RS;
    if( (fCbins != CBinSize) || (RS != RSize))
    {
        std::cerr << "\nERROR scalarC.dat has incorrect sizes:" << fCbins <<' ' << RS<<'\n';
        file.close();
        exit(0);
    }

    // -- read scalarZC file
    std::vector<double> readData;
    double num;
    while (file >> num) {
        readData.push_back(num);
    }
    file.close();

    // -- read flag file

    size_t index = 0;
    size_t indexFlag = 0;
    for(int jx = 0; jx < CBinSize; jx++){
        indexFlag++;
        for(size_t kx = 0; kx < RSize; kx++){
            scalarC[jx][kx]= readData[index];
            index++;
        }
    }

}

/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarZC[Zbin][Cbin][scalar] with equilibrium
 * values, simple linear in C for T ; inits 0D reactor for each Z_bin and steps
 * to equilibrium
 * @author      A.Menon, 20th Nov 2021
 *
 *---------------------------------------------------------------------------*/
//void LEMLINE::initScalarZC_0D(double flamLimit)
//{
//
//    std::cout << "INSIDE 0D Reactor!! ";
//    std::cout << "flammability limit  Z = "<< flamLimit<<"\n";
//
//    //initScalarZC_linear();
//
//    // -- set scalarZC to unburnt values in case of holes and flammability
//    // -- limit
//
//    for(size_t iz = 0; iz < ZBinSize; iz++){
//        double Z  = iz * Zstep;
//        double Tu = Z*lem.fuelStreamT + (1.- Z)*lem.oxStreamT;
//      //std::vector<double> &Yu = lem.YUnbZ[iz];
//        double *Yu = lem.YUnbZ[iz].origin();
//
//        for(size_t jx = 0; jx < CBinSize; jx++){
//             scalarZC[iz][jx][IT]   = Tu;
//             for(size_t kx = 0; kx < nSpecies; kx++){
//                scalarZC[iz][jx][IY+kx] = Yu[kx];
//             }
//        }
//        scalarZC_flag[iz][0] = true; // all 'found'
//        dCdt_flag[iz][0] = true; // all 'found'
//    }
//
//    auto gas = lem.gas;
//    Cantera::ConstPressureReactor& comb = lem.combustor;
//    comb.setInitialVolume(1.0);
//    Cantera::ReactorNet& rNet =  lem.sim;
//
//    // -- for reaction rates
//    double* Ms   = lem.Ms;
//    double* RR   = lem.tmp_s1;
//
//    gas->setState_TPX(300.0, 101325, "H:1.0");
//    Cantera:: Reservoir igniter;
//    igniter.insert(gas);
//
//    // The igniter will use a Gaussian 'functor' object to specify the
//    // time-dependent igniter mass flow rate.
//    double A = 0.01;
//    double FWHM = 0.2;
//    double t0 = 0.05;
//    Cantera::Gaussian igniter_mdot(A, t0, FWHM);
//
//    Cantera::MassFlowController m3;
//    m3.install(igniter, comb);
//    m3.setFunction(&igniter_mdot);
//
//
//
//    std::vector<double> Ynow(nSpecies,-1.0);
//
//    flamLimit = flamLimit > ZMax ? ZMax: flamLimit;
//
//    size_t ZBinLIM = flamLimit/Zstep;
//
//   // run zeroD reactor for all ZBins except extremes
//    for(size_t iz = 1; iz < ZBinLIM; iz++){
//
//        try{
//
//        double Z  = iz * Zstep;
//        double Tu = Z*lem.fuelStreamT + (1.- Z)*lem.oxStreamT;
//
//     // std::vector<double> Yu = lem.YUnbZ[iz];
//        double *Yu = lem.YUnbZ[iz].origin();
//        double Tb = lem.TEquil_z[iz]; // final state
//
//
//
//        gas->setState_TPY(Tu,lem.pressure,Yu);
//        double rho = gas->density();
//
//        // -- set up 0D reactor; (1) ySol
//        ySol[0] = rho;
//        ySol[1] = gas->enthalpy_mass()*rho;
//        gas->getMassFractions(ySol.data()+2);
//
//        // -- set up 0D reactor; (2)
//        comb.updateState(ySol.data());
//        rNet.setInitialTime(0.0);
//
//        // -- conditioning
//        std::vector<int>    hit(CBinSize, 0);
//        std::vector<std::vector<double>>    tempLoc(CBinSize);
//        for (size_t kx = 0; kx < CBinSize;kx++){
//            tempLoc[kx].resize(lem.ns + 2, 0.0);
//        }
//
//        double YBF = 0.0;
//        double YUF = 0.0;
//        for (size_t kx = 0; kx < lem.C_select.size();kx++){
//            YBF += lem.YEquilZ[iz][lem.C_select[kx]];
//            YUF +=   lem.YUnbZ[iz][lem.C_select[kx]];
//        }
//
//        double time = 0.0;
//        double Tnow = comb.temperature();
//        double cEnd = 1.0 - 0.5/CBinSize;
//
//        while(rNet.time() < 300){
//          //t += 1e-5;
//          //rNet.advance(t);
//          //rNet.step();
//            gas->getMassFractions(Ynow.data());
//            lem.kin->getNetProductionRates(RR); // units Kmol/m^3/s
//            Tnow = comb.temperature();
//            std::vector<double> CC1= getProgVble_nonpremixSp(Z,Ynow.data(), Tnow);
//
//            double Clocal = CC1[0];
//            double dCdt = CC1[1];
//            // -- find suitable bin, add
//            int posCy = std::round(Clocal*(CBinSize-1));
//
//            tempLoc[posCy][0] += Tnow;
//            tempLoc[posCy][1] += dCdt;
//
//            for(size_t kx = 0; kx < lem.ns; kx++){
//                tempLoc[posCy][kx + 2] += Ynow[kx];
//            }
//            hit[posCy] ++;
//            rNet.step();
//        }
//
//        // -- exit 0D reactor average over Chits
//
//        int nope = 0;
//        for(size_t jx = 0; jx < CBinSize; jx++){
//            if(hit[jx] > 0){
//                scalarZC[iz][jx][IT]     = tempLoc[jx][0]/hit[jx];
//                scalarZC[iz][jx][ITheta] = tempLoc[jx][1]/hit[jx];
//                scalarZC[iz][jx][IC2]    = tempLoc[jx][1]/hit[jx];
//                for(size_t kx = 0; kx < nSpecies; kx++){
//                    scalarZC[iz][jx][kx+IY]  = tempLoc[jx][kx+2]/hit[jx];
//                }
//                scalarZC_flag[iz][jx]  = true;
//                dCdt_flag[iz][jx]  = true;
//            }
//            else {nope++;
//            scalarZC_flag[iz][jx]  = false;
//            dCdt_flag[iz][jx]  = false;
//            }
//        }
//       std::cout << " Zbin | Zmean " << iz<<" " << Z<<", unfilled Cbins= "<< nope <<" ";
//       std::cout << " max T  = " << Tnow << "\n" ;
//
//        }catch (Cantera::CanteraError& err) {
//        // handle exceptions thrown by Cantera
//            std::cout << err.what() << std::endl;
//            std::cout << "Failed,  trying next Z bin \n";
//        }
//    }
//}

/**----------------------------------------------------------------------------
 * @brief       reads CFD data, cell by cell, fills the scalarZC[][][] sum only
 *
 * @author      A.Menon, 27th Sep 2023
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarZC_fromCFDSum(double Z, double C, double T, double *Y, double dcDt)
{
 // int    posZx, posCy;

 // auto gas = lem.gas;


 // posZx = std::floor(Z/Zstep);
 // posCy = std::floor(C*(CBinSize-1));

 // lem.tempZC[posZx][posCy][IT    ] += T;
 // tempZC[posZx][posCy][IRho  ] += (*it).rho;
 // lem.tempZC[posZx][posCy][ITheta] += dcDt;

 // for(size_t kx = 0; kx < nSpecies; kx++){
 //     lem.tempZC[posZx][posCy][kx + IY] += Y[kx];
 // }

 ///lem.hitCountZC[posZx][posCy] += 1.0;
}
/**----------------------------------------------------------------------------
 * @brief       reads CFD data, cell by cell, fills the scalarZC[][][], divides
 * by hits
 *
 * @author      A.Menon, 27th Sep 2023
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::initScalarZC_fromCFDNorm(){

    // average over number of hits
  //for(size_t ix = 0; ix < ZBinSize; ix++){
  //    for(size_t jx = 0; jx < CBinSize; jx++){
  //        if(lem.hitCountZC[ix][jx] > 0.0){

  //            scalarZC_flag[ix][jx] = true; // update the global flag matrix
  //            dCdt_flag[ix][jx] = true; // update the global flag matrix

  //            // blending + persistence
  //            double Y;
  //            // temp
  //            Y = lem.tempZC[ix][jx][IT]/lem.hitCountZC[ix][jx];
  //            scalarZC[ix][jx][IT] = Y;

  //            // dcdt
  //            Y = lem.tempZC[ix][jx][ITheta]/lem.hitCountZC[ix][jx];
  //            scalarZC[ix][jx][ITheta]=Y;

  //            // species
  //            for(size_t kx = 0; kx < nSpecies; kx++){
  //                Y = lem.tempZC[ix][jx][kx+IY]/lem.hitCountZC[ix][jx];
  //                scalarZC[ix][jx][kx+IY] = Y;
  //            }
  //        }
  //    }
  //}
}

/**----------------------------------------------------------------------------
 * @brief       Initializes the scalarZC[Zbin][Cbin][scalar] with equilibrium
 * values, simple linear in C for T ; inits 0D reactor for each Z_bin and steps
 * to equilibrium
 * @author      A.Menon, 20th Nov 2021
 *
 *---------------------------------------------------------------------------*/
//void LEMLINE::initScalarC_0D()
//{
//
//    initScalarC();
//
//    auto gas = lem.gas;
//    Cantera::ConstPressureReactor& comb = lem.combustor;
//    comb.setInitialVolume(1.0);
//    Cantera::ReactorNet& rNet =  lem.sim;
//
//    // -- conditioning
//    double hit[CBinSize]{};
//
//    // 0-T, 1-dCdT, 2... species
//    double temp[CBinSize][nSpecies+2]{};
//
//
//    // the ingitor for complex mechanisms
//    int indH = gas->speciesIndex("H");
//    Cantera:: Reservoir igniter;
//
//    // The igniter will use a Gaussian 'functor' object to specify the
//    // time-dependent igniter mass flow rate.
//    double A = 0.01;
//    double FWHM = 0.02;
//    double t0 = 0.01;
//    Cantera::Gaussian igniter_mdot(A, t0, FWHM);
//    Cantera::MassFlowController m3;
//
//    if(indH != -1){
//        gas->setState_TPX(300.0, 101325, "H:1.0");
//        igniter.insert(gas);
//        std::cout << "0D Reactor: H radical present, using ignitor..\n";
//        m3.install(igniter, comb);
//        m3.setFunction(&igniter_mdot);
//    }
//
//    // -- for reaction rates
//    double* Ms   = lem.Ms;
//    double* RR   = lem.tmp_s1;
//
//    gas->setState_TPY(lem.fuelStreamT,lem.pressure,Yu_p.data());
//
//    double t = 0;
//    double Tnow = comb.temperature();
//    double rho = gas->density();
//
//    // -- set up 0D reactor; (1) ySol
//    ySol[0] = rho;
//    ySol[1] = gas->enthalpy_mass()*rho;
//    gas->getMassFractions(ySol.data()+2);
//
//    // -- set up 0D reactor; (2)
//    comb.updateState(ySol.data());
//    rNet.setInitialTime(0);
//
//    double Yu = 0.0;
//    double Yb = 0.0;
//    std::vector<double> Ynow(lem.ns,0.0);
//    for (size_t kx = 0; kx < lem.C_select.size();kx++){
//        Yb += Yb_p[lem.C_select[kx]];
//        Yu += Yu_p[lem.C_select[kx]];
//    }
//
//
//    double C= 0.0;
//  //while(C< 0.999){
//    while(Tnow < Tb_p-1e-7){
//  //while(rNet.time() < 30000){
//
//        gas->getMassFractions(Ynow.data());
//        lem.kin->getNetProductionRates(RR); // units Kmol/m^3/s
//        Tnow = comb.temperature();
//
//        double rho = gas->density();
//        double dCdt = 0.0;
//
//        C = 0.0;
//        // -- find progress variable
//        for (size_t kx = 0; kx < lem.C_select.size();kx++){
//            size_t oo = lem.C_select[kx];
//            C    += Ynow[oo];
//            dCdt += RR[oo]*Ms[oo]/gas->density()/(Yb - Yu);
//        }
//    //  std::cout << "0d advanced to " << t << " Y = "<< C <<" Yb =  "<< Yb <<" Yu = " << Yu;
//        C = (C-Yu)/(Yb-Yu);
//   //   std::cout << " C = " << C << "\n";
//
//
//        C = C > 1.0? 1.0:C;
//        C = C < 0.0? 0.0:C;
//
//        // -- find suitable bin, add
//        size_t posCy = std::round(C*(CBinSize-1));
//        temp[posCy][0] += Tnow;
//        temp[posCy][1] += dCdt;
//
//        for(size_t kx = 0; kx < nSpecies; kx++){
//            temp[posCy][kx + 2] += Ynow[kx];
//        }
//        scalarC_flag[posCy] = 1.0;
//        hit[posCy] += 1.0;
//        rNet.step();
//      //t += 0.000005;
//      //rNet.advance(t);
//
//    }
//
//    // -- exit 0D reactor average over Chits
//    int nope = 0;
//    for(size_t jx = 0; jx < CBinSize; jx++){
//        if(hit[jx] > 0){
//            scalarC[jx][IT ] = temp[jx][0]/hit[jx];
//            scalarC[jx][IC ] = temp[jx][1]/hit[jx];
//            scalarC[jx][IC2] = temp[jx][1]/hit[jx]; // stores the initial value of dCdT
//            for(size_t kx = 0; kx < nSpecies; kx++){
//                scalarC[jx][kx+IY]  = temp[jx][kx+2]/hit[jx];
//            }
//        }
//        else nope++;
//    }
//   std::cout << "unfilled c bins= "<< nope <<"\n";
//}


/**---------------------------------------------------------------------------
 * @brief      Conditions all scalars on Z alone scalarZ[Zbin][scalar];
 * conventional averaging, persistence and blending
 *
 * @author      A.Menon, 18th Nov
 *
 * UPDATE - Now being used for integration skips (when no combustion occurs)
 *--------------------------------------------------------------------------*/

void LEMLINE::calcScalarZ()
{
    int pos;
    double hitCount[ZBinSize]{};
    double tempZ[ZBinSize][nSpecies+NSCALAR]{};
    auto gas= lem.gas;

    std::list<CellData>::iterator it;

    // -- initialize temp to Zero
    for(size_t ix = 0; ix < ZBinSize; ix++){

        for(size_t kx = 0; kx < nSpecies + NSCALAR;kx++){
            tempZ[ix][kx] = 0.0;
        }
    }

    // collect Ys for corresponding bin
    for (it = cells.begin(); it != cells.end(); ++it){
        pos = std::round((*it).Z/Zstep);
        hitCount[pos]++;
        tempZ[pos][IT]     += (*it).T;
      //lem.tempZ[pos][IRho]   += (*it).rho;
      //lem.tempZ[pos][ITheta] += (*it).dedt;
      //lem.tempZ[pos][IC]     += (*it).C;

        for(size_t jx = 0; jx < nSpecies; jx++){
            tempZ[pos][IY + jx] += (*it).Y[jx];
        }

        // calc sensible enthalpy for bin
        gas->setState_TPY(298, 101325, (*it).Y);
        double hForm = gas->enthalpy_mass(); //mean formation
        tempZ[pos][Ih ]    += (*it).E - hForm;
    }

    // persistence + blending update for YcondZ[][]
    for(size_t ix = 0; ix < ZBinSize; ix++){
        if(hitCount[ix] > 0){

            scalarZ[ix][IT    ] = bFac*(tempZ[ix][IT    ]/hitCount[ix])+(1.0-bFac)*
            scalarZ[ix][IT    ];

            scalarZ[ix][IRho  ] = bFac*(tempZ[ix][IRho  ]/hitCount[ix])+(1.0-bFac)*
            scalarZ[ix][IRho  ];

            scalarZ[ix][Ih ] = bFac*(tempZ[ix][Ih ]/hitCount[ix])+(1.0-bFac)*
            scalarZ[ix][Ih ];

          //scalarZ[ix][ITheta] = bFac*(lem.tempZ[ix][ITheta]/hitCount[ix])+(1.0-bFac)*
          //scalarZ[ix][ITheta];

            // -- use condProgVble() for now
          //scalarZ[ix][IC    ] = bFac*(lem.tempZ[ix][IC    ]/hitCount[ix])+(1.0-bFac)*
          //scalarZ[ix][IC    ];

            for(size_t jx = 0; jx < nSpecies; jx++)
                scalarZ[ix][IY + jx] = bFac*(tempZ[ix][IY + jx] / hitCount[ix])+(1.0-bFac)*
                scalarZ[ix][IY + jx];
        }
    }

}

/**---------------------------------------------------------------------------
 * @brief      Conditions all scalars on C alone scalarC[CBin][scalar];
 * conventional averaging, persistence and blending
 *
 * @author      A.Menon, 26th Nov 2021
 *
 *--------------------------------------------------------------------------*/
void LEMLINE::calcScalarC()
{
    int pos;
    double hitCount[CBinSize]{};
    auto gas= lem.gas;

    std::list<CellData>::iterator it;

    // -- initialize temp to Zero
    for(size_t ix = 0; ix < CBinSize; ix++)
        for(size_t jx = 0; jx < lem.ns + NSCALAR;jx++)
            tempC[ix][jx] = 0.0;

    // collect Ys for corresponding bin
    for (it = cells.begin(); it != cells.end(); ++it){
        pos = std::round((*it).C/Cstep);
        pos = std::min(pos,CBinSize-1);
        pos = std::max(pos,0);
        hitCount[pos]++;
        tempC[pos][IT]     += (*it).T;
        tempC[pos][IRho]   += (*it).rho;
        tempC[pos][ITheta] += (*it).dedt;
        tempC[pos][IC]     += (*it).dCdt;

        for(size_t jx = 0; jx < nSpecies; jx++){
            (*it).Y[jx] = (*it).Y[jx] < 0? 0.0: (*it).Y[jx];
            tempC[pos][IY+jx] += (*it).Y[jx];
        }

        // calc sensible enthalpy for bin
        gas->setState_TPY(298, 101325, (*it).Y);
        double hForm = gas->enthalpy_mass(); //mean formation
        tempC[pos][Ih ]    += (*it).E - hForm;
    }

    // persistence + blending update for YcondZ[][]
    for(size_t ix = 0; ix < CBinSize; ix++){
        if(hitCount[ix] > 0){

            scalarC_flag[ix] = 1.0; // update the global flag matrix
            scalarC[ix][IT]     = bFac*(tempC[ix][IT]/hitCount[ix])+(1.0-bFac)*
            scalarC[ix][IT];

            scalarC[ix][IRho]   = bFac*(tempC[ix][IRho]/hitCount[ix])+(1.0-bFac)*
            scalarC[ix][IRho];

            scalarC[ix][Ih]     = bFac*(tempC[ix][Ih]/hitCount[ix])+(1.0-bFac)*
            scalarC[ix][Ih];

            scalarC[ix][ITheta] = bFac*(tempC[ix][ITheta]/hitCount[ix])+(1.0-bFac)*
            scalarC[ix][ITheta];

            scalarC[ix][IC]     = bFac*(tempC[ix][IC]/hitCount[ix])+(1.0-bFac)*
            scalarC[ix][IC];

            for(size_t jx = 0; jx < lem.ns; jx++){
                scalarC[ix][IY+jx] = bFac*(tempC[ix][IY+jx]/hitCount[ix])+(1.0-bFac)*
                scalarC[ix][IY+jx];

            }

        }

    }


}


/**---------------------------------------------------------------------------
 * @brief      Conditions specied production rates on C ; conventional
 * averaging, persistence and blending
 *
 * @author      A.Menon, 1st Dec 2021
 *
 *--------------------------------------------------------------------------*/
void LEMLINE::calcSpeciesC()
{
    int pos;
    int hitCount[CBinSize]{};
    double   RR[nSpecies];
    double*  Ms=lem.Ms;

    auto gas= lem.gas;

    std::list<CellData>::iterator it;

    // -- initialize temp to Zero
    for(size_t ix = 0; ix < CBinSize; ix++)
        for(size_t jx = 0; jx < lem.ns ;jx++)
            tempC[ix][jx] = 0;

    // collect Ys for corresponding bin
    for (it = cells.begin(); it != cells.end(); ++it){
        pos = std::round((*it).C*(CBinSize-1));

        // calc sensible enthalpy for bin
        gas->setState_TPY((*it).T, (*it).p, (*it).Y);
        lem.kin->getNetProductionRates(RR);
        // just a test

        hitCount[pos]++;
        for(size_t jx = 0; jx < nSpecies; jx++){
            tempC[pos][jx] += RR[jx]*Ms[jx];
        }

    }

    // -- presist + blend
    for(size_t ix = 0; ix < CBinSize; ix++){
        if(hitCount[ix] > 0){
            for(size_t jx = 0; jx < nSpecies; jx++)
                speciesProdC[ix][jx] = bFac*(tempC[ix][jx]/hitCount[ix])+
                                      (1.0-bFac)* speciesProdC[ix][jx];
        }
    }

}




/**---------------------------------------------------------------------------
 * @brief      Conditions progress variable on Z alone scalarZ[Zbin][IC];
 * conventional averaging, persistence and blending
 *
 * @author      A.Menon, 19th Nov
 *--------------------------------------------------------------------------*/

int LEMLINE::condProgVble()
{
    int pos;
    double temp[ZBinSize]{};
    double hitCount[ZBinSize]{};

    std::list<CellData>::iterator it;

    // collect Ys for corresponding bin
    for (it = cells.begin(); it != cells.end(); ++it){
        pos = (*it).Z/Zstep;
        hitCount[pos]++;
     // temp[pos]   += (*it).dCdt;
        temp[pos]   += (*it).C;

    }

    int hits = 0;

  //double Cdef = getMeanS()[IC];
    // persistence + blending update for scalarZ[Zbin][IC]
    for(size_t ix = 0; ix < ZBinSize; ix++){
        scalarZ[ix][IC] = 0;
        if(hitCount[ix] >0)
            scalarZ[ix][IC] = bFac*(temp[ix]/hitCount[ix])+(1.0-bFac)*
            scalarZ[ix][IC];
        else hits++;
    }

    return hits;
}

/**---------------------------------------------------------------------------
 * @brief       Calculates the matrix scalarZC[ZBin][CBin][scalar] using hybrid
 * persistence-temporal blending
 *
 * @author      A.Menon 19th Nov
 *--------------------------------------------------------------------------*/
void LEMLINE::calcScalarZC()
{

    // -- condition progress variable only on Z
//  condProgVble();

    int    posZx, posCy;

    double hitCountZC[ZBinSize][CBinSize]{};
    double tempZC[ZBinSize][CBinSize][nSpecies+NSCALAR]{};

    auto gas = lem.gas;

    std::list<CellData>::iterator it;

    // increment the "found" buckets
    for (it = cells.begin(); it != cells.end(); ++it){

        posZx = std::round((*it).Z/Zstep);
        posZx = std::min(posZx, ZBinSize-1);

        posCy = std::round((*it).C*(CBinSize-1));
        posCy = std::min(posCy, CBinSize-1);

        tempZC[posZx][posCy][IT    ] += (*it).T;
        tempZC[posZx][posCy][IRho  ] += (*it).rho;

        // -- comment out for Z alone conditioning
        tempZC[posZx][posCy][IC] += (*it).dCdt;
        tempZC[posZx][posCy][ITheta] += (*it).dedt;

        for(size_t kx = 0; kx < nSpecies; kx++){
            tempZC[posZx][posCy][kx + IY] += (*it).Y[kx];
        }

        // calc sensible enthalpy for bin
      //gas->setState_TPY(298, 101325, (*it).Y);
      //double hForm = gas->enthalpy_mass(); //mean formation
      //lem.tempZC[posZx][posCy][Ih ]  += (*it).E - hForm;
        hitCountZC[posZx][posCy] += 1.0;

    }

    // average over number of hits
    for(size_t ix = 0; ix < ZBinSize; ix++){
        for(size_t jx = 0; jx < CBinSize; jx++){
            if(hitCountZC[ix][jx] > 0.0){

                scalarZC_flag[ix][jx] = true; // update the global flag matrix
                dCdt_flag[ix][jx] = true; // update the global flag matrix

                // blending + persistence
                double Y;
                Y = tempZC[ix][jx][IT    ]/hitCountZC[ix][jx];

                scalarZC[ix][jx][IT    ] = bFac*Y + (1.-bFac)*
                scalarZC[ix][jx][IT    ];

                Y = tempZC[ix][jx][IRho  ]/hitCountZC[ix][jx];
                scalarZC[ix][jx][IRho  ] = bFac*Y + (1.-bFac)*
                scalarZC[ix][jx][IRho  ];

              //Y = lem.tempZC[ix][jx][Ih ]/lem.hitCountZC[ix][jx];
              //scalarZC[ix][jx][Ih ] = bFac*Y + (1.-bFac)*
              //scalarZC[ix][jx][Ih ];


                Y = tempZC[ix][jx][IC    ]/hitCountZC[ix][jx];
                scalarZC[ix][jx][IC    ] = bFac*Y + (1.-bFac)*
                scalarZC[ix][jx][IC];

                Y = tempZC[ix][jx][ITheta]/hitCountZC[ix][jx];
                scalarZC[ix][jx][ITheta] = bFac*Y + (1.-bFac)*
                scalarZC[ix][jx][ITheta];

                for(size_t kx = 0; kx < nSpecies; kx++){
                    Y = tempZC[ix][jx][kx+IY]/hitCountZC[ix][jx];
                    scalarZC[ix][jx][kx+IY] = bFac*Y + (1.-bFac)*
                    scalarZC[ix][jx][kx+IY];
                }
            }
        }
    }
}

/**---------------------------------------------------------------------------
 * Place scalarC into appropriate Re_t bin. Re-using scalarZC for this.
 *
 *--------------------------------------------------------------------------*/
void LEMLINE::calcScalarReC(double Re_t)
{

    int posRex;
    int RSize = scalarC[0].size();
    posRex = std::round(Re_t/Zstep);
    for(size_t jx = 0; jx < CBinSize; jx++){
        std::memcpy(scalarZC[posRex][jx].origin(), scalarC[jx].origin(), RSize*sizeof(double));
        scalarZC_flag[posRex][jx] = true;
    }

}

/**---------------------------------------------------------------------------
 * Interpolate scalarZC (Z is Re_) to compute scalarC that can be used for
 * integration
 *
 *--------------------------------------------------------------------------*/
void LEMLINE::interpCRe(double Re_t)
{
    int posRex;
    int RSize = scalarC[0].size();
    posRex = std::floor(Re_t/Zstep);
    posRex = std::max(0,posRex);
    posRex = std::min(ZBinSize-1,posRex);
    int posRexIncr = std::min(ZBinSize-1, posRex+1);

    double diffRe = Re_t - posRex*Zstep;

    for(size_t ix = 0; ix < CBinSize; ix++){
        for(size_t jx = 0; jx < RSize; jx++){
            double m = (scalarZC[posRexIncr][ix][jx] - scalarZC[posRex][ix][jx])/Zstep;
            double sum = scalarZC[posRex][ix][jx] + m*diffRe;
            scalarC[ix][jx] = sum;

        }
    }

}


/**----------------------------------------------------------------------------
 * @brief       Sets the progress variable for each cell using
 *              adiabatic/equilibrium flame temperature
 *
 *              OBS linear interpolation flag added
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/
void LEMLINE::setProgressVariable_T(bool flagLinear, double dt)
{
    std::list<CellData>::iterator it;
    auto gas = lem.gas;
    for(it = cells.begin(); it != cells.end(); it++)
    {
        double T = (*it).T;
        double Z = (*it).Z;

        // This seems to prevent the false C values affected by small diffusive
        // heat fluxes in singularity regions
        if(Z < 1E-9 || Z > 0.9999999999) { (*it).C = 0.0; continue;}

        int     pos  = Z/Zstep;
        double  TUnb = Z*lem.fuelStreamT + (1-Z)*lem.oxStreamT;
        double  TB;

        if(flagLinear){
            TB = lem.TEquil_z[pos] + std::fmod(Z, Zstep)*(lem.TEquil_z[pos +1] -
                 lem.TEquil_z[pos]) / Zstep;
        }

        else{
            TB = lem.TEquil_z[pos];
        }

        double C = (T - TUnb) / (TB - TUnb);

        // Limits
        if (C > 1.0 ) C = 1.0;
        if (C < 0.0 ) C = 0.0;

        double crate =  (C - (*it).C)/dt;
        (*it).C = C;
        (*it).dCdt = crate;

        // some trial.. possible instability
     // gas->setState_TPY((*it).T, (*it).p, (*it).Y);
     // double cv = gas->cv_mass();
     // (*it).dCdt = -1 * (*it).dedt/(TB-TUnb)/cv;
    }
}

/**----------------------------------------------------------------------------
 * @brief       Sets the dCdt for every cell, this is only valid if the
 *              progress variable is based on temperature,
 *              dCdt = dedt/cp/(Tb- Tu)
 *
 *              needs dedt to be calculated first
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/

void LEMLINE::setdCdt(){

    std::list<CellData>::iterator it;
    for(it = cells.begin(); it != cells.end(); it++){

    }
}


/**----------------------------------------------------------------------------
 * @brief       Sets the progress variable for each cell using  mass fraction
 *              of CO2 and H2O (Pierce and Moin)
 *
 * @author      A.Menon 22nd Dec
 *---------------------------------------------------------------------------*/
void LEMLINE::setProgressVariable_simple(double dt)
{

    std::list<CellData>::iterator it;
    auto gas = lem.gas;
    int indCO2 = gas->speciesIndex("CO2");
    int indH2O = gas->speciesIndex("H2O");

    if(indCO2 != -1 && indH2O != -1)
    for(it = cells.begin(); it != cells.end(); it++)
    {
        double Z = (*it).Z;

        // This seems to prevent the false C values affected by small diffusive
        // heat fluxes in singularity regions
        if(Z < 1E-9 || Z > 0.9999999999) { (*it).C = 0.0; continue;}

        int     pos  = Z/Zstep;

        double Y = (*it).Y[indCO2] + (*it).Y[indH2O];
        double YBCO2, YBH2O;

        if(pos < ZBinSize -1){
            YBCO2 = lem.YEquilZ[pos][indCO2] + std::fmod(Z, Zstep)*
                (lem.YEquilZ[pos +1][indCO2] - lem.YEquilZ[pos][indCO2]) / Zstep;
            YBH2O = lem.YEquilZ[pos][indH2O] + std::fmod(Z, Zstep)*
                (lem.YEquilZ[pos +1][indH2O] - lem.YEquilZ[pos][indH2O]) / Zstep;
        }
        else{
            YBCO2 = lem.YEquilZ[ZBinSize -1][indCO2];
            YBH2O = lem.YEquilZ[ZBinSize -1][indH2O];
        }

        double YB = YBCO2 + YBH2O;
        double C = YB>0? Y/YB:0.0;

        // Limits
     // if (C > 1.0 ) C = 1.0;
     // if (C < 0.0 ) C = 0.0;

        (*it).C = C;

        // since this is called after advanceSourceTerm
        double rate = (C - (*it).C)/(dt);
   //   double rate = (C)/dt;
        (*it).dCdt = rate;

    }
    else if (indCO2 == -1){
    for(it = cells.begin(); it != cells.end(); it++)
    {
        double Z = (*it).Z;

        // This seems to prevent the false C values affected by small diffusive
        // heat fluxes in singularity regions
        if(Z < 1E-9 || Z > 0.9999999999) { (*it).C = 0.0; continue;}

        int     pos  = Z*(ZBinSize-1);

        double Y = (*it).Y[indH2O];
        double YBH2O;

        if(pos < ZBinSize -1){
            YBH2O = lem.YEquilZ[pos][indH2O] + std::fmod(Z, Zstep)*
                (lem.YEquilZ[pos +1][indH2O] - lem.YEquilZ[pos][indH2O]) / Zstep;
        }
        else{
            YBH2O = lem.YEquilZ[ZBinSize -1][indH2O];
        }

        double YB = YBH2O;
        double C  = YB>0? Y/YB:0.0;

        // Limits
        if (C > 1.0 ) C = 1.0;
        if (C < 0.0 ) C = 0.0;

        (*it).C = C;

        // since this is called after advanceSourceTerm
        double rate = (C - (*it).C)/(dt);
    /// double rate = (C)/dt;
        (*it).dCdt = rate;

    }

    }
}


void LEMLINE::setProgVble_premix(double dt)
{

    std::list<CellData>::iterator it;
    auto gas = lem.gas;
    int indCO2 = gas->speciesIndex("CO2");
    int indH2O = gas->speciesIndex("H2O");
    int indO2  = gas->speciesIndex("O2");

    if(indCO2 != -1 && indH2O != -1){

      //double Yu = Yu_p[indCO2] + Yu_p[indH2O];
      //double Yb = Yb_p[indCO2] + Yb_p[indH2O];
      //double Yu = Yu_p[indCO2] ;
      //double Yb = Yb_p[indCO2] ;
        double Yu = Yu_p[indO2] ;
        double Yb = Yb_p[indO2] ;

        for(it = cells.begin(); it != cells.end(); it++)
        {
         // double    Y = (*it).Y[indCO2] + (*it).Y[indH2O];
         // double    Y = (*it).Y[indCO2] ;
            double    Y = (*it).Y[indO2] ;
            double    C = (Y-Yu)/(Yb - Yu);

            double rate = (C - (*it).C)/(dt);

            (*it).C    = C;
            (*it).dCdt = rate;

        }
    }
    else if (indH2O != -1){

        double Yu = Yu_p[indH2O];
        double Yb = Yb_p[indH2O];

        for(it = cells.begin(); it != cells.end(); it++)
        {

            double    Y = (*it).Y[indH2O];
            double    C = (Y-Yu)/(Yb - Yu);
            double rate = (C - (*it).C)/(dt);

            (*it).C    = C;
            (*it).dCdt = rate;
        }
    }

  //else{

  //    std::cout << "ERROR, setProgVble_premix, exiting";
  //    exit(0);
  //}
}



/**----------------------------------------------------------------------------
 * @brief       Sets the progress variable for each cell using  mass fraction
 * of species O2
 *
 * @author      A.Menon 27 Jan 2022
 *---------------------------------------------------------------------------*/
void LEMLINE::setProgVble_premixO2()
{

    std::list<CellData>::iterator it;
    auto gas = lem.gas;

    int indO2 = gas->speciesIndex("O2");
    double* Ms   = lem.Ms;
    double* RR   = lem.tmp_s1;

    double Yu = 0, Yb = 0;
    Yu += Yu_p[indO2];
    Yb += Yb_p[indO2];


    for(it = cells.begin(); it != cells.end(); it++)
    {
        double Y = (*it).Y[indO2];
        double C = (Y-Yu)/(Yb - Yu);

        C = C <0.0? 0.0:C;
        C = C >1.0? 1.0:C;
        (*it).C    = C;


        gas->setState_TPY((*it).T, (*it).p, (*it).Y);
        lem.kin->getNetProductionRates(RR);// units Kmol/m^3/s
        double  rho = gas->density();
        double dYdt = RR[indO2]*Ms[indO2]/rho/(Yb - Yu);

        (*it).dCdt  = abs(dYdt);

    }
}

/**----------------------------------------------------------------------------
 * @brief       Sets the progress variable for each cell using  mass fraction
 * of species in C_select
 *
 * @author      A.Menon 4 Jan 2023
 *---------------------------------------------------------------------------*/
void LEMLINE::setProgVble_premixSp()
{

    std::list<CellData>::iterator it;
    auto gas = lem.gas;

    double* Ms   = lem.Ms;
    double* RR   = lem.tmp_s1;

    double Yu = 0.0, Yb = 0.0;
    for(size_t ix = 0; ix < lem.C_select.size(); ix++){
        Yu += Yu_p[lem.C_select[ix]];
        Yb += Yb_p[lem.C_select[ix]];
    }


    for(it = cells.begin(); it != cells.end(); it++)
    {
        double Y = 0.0;

        for(size_t ix = 0; ix < lem.C_select.size(); ix++){
            Y += (*it).Y[lem.C_select[ix]];
        }

        double C = (Y-Yu)/(Yb - Yu);

        C = C <0.0? 0.0:C;
        C = C >1.0? 1.0:C;
        (*it).C    = C;

        gas->setState_TPY((*it).T, (*it).p, (*it).Y);
        lem.kin->getNetProductionRates(RR);// units Kmol/m^3/s
        double  rho = gas->density();
        double dYdt = 0.0;
        for(size_t ix = 0; ix < lem.C_select.size(); ix++){
            dYdt += RR[lem.C_select[ix]]*Ms[lem.C_select[ix]];
        }

        dYdt = dYdt/rho/(Yb - Yu);
        (*it).dCdt  = dYdt;

    }
}

/**----------------------------------------------------------------------------
 * @brief      returns the progress variable for each cell using  mass fractions and Z
 *
 * @author      A.Menon 4 Jan 2023
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getProgVble_nonpremixSp(double Z, double* Y, double T)
{
    std::vector<double> C_and_Source (2,0.0);
    if(Z <= Zstep/2 || Z >= ZMax -Zstep/2) {
        C_and_Source[0] = 0.0;
        C_and_Source[1] = 0.0;
        return C_and_Source;
    }


    auto gas = lem.gas;
    double Yu;
    double Yb;
    double* Ms   = lem.Ms;
    double* RR   = lem.tmp_s2;


    double Yu1 = 0, Yb1 = 0;
    double Y1 = 0.0;

    // --- unburnt
    int Zpos     = std::floor(Z/Zstep);
    Zpos = Zpos < 0 ? 0: Zpos;
    Zpos = std::min(ZBinSize-1, Zpos);

    int ZposIncr = std::min(ZBinSize-1, Zpos+1);
    double diffZ = Z - Zpos*Zstep;

    for(int cs:lem.C_select){
        Yu= Z*lem.fuelY[cs] + (1-Z)*lem.oxY[cs];

        // interpolate  to find yb, 1st order accuracy
        double m = (lem.YEquilZ[ZposIncr][cs] - lem.YEquilZ[Zpos][cs])/Zstep;
        Yb = lem.YEquilZ[Zpos][cs] + m*diffZ;

        Yu1 += Yu;
        Yb1 += Yb;
        Y1  += Y[cs];
    }

    double C = (Y1-Yu1)/(Yb1 - Yu1);


    C = C <0.0? 0.0:C;
    C = C >1.0? 1.0:C;

    gas->setState_TPY(T, lem.pressure,Y);
    lem.kin->getNetProductionRates(RR);// units Kmol/m^3/s

    double dYdt = 0.0;
    double rho  = gas->density();

    for(size_t ix = 0; ix < lem.C_select.size(); ix++){

        dYdt += RR[lem.C_select[ix]]*Ms[lem.C_select[ix]];
    }

    dYdt = dYdt/(Yb1 - Yu1)/rho;
    dYdt = dYdt < 0. ? 0.: dYdt;

    C_and_Source[0] = C;
    C_and_Source[1] = dYdt;

    return C_and_Source;


}

/**----------------------------------------------------------------------------
 * @brief       Sets the progress variable for each cell using  mass fraction
 * of species in C_select, for non premixed case, dCdt set based on dYdt  of
 * C_select
 *
 * @author      A.Menon 4 Jan 2023
 *---------------------------------------------------------------------------*/
void LEMLINE::setProgVble_nonpremixSp()
{

    std::list<CellData>::iterator it;
    auto gas = lem.gas;

    double* RR   = lem.tmp_s1;

    for(it = cells.begin(); it != cells.end(); it++)
    {
        double  Z = (*it).Z;
        double *Y = (*it).Y;
        double T = (*it).T;
        std::vector<double> C_DcDt = getProgVble_nonpremixSp(Z,Y,T);
        (*it).C    = C_DcDt[0];
        (*it).dCdt = C_DcDt[1];
    }
}

/**----------------------------------------------------------------------------
 * @brief       Sets indices of species in C_select
 * \input       vector of species names
 *
 * @author      A.Menon 4 Jan 2023
 *---------------------------------------------------------------------------*/
void LEMLINE::setC_select(std::vector<std::string> C_list){

    auto gas = lem.gas;
    lem.C_select.resize(C_list.size());
    for(size_t ix=0; ix < C_list.size(); ix++){
        int index = gas->speciesIndex(C_list[ix]);
        if (index == -1){
            std::cout << "ERROR, species" << C_list[ix] << " not found !!! exiting \n";
            exit(0);
         }
        lem.C_select[ix] = index;
     // std::cout << "index for " << C_list[ix]<< "=" << index <<"\n";
    }

    if(lem.C_select.size()<1){
        std::cout << "ERROR, species list for C is empty !!! exiting \n";
        exit(0);
    }
}

/**----------------------------------------------------------------------------
 * @brief       Gets the progress variable for each cell using  mass fraction
 * of species O2
 *
 * @author      A.Menon Mar 26 2022
 *---------------------------------------------------------------------------*/
double LEMLINE::getProgVble_premixO2(double Y_O2)
{

    auto gas = lem.gas;
    int indO2 = gas->speciesIndex("O2");

    double Yu = Yu_p[indO2];
    double Yb = Yb_p[indO2];

    double C = (Y_O2-Yu)/(Yb - Yu);

    C = C < 0? 0:C;
    C = C > 1? 1:C;

    return  C;
}



/**----------------------------------------------------------------------------
 * @brief       Returns the value \int_{0}^{1} P(Z|zMean, zVar)Phi(Z)dz. Uses
 * the scalars conditioned on Z. On the fly PDF calculation, so it can be slow.
 * Returns a vector with  [T, DivU, dedt,Y....]
 *
 *
 * @param      zMean    \input  Favre averaged Mixture fraction from CFD
 * @param      zVar     \input  Favre averaged mixture fraction variance from
 *                              CFD
 *
 * @author      A.Menon, 17th August 2020
 *
 * UPDATE now using scalarZ[Zbin][scalar]
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getPDFMeanSZ(double zMean, double zVar, double press){


    // --- random variable P(Z|zMean, zVar) * phi(Z)
    size_t RSize = nSpecies + NSCALAR;
    std::vector<double> fx(ZBinSize, 0.0);
    std::vector<double> sum(RSize, 0.0);
    double gamma, alpha, beta;


    // --- coefficients and make beta distributioZBinSize -1 object
    gamma = zMean * (1. - zMean)/zVar - 1.;
    alpha = gamma * zMean;
    beta  = gamma * (1. - zMean);


    // -- generate PDF integrals, single call for ibeta
    for(size_t ix = 0; ix < ZBinSize; ix++){
        double Z;
        Z  = (ix) * Zstep;
        try{
            fx[ix] = boost::math::ibeta(alpha, beta, Z);
        }
        catch(boost::exception const& err){
            std::cout << "Boost ERROR for Zmean,Zvar,Z = "
                    << zMean << " "
                    << zVar << " "
                    << Z << "\n";
        }
    }


    // -- sum up
    for(size_t ix = 0; ix < ZBinSize -1; ix++){
        double diff = fx[ix+1] - fx[ix];
        for(size_t jx = 0; jx < RSize; jx++){
            double ag = (scalarZ[ix][jx] + scalarZ[ix+1][jx]);
            sum[jx] += 0.5*ag*diff;
        }
    }

    // -- sanity check
    double T = sum[IT];
    if(T < 100 || T > 4000){
          return  std::vector<double> (nSpecies+ NSCALAR, 0.0);
    }

    // -- fix remaining scalars for filtered state
    auto gas = lem.gas;
    gas->setState_TPY(298.15, 101325, sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
    double rho  = gas->density();

    sum[IRho]   = rho;
    sum[Ih]     = gas->enthalpy_mass()- hForm;

    return sum;         // returns the integrated value
}

/**----------------------------------------------------------------------------
 * @brief       Returns the value \int_{0}^{1} P(C|cMean)Phi(C)dz.
 *              as a vector with  [T, DivU, dedt,Y....] {check scalar
 *              positions}
 *
 *
 * @param      cMean    \input  Mean of progress variable
 * @param      press    \input  pressure
 *
 * @author      A.Menon, 4th Dec 2021
 *
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getPDFMeanSC(double cMean, double press, std::vector<int> pList){


    size_t RSize = nSpecies + NSCALAR;
    size_t pSize = pList.size();
    std::vector<double> fxC(CBinSize, 0.0);
    std::vector<double> sum(RSize+pSize, 0.0);
    std::vector<double>  &H = scalarC_flag; // the 'non-holes' matrix
    auto gas = lem.gas;


    cMean = cMean < 1e-7? 1e-7: cMean; // prevent singularities for now
    cMean = cMean > 1-1e-7? 1-1e-7: cMean; // prevent singularities for now

    // -- generate the C PDF, Alan's step function
    for(size_t ix = 0; ix < CBinSize; ix++){
        double C;
        C  = (ix) * Cstep;

        // integral areas like the incomplete beta function
        if(C < cMean)  fxC[ix] = (1.-cMean)/cMean*C;
        else fxC[ix] = (1-cMean)+cMean/(1.-cMean)*(C-cMean);

    }

    //-- integrate the scalars plus 'non-holes' scaling factor
    double S = 0;

    for(size_t jx = 0; jx < CBinSize-1; jx++){
        double diffC = fxC[jx+1] - fxC[jx];

        S += 0.5*(H[jx]+H[jx+1])*diffC;

        //-- all the usual scalars
        for(size_t kx = 0; kx < RSize;kx++){
            sum[kx] += 0.5*(scalarC[jx][kx]+scalarC[jx+1][kx])*diffC;
        }

        // -- the selected production rates
        for(size_t kx = 0; kx < pSize;kx++){
            size_t sx = pList[kx];
            sum[kx+RSize] += 0.5*(speciesProdC[jx][sx]+speciesProdC[jx+1][sx])*diffC;
        }
    }

    // -- normalize using S
    if(S>0)
    for(size_t kx = 0; kx < RSize+pSize; kx++){
        sum[kx] /= S;
    }

    // -- sanity check
    double T = sum[IT];
    if(T < 100 || T > 4000){
        std::cout << "conditioned T's \n";
        for(size_t ix = 0; ix < CBinSize; ix++){
            std::cout << scalarC[ix][IT]<< " ";
        }
        return  std::vector<double> (nSpecies+ NSCALAR, 0.0);
    }

    // -- fix remaining scalars for filtered state
    gas->setState_TPY(298.15, 101325, sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    gas->setState_TPY(T, press, sum.data()+IY);
    double rho  = gas->density();

    sum[IRho]   = rho;
    sum[Ih]     = gas->enthalpy_mass()- hForm;


    return sum;         // returns the integrated value
}

/**----------------------------------------------------------------------------
 * @brief       Returns the value \int_{0}^{1} P(C|cMeanr)Phi(C)dz.
 *              as a vector with  [T, DivU, dedt,Y....] {check scalar
 *              positions}. Uses Beta PDF instead of the step PDF
 *
 *
 * @param      cMean    \input  Mean of progress variable
 * @param      press    \input  pressure
 *
 * @author      A.Menon, 4th Dec 2021
 *
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getPDFMeanSC_Beta(double cMean, double cVar, std::vector<int> pList){


    size_t RSize = nSpecies + NSCALAR;
    size_t pSize = pList.size();
    std::vector<double> sum(RSize+pSize, 0.0);
    auto gas = lem.gas;

    std::vector<double> fxC=interpBeta(cMean,cVar);

    for(size_t jx = 0; jx < CBinSize-1; jx++){
        double diffC = fxC[jx+1] - fxC[jx];

        //-- all the usual scalars
        for(size_t kx = 0; kx < RSize;kx++){
            sum[kx] += 0.5*(scalarC[jx][kx]+scalarC[jx+1][kx])*diffC;
        }

        // -- the selected production rates
        for(size_t kx = 0; kx < pSize;kx++){
            size_t sx = pList[kx];
            sum[kx+RSize] += 0.5*(speciesProdC[jx][sx]+speciesProdC[jx+1][sx])*diffC;
        }
    }

    // -- sanity check
    double T = sum[IT];

    // -- fix remaining scalars for filtered state
    gas->setState_TPY(298.15, 101325.0, sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
    double rho  = gas->density();

    sum[IRho]   = rho;
    sum[Ih]     = gas->enthalpy_mass()- hForm;

    return sum;         // returns the integrated value
}
/**----------------------------------------------------------------------------
 * @brief       Returns the value \int_{0}^{1} P(C|cMean)Phi(C)dz.
 *              as a vector with  [T, DivU, dedt,Y....] {check scalar
 *              positions}. Uses Beta PDF instead of the step PDF
 *
 *
 * @param      cMean    \input  Mean of progress variable
 * @param      press    \input  pressure
 *
 * @author      A.Menon, 4th Dec 2021
 *
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getPDFMeanSC_TH(double cMean, double cVar, std::vector<int> pList){


    size_t RSize = nSpecies + NSCALAR;
    size_t pSize = pList.size();
    std::vector<double> sum(RSize+pSize, 0.0);
    auto gas = lem.gas;


    std::vector<double> fxC = TopHat2(cMean, cVar, Cstep, CBinSize);

    for(size_t jx = 0; jx < CBinSize-1; jx++){
        double diffC = fxC[jx+3] - fxC[jx+2];

        //-- all the usual scalars
        for(size_t kx = 0; kx < RSize;kx++){
            sum[kx] += 0.5*(scalarC[jx][kx]+scalarC[jx+1][kx])*diffC;
        }

        // -- the selected production rates
        for(size_t kx = 0; kx < pSize;kx++){
            size_t sx = pList[kx];
            sum[kx+RSize] += 0.5*(speciesProdC[jx][sx]+speciesProdC[jx+1][sx])*diffC;
        }
    }

    double w0 = fxC[0];
    double w1 = fxC[1];
    for(size_t kx = 0; kx < RSize+pSize;kx++){
        sum[kx] += (w0*scalarC[0][kx] + w1*scalarC[CBinSize-1][kx]);
    }

    // -- sanity check
    double T = sum[IT];

    // -- fix remaining scalars for filtered state
    gas->setState_TPY(298.15, 101325.0, sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
    double rho  = gas->density();

    sum[IRho]   = rho;
    sum[Ih]     = gas->enthalpy_mass()- hForm;

    return sum;         // returns the integrated value
}

/**----------------------------------------------------------------------------
 * @brief Reads T, P and Y from CFD to give an implicit source term for C mean
 * based on the production rates of some species (O2 in this case)
 *
 *---------------------------------------------------------------------------*/

double LEMLINE::explicitCsource(double *dYdt, double Z){

     if(Z <= Zstep/20 || Z >= ZMax -Zstep/20) {
        return 0.0;
    }

    int Zpos     = std::floor(Z/Zstep);
    int ZposIncr = std::min(ZBinSize-1, Zpos+1);
    double diffZ = Z - Zpos*Zstep;

    double Yu;
    double Yb;
    double Yu1 = 0, Yb1 = 0, dYdt_now = 0;


    for(size_t ix = 0; ix < lem.C_select.size(); ix++){
        size_t cs = lem.C_select[ix];

        Yu= Z*lem.fuelY[cs] + (1-Z)*lem.oxY[cs];

        // interpolate  to find yb, 1st order accuracy
        double m = (lem.YEquilZ[ZposIncr][cs] - lem.YEquilZ[Zpos][cs])/Zstep;
        Yb = lem.YEquilZ[Zpos][cs] + m*diffZ;

        Yu1 += Yu;
        Yb1 += Yb;
        dYdt_now += dYdt[cs];
    }

    dYdt_now/= (Yb1-Yu1);
    return dYdt_now;

}

/**----------------------------------------------------------------------------
 * @brief      Returns the integral phi(Z,C)P(C|Cmean)P(Z|Zmean,Zvar) over the
 * C = [0,1]; Z = [0,1]. Beta PDF for Z and new Modified PDF for C, holes are
 * ignored by implicitly rescaling the PDF for the found values on the line.
 *
 * OBS - This requires scalarZC[][][] to beinitialized to 0 and foundMat[][]
 * to be updated during conditioning.
 *
 * @param      zMean    \input  Favre averaged Mixture fraction from CFD
 * @param      CMean    \input  Favre averaged progress variable from CFD
 * @param      zVar     \input  Favre averaged mixture fraction variance from
 *                              CFD
 *
 *---------------------------------------------------------------------------*/
//std::vector<double> LEMLINE::getPDFMeanS4(double zMean, double zVar, double cMean, double press){
//
//    int RSize = lem.ns + NSCALAR;
//    std::vector<double>  sum(RSize, 0.0);
//    std::vector<double>  fxZ(ZBinSize, 0.0);
//    std::vector<double>  fxZ_d(ZBinSize, 0.0);
//    std::vector<double>  fxC(CBinSize, 0.0);
//    std::vector<double>  sumZ(RSize, 0.0);
//    std::vector<std::vector<double>>  &H = scalarZC_flag; // the 'non-holes' matrix
//    double gamma, alpha, beta;
//
//
//    // --- coefficients and make beta PDF
//    gamma = zMean * (1. - zMean)/zVar - 1.;
//    alpha = gamma * zMean;
//    beta  = gamma * (1. - zMean);
//
//
//    // -- generate PDF integrals, single call for ibeta
//    for(size_t ix = 0; ix < ZBinSize; ix++){
//        double Z;
//		Z = ix* Zstep;
//        fxZ[ix + 0] = boost::math::ibeta(alpha, beta, Z);
//    }
//
//    cMean = cMean == 0.0? 1e-5: cMean; // prevent singularities for now
//    cMean = cMean == 1.0? 1-1e-5: cMean;
//
//    // -- generate the C PDF
//    for(size_t ix = 0; ix < CBinSize; ix++){
//        double C;
//     // C  = ix * 1./CBinSize + 0.5/CBinSize;
//		C = ix* Cstep;
//
//        // integral areas like the incomplete beta function
//        if(C < cMean)  fxC[ix] = (1.-cMean)/cMean*C;
//        else fxC[ix] = (1-cMean)+cMean/(1.-cMean)*(C-cMean);
//
//    }
//
//    //-- integrate the scalars plus 'non-holes' scaling factor
//    double S = 0;
//    double SZ[ZBinSize]{}; // for the C source term integration
//
//    for(size_t ix = 0; ix < ZBinSize-1; ix++){
//        double diffZ = fxZ[ix+1] - fxZ[ix];
//		fxZ_d[ix] = diffZ;
//        for(size_t jx = 0; jx < CBinSize-1; jx++){
//            double diffC = fxC[jx+1] - fxC[jx];
//
//            S += 0.25*(H[ix  ][jx]+H[ix  ][jx+1]+
//                       H[ix+1][jx]+H[ix+1][jx+1])*diffZ*diffC;
//
//            SZ[ix] = H[ix][jx] > 0.0 ? 1.0 : SZ[ix];
//
//            for(size_t kx = 0; kx < RSize;kx++){
//                sum[kx] += 0.25*(scalarZC[ix  ][jx][kx] + scalarZC[ix  ][jx+1][kx] +
//                                 scalarZC[ix+1][jx][kx] + scalarZC[ix+1][jx+1][kx])
//                            *diffC*diffZ;
//            }
//
//        }
//
//        SZ[ix] = H[ix][CBinSize-1] > 0.0? 1.0 :SZ[ix];
//    }
//
//
//    // -- normalize using S
//    for(size_t kx = 0; kx < RSize; kx++){
//        sum[kx] /= S;
//    }
//
//    // --- now integrate C source term only on Z
//    S = 0;
//    sum [IC] = 0;
//    double sumC = 0;
//    for(size_t ix = 0; ix < ZBinSize-1; ix++){
//        double diffZ = fxZ[ix+1] - fxZ[ix];
//        sumC += (scalarZ[ix][IC] + scalarZ[ix+1][IC])*0.5*diffZ;
//        S += diffZ;
//    }
//
//    sum[IC] = sumC;
//    auto gas = lem.gas;
//    double T = sum[IT];
//
//    gas->setState_TPY(298, 101325, sum.data()+IY);
//    double hForm = gas->enthalpy_mass(); //mean formation
//
//    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
//    double rho = gas->density();
//    double hIn = gas->enthalpy_mass(); // Joule/kg
//    hIn -= hForm;
//
//    sum[IRho]  = rho;
//    sum[Ih]    = hIn; // sensible enthalpy
//
//    return sum;         // returns the integrated value
//}

/**----------------------------------------------------------------------------
 * @brief      Returns the integral phi(Z,C)P(C|Cmean)P(Z|Zmean,Zvar) over the
 * C = [0,1]; Z = [0,1]. Beta PDF for Z and delta PDF for C, reworked
 * getPDFMeanS_2 for using delta PDF areas instead of interpolation.
 *
 * OBS- required matrix intialization
 *
 *
 * @param      zMean    \input  Favre averaged Mixture fraction from CFD
 * @param      CMean    \input  Favre averaged progress variable from CFD
 * @param      zVar     \input  Favre averaged mixture fraction variance from
 *                              CFD
 *
 * UPDATE:: Too slow!! too many 'if' conditions in inner loop
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getPDFMeanS5(double zMean, double zVar, double cMean, double press){

    int RSize = nSpecies + NSCALAR;
    std::vector<double>  sum(RSize, 0.0);
    std::vector<double>  fxZ(ZBinSize, 0.0);
    std::vector<double>  fxC(CBinSize, 0.0);
    double gamma, alpha, beta;


    // --- coefficients and make beta PDF
  //gamma = zMean * (1. - zMean)/zVar - 1.;
  //alpha = gamma * zMean;
  //beta  = gamma * (1. - zMean);


  //for(size_t ix = 0; ix < ZBinSize; ix+=2){
  //    double Z;
  //    Z  = (ix) * Zstep;
  //    fxZ[ix + 0] = boost::math::ibeta(alpha, beta, Z);
  //}



    // -- generate the C PDF: rectangle approximation of Dirac funtion PDF
    double S=0;
    size_t cStart = std::floor(cMean/Cstep);

    for(size_t ix = cStart; ix < CBinSize; ix++){
        double C;
        fxC[ix] = 1.0;
        S += 1.0;
    }

    if(S <=0) {std::cout << "C fx error, exiting\n"; exit(0);}
    S = 0;

    for(size_t ix = 0; ix < ZBinSize-1; ix++){
        double diffZ = fxZ[ix+1] - fxZ[ix];
        for(size_t jx = 0; jx < CBinSize-1; jx++){
            double diffC = fxC[jx+1] - fxC[jx];

            S += diffC*diffZ;

            for(size_t kx = 0; kx < RSize;kx++){
                sum[kx] += 0.25*(scalarZC[ix  ][jx][kx] + scalarZC[ix  ][jx+1][kx] +
                                 scalarZC[ix+1][jx][kx] + scalarZC[ix+1][jx+1][kx])
                            *diffC*diffZ;
            }

        }
    }

    // --- now integrate C source term only on Z
    double sumC = 0;
    for(size_t ix = 0; ix < ZBinSize-1; ix++){
        double diffZ = fxZ[ix+1] - fxZ[ix];
        sumC += (scalarZ[ix][IC] + scalarZ[ix+1][IC])*0.5*diffZ;
    }

    sum[IC] = sumC;

//  if(sum[IC]>1.0){

//      std::cout<<"PDF ERROR!!\n\n"
//               <<" PDF sum = " << S <<"\n"
//               <<" C_PDF sum = " << sum[IC] <<"\n"
//               << "C(Z) = \n";
//      for(size_t ix = 0; ix < ZBinSize; ix++)
//          std::cout << scalarZ[ix][IC] << "\t"<<fxZ[ix+1] - fxZ[ix]<<"\n";

//      exit(0);
//  }


    auto gas = lem.gas;

    double T = sum[IT];
//  if(T < 100 || T > 4000){
//      std::cout<< "PDF integration error\n";
//      exit(0);
//  }

    gas->setState_TPY(298, 101325, sum.data()+IY);
    double hForm = gas->enthalpy_mass(); //mean formation

    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
    double rho = gas->density();
    double hIn = gas->enthalpy_mass(); // Joule/kg
    hIn -= hForm;

    sum[IRho]  = rho;
    sum[Ih]    = hIn; // sensible enthalpy

    return sum;         // returns the integrated value
}


/**
 *@brief        returns interplated T and mass fractions, density and sensible enthalpy set after
 *
 * \param       zMean
 * \param       cMean
 * \param       press   pressure
 * \param       uFlag   unburnt flag
 *
 *@author       A. Menon
 */

std::vector<double> LEMLINE::getInterpS(double zMean,double cMean,double press, bool uFlag) {
    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize, 0.0);
    auto gas = lem.gas;

    if(uFlag){
        for(size_t ix = 0; ix < nSpecies;ix++){
            sum[ix + IY] = zMean*lem.fuelY[ix] + (1. - zMean)*lem.oxY[ix];
        }
        sum[IT] = zMean*lem.fuelStreamT + (1. - zMean)*lem.oxStreamT;
        return sum;
    }

    int Zpos = std::floor(zMean/Zstep);
    int Cpos = std::floor(cMean/Cstep);

        Zpos = std::min(ZBinSize-1, Zpos);
        Cpos = std::min(CBinSize-1, Cpos);

    int ZposIncr = std::min(ZBinSize-1, Zpos+1);
    int CposIncr = std::min(CBinSize-1, Cpos+1);

    // -- interpolate only filtered state, set the rest with cantera
    double diffZ = zMean - Zpos*Zstep;
    double diffC = cMean - Cpos*Cstep;


    for(size_t ix = 0; ix < RSize;ix++){
        // interpolate along Z for Cpos
        double m = (scalarZC[ZposIncr][Cpos][ix] - scalarZC[Zpos][Cpos][ix])/Zstep;
        double sum1 = scalarZC[Zpos][Cpos][ix] + m*diffZ;

        // interpolate along Z for Cpos+1
        m = (scalarZC[ZposIncr][CposIncr][ix] - scalarZC[Zpos][CposIncr][ix])/Zstep;
        double sum2 = scalarZC[Zpos][CposIncr][ix] + m*diffZ;
        // interpolate along C
        //
        m = (sum2-sum1)/Cstep;
        sum[ix] = sum1 + m*diffC;
    }

    // -- formation enthalpy
    gas->setState_TPY(298.15,101325,sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    // -- sensible enthalpy
    gas->setState_TPY(sum[IT],lem.pressure, sum.data()+IY);
    double hs = gas->enthalpy_mass() -hForm;

    sum[IRho] = gas->density();
    sum[Ih]   = hs;
    return sum;

}

/**---------------------------------------------------------------------------
 * @brief        returns interplated scalars, Interpolate only in Z (should not
 * change elemental mass fractions), choose to nearest filled Cbins for given Z
 *
 * \param       zMean
 * \param       cMean
 *---------------------------------------------------------------------------*/

std::vector<double> LEMLINE::getInterpS_fillC(double zMean,double cMean) {
    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize, 0.0);
    auto gas = lem.gas;

    int Zpos = std::floor(zMean/Zstep);
    int ZposIncr = std::min(ZBinSize-1, Zpos+1);

    int Cpos1_start = std::floor(cMean/Cstep);
    int Cpos2_start = std::floor(cMean/Cstep);

    int up, down;
    bool upF =false, downF =false;

    for (up = Cpos1_start; up < CBinSize; up++){
        if (scalarZC_flag[Zpos][up]) {upF = true; break;}
    }

    for (down = Cpos1_start; down > 0; down--){
        if (scalarZC_flag[Zpos][down] ) {downF = true; break;}
    }
    int upDist = up - Cpos1_start;
    int downDist = Cpos1_start-down;

    if(!upF && !downF) Cpos1_start = 0;
    else if(!downF) Cpos1_start = up;
    else if(!upF) Cpos1_start = down;
    else Cpos1_start = (upDist < downDist)? up: down;

    upF =false; downF =false;
    for (up = Cpos2_start; up < CBinSize; up++){
        if (scalarZC_flag[ZposIncr][up]) {upF = true; break;}
    }

    for (down = Cpos2_start; down > 0; down--){
        if (scalarZC_flag[ZposIncr][down]) {downF = true; break;}
    }
    upDist = up - Cpos1_start;
    downDist = Cpos1_start-down;

    if(!upF && !downF) Cpos2_start = 0;
    else if(!downF) Cpos2_start = up;
    else if(!upF) Cpos2_start = down;
    else Cpos2_start = (upDist < downDist)? up: down;



    // -- interpolate onl in Z, with Cpos2_start and Cpos1_start
    double diffZ = zMean - Zpos*Zstep;

    for(size_t ix = 0; ix < RSize;ix++){
        // interpolate along Z for Cpos
        double m = (scalarZC[ZposIncr][Cpos2_start][ix] - scalarZC[Zpos][Cpos1_start][ix])/Zstep;
        double sum1 = scalarZC[Zpos][Cpos1_start][ix] + m*diffZ;
        sum[ix] = sum1;
    }

    // -- formation enthalpy
    gas->setState_TPY(298.15,101325,sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    // -- sensible enthalpy
    gas->setState_TPY(sum[IT],lem.pressure, sum.data()+IY);
    double hs = gas->enthalpy_mass() -hForm;

    sum[IRho] = gas->density();
    sum[Ih]   = hs;
    sum[ITheta] = getdCdt_fillC(zMean, cMean);

    return sum;

}

/**---------------------------------------------------------------------------
 * @brief        returns interplated scalars, Interpolate only in Z (should not
 * change elemental mass fractions), choose to nearest filled Cbins for given Z
 *
 * \param       zMean
 * \param       cMean
 *---------------------------------------------------------------------------*/

double LEMLINE::getdCdt_fillC(double zMean,double cMean) {
    double  sum = 0;
    auto gas = lem.gas;

    int Zpos = std::floor(zMean/Zstep);
    int ZposIncr = std::min(ZBinSize-1, Zpos+1);

    int Cpos1_start = std::floor(cMean/Cstep);
    int Cpos2_start = std::floor(cMean/Cstep);

    int up, down;
    bool upF =false, downF =false;

    for (up = Cpos1_start; up < CBinSize; up++){
        if (dCdt_flag[Zpos][up]) {upF = true; break;}
    }

    for (down = Cpos1_start; down > 0; down--){
        if (dCdt_flag[Zpos][down] ) {downF = true; break;}
    }
    int upDist = up - Cpos1_start;
    int downDist = Cpos1_start-down;

    if(!upF && !downF) Cpos1_start = 0;
    else if(!downF) Cpos1_start = up;
    else if(!upF) Cpos1_start = down;
    else Cpos1_start = (upDist < downDist)? up: down;

    upF =false; downF =false;
    for (up = Cpos2_start; up < CBinSize; up++){
        if (dCdt_flag[ZposIncr][up]) {upF = true; break;}
    }

    for (down = Cpos2_start; down > 0; down--){
        if (dCdt_flag[ZposIncr][down]) {downF = true; break;}
    }
    upDist = up - Cpos1_start;
    downDist = Cpos1_start-down;

    if(!upF && !downF) Cpos2_start = 0;
    else if(!downF) Cpos2_start = up;
    else if(!upF) Cpos2_start = down;
    else Cpos2_start = (upDist < downDist)? up: down;



    // -- interpolate onl in Z, with Cpos2_start and Cpos1_start
    double diffZ = zMean - Zpos*Zstep;
        // interpolate along Z for Cpos
    double m = (scalarZC[ZposIncr][Cpos2_start][ITheta] - scalarZC[Zpos][Cpos1_start][ITheta])/Zstep;
    sum = scalarZC[Zpos][Cpos1_start][ITheta] + m*diffZ;

    return sum;

}

/**---------------------------------------------------------------------------
 * @brief       For each Z bin, scan upward and fill missing values with
 * previous filled (in C space) values. Does not update scalarZC_flag. For use
 * with PDFs
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::fillZC_biased()
{
    for(size_t ix = 0 ; ix < ZBinSize; ix++){
        for(size_t jx = 1; jx < CBinSize; jx ++){
            if(!scalarZC_flag[ix][jx]){
            for(size_t kx = 0; kx < NSCALAR + nSpecies; kx ++){
                scalarZC[ix][jx][kx] = scalarZC[ix][jx-1][kx];
            }
            }
          //if(!dCdt_flag[ix][jx]){
          //    scalarZC[ix][jx][ITheta] = scalarZC[ix][jx-1][ITheta];
          //}
        }
    }
}




/**---------------------------------------------------------------------------
 *@brief        returns interplated T and mass fractions, density and sensible
                enthalpy set after, only C conditioning
 *
 * \param       zMean
 * \param       press   pressure
 * \param       uFlag   unburnt flag
 *
 *@author       A. Menon
 *---------------------------------------------------------------------------*/

std::vector<double> LEMLINE::getInterpSC(double cMean,double CFDpress, std::vector<int> pList) {
    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize+pList.size(), 0.0);
    auto gas = lem.gas;

    int Cpos = std::floor(cMean/Cstep);
    Cpos = Cpos < 0? 0 :Cpos;
    Cpos = std::min(CBinSize-1, Cpos);
    int CposIncr = std::min(CBinSize-1, Cpos+1);

    double diffC = cMean - Cpos*Cstep;

    for(size_t ix = 0; ix < RSize;ix++){
        // interpolate along C for Cpos
        double m = (scalarC[CposIncr][ix] - scalarC[Cpos][ix])/Cstep;
        sum[ix]  = scalarC[Cpos][ix] + m*diffC;
    }

    //-- interpolate production rates
    for(size_t ix = 0; ix < pList.size();ix++){
        size_t pPos = pList[ix];
        double m = (speciesProdC[CposIncr][pPos] - speciesProdC[Cpos][pPos])/Cstep;
        sum[RSize+ix] = speciesProdC[Cpos][pPos] + m*diffC;
    }

    // -- formation enthalpy
    gas->setState_TPY(298.15,101325,sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    // -- sensible enthalpy
    gas->setState_TPY(sum[IT],CFDpress, sum.data()+IY);
    double hs = gas->enthalpy_mass() -hForm;

    sum[IRho] = gas->density();
    sum[Ih]   = hs;
    return sum;

}

/**---------------------------------------------------------------------------
 *@brief        returns interpolated T and mass fractions, density and sensible
                enthalpy set after, only C conditioning when scalarC or
                speciesProdC are initialized as zeros.
                Essentially reduces to a delta pdf
 *
 * \param       cMean \input favre averaged C mean
 * \param       pList
 *
 *@author       A. Menon
 *---------------------------------------------------------------------------*/

std::vector<double> LEMLINE::getInterpSC_incomplete(double cMean, std::vector<int> pList) {
    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize+pList.size(), 0.0);
    auto gas = lem.gas;

    int Cpos1 = cMean/(Cstep);
    int Cpos2 = Cpos1+1;

    while((scalarC_flag[Cpos1]!= 1) && (Cpos1 > -1)){
        Cpos1--;
    }
    while(scalarC_flag[Cpos2]!= 1 && Cpos2 < CBinSize){
        Cpos2++;
    }

    // --
    if(Cpos1 >= 0){
        double diffC = cMean - Cpos1*Cstep;
        if(Cpos2 < CBinSize){
            // mass fractions and T
            for(size_t ix = 0; ix < RSize;ix++){
                // interpolate along C for Cpos
                double m = (scalarC[Cpos2][ix] - scalarC[Cpos1][ix])/Cstep/(Cpos2-Cpos1);
                sum[ix] = scalarC[Cpos1][ix] + m*diffC;
            }

            //-- interpolate production rates
            for(size_t ix = 0; ix < pList.size();ix++){
                // interpolate along Z for Cpos
                size_t pPos = pList[ix];
                double m = (speciesProdC[Cpos2][pPos] - speciesProdC[Cpos1][pPos]);
                       m = m/(Cpos2-Cpos1)/Cstep;
                sum[RSize+ix] = speciesProdC[Cpos1][pPos] + m*diffC;
            }
        }
        else{ // run out on burnt side

            for(size_t ix = 0; ix < RSize;ix++)
                sum[ix] = scalarC[Cpos1][ix];

            for(size_t ix = 0; ix < pList.size();ix++)
                sum[RSize+ix] = speciesProdC[Cpos1][pList[ix]];
        }
    }
    else{
        // run out on un-burnt side

        for(size_t ix = 0; ix < RSize;ix++)
            sum[ix] = scalarC[Cpos2][ix];

        for(size_t ix = 0; ix < pList.size();ix++)
            sum[RSize+ix] = speciesProdC[Cpos2][pList[ix]];
    }

    // -- formation enthalpy
    gas->setState_TPY(298.15,101325,sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    // -- sensible enthalpy
    gas->setState_TPY(sum[IT],lem.pressure, sum.data()+IY);
    double hs = gas->enthalpy_mass() -hForm;

    sum[IRho] = gas->density();
    sum[Ih]   = hs;
    return sum;

}

/**---------------------------------------------------------------------------
 *@brief        returns interplated T and mass fractions, density and sensible
                enthalpy set after, only Z conditioning
 *
 * \param       zMean
 * \param       press   pressure
 * \param       uFlag   unburnt flag
 *
 *@author       A. Menon
 *---------------------------------------------------------------------------*/

std::vector<double> LEMLINE::getInterpSZ(double zMean,double press) {
    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize, 0.0);
    auto gas = lem.gas;

    size_t Zpos = zMean*(ZBinSize-2)/ZMax;

    // -- interpolate only filtered state, set the rest with cantera
    double diffZ = zMean - Zpos*Zstep;

    for(size_t ix = 0; ix < RSize;ix++){
        // interpolate along Z for Cpos
        double m = (scalarZ[Zpos+1][ix] - scalarZ[Zpos][ix])/Zstep;
        sum[ix] = scalarZ[Zpos][ix] + m*diffZ;
    }

    // -- formation enthalpy
    gas->setState_TPY(298.15,101325,sum.data()+IY);
    double hForm = gas->enthalpy_mass();
    // -- sensible enthalpy
    gas->setState_TPY(sum[IT],lem.pressure, sum.data()+IY);
    double hs = gas->enthalpy_mass() -hForm;

    sum[IRho] = gas->density();
    sum[Ih]   = hs;
    return sum;

}

/**----------------------------------------------------------------------------
 * @brief      Returns the integral phi(Z,C)P(C|Cmean)P(Z|Zmean,Zvar) over the
 * C = [0,1]; Z = [0,1]. P(C) assumes Dirac delta funtion P(Z) is assumed to be
 * a beta PDF
 *
 * @param      zMean    \input  Favre averaged Mixture fraction from CFD
 * @param      CMean    \input  Favre averaged progress variable from CFD
 * @param      zVar     \input  Favre averaged mixture fraction variance from
 *                              CFD
 * @author      A.Menon 19th Nov
 *
 * @author      A.Menon 8th April 2021
 *
 * EDIT: tStep no longer functional, perhaps merge this with earlier versions
 * or remove entirely
 *
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getPDFMeanS2(double zMean, double zVar, double cMean, double press){


    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize, 0.0);
    std::vector<double>  fx(ZBinSize, 0.0);
    double phiFx[ZBinSize][RSize];
    double gamma, alpha, beta;


    int cPos = cMean/Cstep;
    cPos = std::min(cPos,CBinSize-1);
    int cPosIncr = std::min(cPos+1,CBinSize-1);

    double m = cMean - cPos*Cstep;

    for(size_t ix = 0; ix < ZBinSize; ix++){
        double *S1 = scalarZC[ix][cPos].origin();
        double *S2 = scalarZC[ix][cPosIncr].origin();
        for(size_t jx = 0; jx < RSize; jx++){
            phiFx[ix][jx] = S1[jx] + (S2[jx]-S1[jx])*m;// y = {m}x + C
        }
    }

    // --- coefficients and make beta PDF
    gamma = zMean * (1. - zMean)/zVar - 1.;
    alpha = gamma * zMean;
    beta  = gamma * (1. - zMean);


    // -- generate PDF integrals, single call for ibeta
    for(size_t ix = 0; ix < ZBinSize; ix++){
        double Z;
        Z  = (ix) * Zstep;
        fx[ix + 0] = boost::math::ibeta(alpha, beta, Z);
    }


    // -- sum up
    for(size_t ix = 0; ix < ZBinSize -1; ix++){
        double diff = fx[ix+1] - fx[ix];
        for(size_t jx = 0; jx < RSize; jx++){
            double ag = 0.5 * (phiFx[ix][jx] + phiFx[ix+1][jx]);
            sum[jx] += ag*diff;
        }
    }


    auto gas = lem.gas;
    double T = sum[IT];
    gas->setState_TPY(298, 101325, sum.data()+IY);
    double hForm = gas->enthalpy_mass(); //mean formation

    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
    double rho = gas->density();
    double hIn = gas->enthalpy_mass(); // Joule/kg
    hIn -= hForm;

    sum[IRho]  = rho;
    sum[Ih]    = hIn; // sensible enthalpy

    return sum;         // returns the integrated value
}

/**
 * Same as above but using the top-hat2 PDF for Z
 *
 */
std::vector<double> LEMLINE::getPDFMeanS2_TH(double zMean, double zVar, double cMean){


    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize, 0.0);
    double phiFx[ZBinSize][RSize];
    double gamma, alpha, beta;


    int cPos = cMean/Cstep;
    cPos = std::min(cPos,CBinSize-1);
    int cPosIncr = std::min(cPos+1,CBinSize-1);

    double m = cMean - cPos*Cstep;

    for(size_t ix = 0; ix < ZBinSize; ix++){
        double *S1 = scalarZC[ix][cPos].origin();
        double *S2 = scalarZC[ix][cPosIncr].origin();
        for(size_t jx = 0; jx < RSize; jx++){
            phiFx[ix][jx] = S1[jx] + (S2[jx]-S1[jx])*m;// y = {m}x + C
        }
    }

    std::vector<double> fx = TopHat2(zMean, zVar, Zstep, ZBinSize);

    // -- sum up
    for(size_t ix = 0; ix < ZBinSize -1; ix++){
        double diff = fx[ix+3] - fx[ix+2];
        for(size_t jx = 0; jx < RSize; jx++){
            double ag = 0.5 * (phiFx[ix][jx] + phiFx[ix+1][jx]);
            sum[jx] += ag*diff;
        }
    }

    double w0 = fx[0];
    double w1 = fx[1];
    for(size_t jx = 0; jx < RSize; jx++){
        sum[jx] += w0*phiFx[0][jx] + w1*phiFx[ZBinSize-1][jx];
    }

    auto gas = lem.gas;
    double T = sum[IT];
    gas->setState_TPY(298, 101325, sum.data()+IY);
    double hForm = gas->enthalpy_mass(); //mean formation

    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
    double rho = gas->density();
    double hIn = gas->enthalpy_mass(); // Joule/kg
    hIn -= hForm;

    sum[IRho]  = rho;
    sum[Ih]    = hIn; // sensible enthalpy

    return sum;         // returns the integrated value
}

/**----------------------------------------------------------------------------
 * @brief      Returns the integral phi(Z,C)P(C|Cmean)P(Z|Zmean,Zvar) over the
 * C = [0,1]; Z = [0,1]. Bivariate Beta PDFs, if the variance is too small,
 * switch to the step fuction 'Heavyside'
 *
 * @param      zMean    \input  Favre averaged Mixture fraction from CFD
 * @param      CMean    \input  Favre averaged progress variable from CFD
 * @param      zVar     \input  Favre averaged mixture fraction variance from
 * @param      cVar     \input  Favre averaged mixture fraction variance from
 *                              CFD
 *
 * @author      A.Menon 30 Nov 2023
 *
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getPDFMeanSBV(double zMean, double zVar, double cMean,double cVar){


    int RSize = nSpecies + NSCALAR;
    std::vector<double> sum(RSize, 0.0);
    std::vector<double>  fxZ(ZBinSize, 0.0);
    std::vector<double>  fxC(CBinSize, 0.0);
    double gamma, alpha, beta;
    double gammaC, alphaC, betaC;

    double Zvm = zMean*(1.- zMean);
    double Zmin = 1./10*std::min(zMean*zMean, (1.-zMean)*(1.- zMean));


    // --regulalrized CDF
 // if(zVar <= Zmin || zVar < 1e-5){
    if(zVar < 1e-5){
        fxZ = TopHat(zMean,zVar, Zstep, ZBinSize);
    }
    else{
        // --- Beta PDFs for Z
        zVar = zVar >= Zvm? Zvm*0.9999:zVar;
        gamma = zMean * (1. - zMean)/zVar - 1.;
        alpha = gamma * zMean;
        beta  = gamma * (1. - zMean);
        for(int ix = 0; ix < ZBinSize; ix++){
            double Z;
            Z  = ix*Zstep;
            fxZ[ix] = boost::math::ibeta(alpha, beta, Z);
        }
    }

    double Cvm = cMean*(1.- cMean);
    double Cmin = 1./10*std::min(cMean*cMean, (1.-cMean)*(1.- cMean));

 // if(cVar <= Cmin|| cVar < 1e-5){
    if(cVar < 1e-5){
        fxC = TopHat(cMean,cVar,Cstep, CBinSize);
    }
    else{
    // --- Beta PDFs for C
        cVar = cVar >= Cvm?  Cvm*0.99999: cVar;
        gammaC = cMean * (1. - cMean)/cVar - 1.;
        alphaC = gammaC * cMean;
        betaC  = gammaC * (1. - cMean);
        for(int ix = 0; ix < CBinSize; ix++){
            double C;
            C  = ix*Cstep;
            fxC[ix] = boost::math::ibeta(alphaC, betaC, C);
        }
    }



    // -- Integrate
    for(size_t ix = 0; ix < ZBinSize-1; ix++){
        double diffZ = 0.5*(fxZ[ix+1] - fxZ[ix]);
        for(size_t jx = 0; jx < CBinSize-1; jx++){
            double diffC = 0.5*(fxC[jx+1] - fxC[jx]);
            for(size_t kx = 0; kx < RSize;kx++){
                double s1 = diffC*(scalarZC[ix][jx][kx] + scalarZC[ix][jx+1][kx]);
                double s2 = diffC*(scalarZC[ix+1][jx][kx] + scalarZC[ix+1][jx+1][kx]);

                sum[kx] += diffZ*(s1 + s2);

              //sum[kx] += 0.25*(scalarZC[ix  ][jx][kx] + scalarZC[ix  ][jx+1][kx] +
              //                 scalarzc[ix+1][jx][kx] + scalarZC[ix+1][jx+1][kx])
              //            *diffC*diffZ;

            }
        }
    }


    auto gas = lem.gas;
    double T = sum[IT];
    gas->setState_TPY(298, 101325, sum.data()+IY);
    double hForm = gas->enthalpy_mass(); //mean formation

    gas->setState_TPY(T, lem.pressure, sum.data()+IY);
    double rho = gas->density();
    double hIn = gas->enthalpy_mass(); // Joule/kg
    hIn -= hForm;

    sum[IRho]  = rho;
    sum[Ih]    = hIn;

    return sum;
}

/**----------------------------------------------------------------------------
 * @brief      returns averaged (density weighted) scalars on the line as a vector
 *
 * @author      A.Menon 11th Aug 2021
 *
 *---------------------------------------------------------------------------*/
std::vector<double> LEMLINE::getMeanS(std::vector<int>pList){

    size_t pSize = pList.size();
    int    RSize = nSpecies + NSCALAR;
    double* RR   = lem.tmp_s1;
    double* Ms   = lem.Ms;  //  trying Konduri's units

    std::vector<double> sum(RSize + pSize, 0.0);
    double sumDensity = 0.0;
    double nLEM       = 0;
    auto gas= lem.gas;

    if(cells.size()<1){
        std::cout << "\n\tNO CELLS IN LEM LINE!! Exiting";
        exit(0);
    }

    // -- scan the line
    std::list<CellData>::iterator it;
    for (it = cells.begin(); it != cells.end(); ++it){

        double rhoLoc;

        rhoLoc       = (*it).rho;
     // rhoLoc       = 1.0;        // volume weighted might be more accurate?
        sum[IT]     += (*it).T   *rhoLoc;
        sum[IRho]   += (*it).rho *rhoLoc;
        sum[ITheta] += (*it).dCdt*rhoLoc;
        sum[IC]     += (*it).C   *rhoLoc;

        // find sensible enthalpy
      //gas->setState_TPY(298.15, lem.pressure, (*it).Y);
      //E = (*it).E - gas->enthalpy_mass();

      //sum[Ih]     += E*rhoLoc;

        // -- find production rates for averaging select ones
        gas->setState_TPY((*it).T, lem.pressure, (*it).Y);
        lem.kin->getNetProductionRates(RR);// units Kmol/m^3/s

        for(size_t jx = 0; jx < nSpecies; jx++){
            sum[IY + jx] += (*it).Y[jx] *rhoLoc;
        }

        for(size_t kx = 0; kx < pSize;kx++){
            size_t sx = pList[kx];
            sum[kx+RSize] += RR[sx]*Ms[sx];
        }

        nLEM        +=1.0;
        sumDensity  += rhoLoc;
    }

    for(size_t ix =0; ix < RSize;ix++){
        sum[ix] /= sumDensity;
    }

    for(size_t ix =RSize-1; ix < RSize+pSize;ix++){
        sum[ix] /= nLEM;
    }

    gas->setState_TPY(sum[IT], lem.pressure, sum.data()+IY);
    double E = gas->enthalpy_mass();
    gas->setState_TPY(289.15, 101325, sum.data()+IY);
    E -= gas->enthalpy_mass();

    sum[Ih] = E;


    return sum;
}

/**----------------------------------------------------------------------------
 * @brief      returns averaged (density weighted) scalars on the line as a vector
 *
 * @author      A.Menon 11th Aug 2021
 *
 *---------------------------------------------------------------------------*/
double LEMLINE::getCvar(){

    double mean = getCMean();
    // Calculate the variance.
    double variance = 0.0;
    for (const CellData& cell1 : cells) {
        variance += (cell1.C - mean) * (cell1.C - mean);
    }
    variance /= cells.size();
    return variance;
}

/**----------------------------------------------------------------------------
 * @author      A.Menon 12 Feb, 2024
 *
 *---------------------------------------------------------------------------*/
double LEMLINE::getCMean(){

    // Calculate the mean (average) of the elements in the list.
    double sum = 0.0;
    for (const CellData& cell1 : cells) {
        sum += cell1.C;
    }
    double mean = sum / cells.size();
    return mean;
}



/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       calculates heat release for each cell, reused some code from
 *              getMeanHeatProductionRate()
 *
 *              Reworked to correct units problem
 * @author      A.Menon
 *---------------------------------------------------------------------------*/

void LEMLINE::setHeatRelease(){
    std::vector<double>& hFormation = lem.hFormation;
    auto gas = lem.gas;
    std::list<CellData>::iterator it;
    double* RR   = lem.tmp_s1;
    double* H_RT = lem.tmp_s2;
    double* Ms   = lem.Ms;

    for(it = cells.begin(); it !=cells.end();it++){
        double hSource = 0.0;
        double T = (*it).T;
        gas->setState_TPY((*it).T, (*it).p, (*it).Y);
        lem.kin->getNetProductionRates(RR);// units Kmol/m^3/s
     // gas->getEnthalpy_RT(H_RT);     // non dimentional
        gas->getIntEnergy_RT(H_RT);     // non dimentional
        double  rho = gas->density();

        for(size_t kx = 0; kx < nSpecies; kx++){

            // Units 0 * J/Kmol/K * K = J/kmol.;
            double h_s = H_RT[kx] * T * Cantera::GasConstant; // Joule/kmol
            h_s /= Ms[kx]; // Joule/kg

            double dYdt = RR[kx]*Ms[kx]/rho;
          // kmol/m^3/s * Kg/kmol * m^3/kg = s^-1
            hSource += -h_s * dYdt; // J/kg/s

         // hSource += -hFormation[kx] * dYdt; // J/kg/s

        }

        (*it).dedt = hSource;

    }
}

/**----------------------------------------------------------------------------
 *
 *---------------------------------------------------------------------------*/

double LEMLINE::returnHsource(double *dYdt){

    double hSource = 0.0;
    std::vector<double>& hFormation = lem.hFormation;

    for(size_t kx = 0; kx < nSpecies; kx++){

      //hSource += -hFormation[kx]* //[   J/kg      ]
      //            RR[kx]*         //[kmol/m^3 /s  ]
      //            Ms[kx];         //[  kg/kmol    ]
                                    //[   J/m^3 /s  ]


       hSource += -hFormation[kx]* //[   J/kg      ]
                   dYdt[kx];       //[   1/s       ]
                                   //[   J/kg/s    ]

    }

    return hSource;

}

/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       calculates enthalpy for each cell,
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/

void LEMLINE::setEnergy(){
    auto gas = lem.gas;
    std::list<CellData>::iterator it;

    Energy = 0;
    for(it = cells.begin(); it !=cells.end();it++){
        gas->setState_TPY((*it).T, (*it).p, (*it).Y);
        double  E = gas->intEnergy_mass()*(*it).m;
        (*it).E = E;
        Energy += E;
    }
}
/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       Calculates the divU for each cell, considered here are terms
 * from the chemical heat release, temperature diffusion, enthalpy diffusion
 * and Fickian diffusion
 *
 * See equation (16) from Hcci notes by Prof. M.Oevermann
 *
 * @author      A.Menon
 *---------------------------------------------------------------------------*/

//void LEMLINE::setDivU(){
//
//    double   RR[nSpecies];
//    double H_RT[nSpecies];
//
//    int nCells = cells.size();
//
//    auto gas = lem.gas;
//
//    std::list<CellData>::iterator it;
//
//
//    // scalar arrays for gradient calculation
//
//    std::vector<double> T       (nCells, 0.0);
//    std::vector<double> dx      (nCells, 0.0);
//    std::vector<double> rho     (nCells, 0.0);
//    std::vector<double> M       (nCells, 0.0);
//    std::vector<double> cp      (nCells, 0.0);
//    std::vector<double> lambda  (nCells, 0.0);
//    std::vector<double> Q       (nCells, 0.0);
//
//    double Y [nSpecies][nCells];
//    double hs[nSpecies][nCells];
//
//
//    // scan line and calculate divU1[]
//    size_t ix = 0;
//    std::vector<double> divU1(nCells, 0.0); // heat release contribution
//    for(it = cells.begin(); it !=cells.end();it++){
//
//        // -- scan
//             T[ix] = (*it).T;
//           rho[ix] = (*it).rho;
//            dx[ix] = (*it).dx;
//        lambda[ix] = (*it).lambda;
//
//
//        gas->setState_TPY(T[ix], (*it).p, (*it).Y);
//        lem.kin->getNetProductionRates(RR);
//        gas->getEnthalpy_RT(H_RT);
//
//        // -- part of scan
//         M[ix] = gas->meanMolecularWeight();
//        cp[ix] = gas->cp_mass();
//
//        double divU1_temp = 0.0;
//        for(size_t kx = 0; kx < nSpecies; kx++){
//
//            //--- omega_dot, species mass production rate
//            double speciesProd = lem.Ms[kx] * RR[kx];
//
//            //--- species enthalpy
//            double hs_temp   = H_RT[kx]  * T[ix] * lem.Rs[kx];
//
//            divU1_temp +=
//                (M[ix]/lem.Ms[kx] - hs_temp/cp[ix]/T[ix]) *speciesProd;
//
//            // -- part of scan
//            hs[kx][ix] = hs_temp;
//             Y[kx][ix] = (*it).Y[kx];
//
//        }
//        divU1[ix] = divU1_temp/rho[ix];
//        ix++;
//    }
//
//    // - compute gradients
//    std::vector<double> dTdx = gradient(T,dx);
//
//    double  dYdx[nSpecies][nCells];
//    double dhsdx[nSpecies][nCells];
//
//    for(size_t kx = 0; kx < nSpecies; kx++){
//
//        std::vector<double> Ysx(nCells);
//        // --- copy species kx along line
//        memcpy(Ysx.data(),Y[kx],nCells*sizeof(double));
//
//        // gradient for species kx
//        std::vector<double> dYdxRow = gradient(Ysx, dx);
//
//        // --- copy-store to matrix
//        memcpy(dYdx[kx],dYdxRow.data(), nCells*sizeof(double));
//
//        // --- hs gradients, reusing variables
//        memcpy(Ysx.data(), hs[kx],nCells*sizeof(double));
//        dYdxRow = gradient(Ysx, dx);
//        memcpy(dhsdx[kx], dYdxRow.data(), nCells*sizeof(double));
//    }
//
//    //-- heat conduction
//    for(ix = 0; ix < nCells; ix++){
//        Q[ix]  = -lambda[ix] * dTdx[ix];
//    }
//
//    double js[nSpecies][nCells];        // species diffusion
//    for(size_t kx = 0; kx < nSpecies; kx++){
//        for(ix = 0; ix < nCells; ix++){
//            // ---- Assuming unity Lewis number
//        //  double Diff = lambda[ix]/rho[ix]/cp[ix];
//            double Diff = lambda[ix]        /cp[ix];
//            // ---- diffusive velocity, uncorrected for now
//            double vs = -Diff * dYdx[kx][ix];
//       //   js[kx][ix] =  vs * rho[ix] * Y[kx][ix];
//            js[kx][ix] =  vs           * Y[kx][ix];
//        }
//    }
//
//    // -- Gradients for Q and js
//    std::vector<double> dQdx = gradient(Q, dx);
//    double djsdx[nSpecies][nCells];
//
//    // -- same as before, reusing variable names
//    for(size_t kx = 0; kx < nSpecies; kx++){
//
//        std::vector<double> Ysx(nCells);
//        memcpy(Ysx.data(),js[kx],nCells*sizeof(double));
//        std::vector<double> dYdxRow = gradient(Ysx, dx);
//        memcpy(djsdx[kx],  dYdxRow.data(), nCells*sizeof(double));
//    }
//
//
//    // compute total divU and update line
//    ix = 0;
//    double divU2, divU3;
//    for(it = cells.begin(); it !=cells.end();it++){
//
//        double pre1 = 1./rho[ix]/cp[ix]/T[ix];
//        double sum = 0.0;
//        double sum2 = 0.0;
//        // column major traversal, cannot be avoided
//        for(size_t kx = 0; kx < nSpecies; kx++){
//            sum += js[kx][ix] *dhsdx[kx][ix];
//            sum2+= M[ix]/lem.Ms[kx]*djsdx[kx][ix];
//        }
//        divU2 = -pre1*(dQdx[ix] + sum);  // heat + species enthalpy diffusion
//        divU3 = -sum2/rho[ix];           // Fickian diffusion
//
//  ///   if(divU2 < -1) {
//  ///       std::cout << divU2 << "\n";
//  ///       std::cout << "dQdx " << dQdx[ix] << "\n"
//  ///                 << " pre1 " << pre1 << "\n"
//  ///                 << "sum " << sum << "\n"
//  ///                 << "divU1" << divU1[ix]<< "\n";
//
//
//
//  ///   }
//
//        (*it).divU = divU1[ix] + divU2 + divU3;
//  //    (*it).divU =  divU2 + divU3;
//        ix++;
//    }
//
//}


/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       writes state of this LEM line to file
 * @author      A.menon
 *---------------------------------------------------------------------------*/

void LEMLINE::writeProps()
{

    std::list<CellData>::iterator it;

    tfile_lem << "#--- \t  \t LEM Line data  --- " << "\n"
        << "#(1) wafer"<< std::setw(19)
        << "(2) x"   << std::setw(19)
        << "(3) dx"   << std::setw(19)
        << "(4) m"   << std::setw(19)
        << "(5) p"    << std::setw(19)
        << "(6) rho"  << std::setw(19)
        << "(7) m"    << std::setw(19)
        << "(8) T"    << std::setw(19)
        << "(9) Z"    << std::setw(19)
        << "(10) C"    << std::setw(19)
        << "(11) dCdt" << std::setw(19)
        << "(12) heatRel";

    for(int ix = 0; ix < nSpecies; ix++)
    {
        tfile_lem << std::setw(19)
                  << "(" << ix + 13 << ") Y_"
                  << lem.gas->speciesNames()[ix];
    }
    tfile_lem << "\n";

    it = cells.begin();
    double x = -(*it).dx/2;

    for(it = cells.begin(); it != cells.end(); ++it){
        x+= (*it).dx;
        tfile_lem << (double) std::distance(cells.begin(), it)
                  << std::setw(19)
                  << x << std::setw(19)
                  << (*it).dx << std::setw(19)
                  << (*it).m << std::setw(19)
                  << (*it).p  << std::setw(19)
                  << (*it).rho<< std::setw(19)
                  << (*it).m  << std::setw(19)
                  << (*it).T  << std::setw(19)
                  << (*it).Z  << std::setw(19)
                  << (*it).C  << std::setw(19)
                  << (*it).dCdt << std::setw(19)
                  << (*it).dedt;

        for(int ix = 0; ix < nSpecies; ix++)
        {
            tfile_lem << std::setw(19) << (*it).Y[ix] ;
        }

        tfile_lem << "\n";
    }
}

/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       writes scalars conditioned on Z to file
 * @author      A.menon 7th Sep, 2020
 *---------------------------------------------------------------------------*/
void LEMLINE::writeZProps()
{

    tfile_lem << "#--- \t  \t LEM Line data  --- " << "\n"
        << "#(1) Z"     << std::setw(19)
        <<  "(2) rho"   << std::setw(19)
        <<  "(3) T"     << std::setw(19)
        <<  "(4) dedt" << std::setw(19)
        <<  "(5) DivU";

    for(size_t ix = 0; ix < nSpecies; ix++)
    {
        tfile_lem << std::setw(19)
                  << "(" << ix + 6 << ") Y_"
                  << lem.gas->speciesNames()[ix];
    }
    tfile_lem << "\n";

    for(size_t ix = 0;ix < ZBinSize; ix++){
        tfile_lem <<    Zstep *(ix) << std::setw(19)
                  << scalarZ[ix][IRho] << std::setw(19)
                  << scalarZ[ix][IT] << std::setw(19)
                  << scalarZ[ix][ITheta] << std::setw(19)
                  << scalarZ[ix][Ih];

        for(size_t jx = 0; jx < nSpecies; jx++)
        {
            tfile_lem << std::setw(19) << scalarZ[ix][IY+jx] ;
        }
        tfile_lem << "\n";
    }
}

/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       writes scalars conditioned on C to file for given bin
 * @author      A.menon 29th Nov 2021
 *---------------------------------------------------------------------------*/
void LEMLINE::writeCProps()
{

    tfile_lem << "#--- \t  \t LEM Line data  --- " << "\n"
        << "#(1) C"     << std::setw(19)
        <<  "(2) rho"   << std::setw(19)
        <<  "(3) T"     << std::setw(19)
        <<  "(4) h"     << std::setw(19)
        <<  "(5) dCdT init"  << std::setw(19)
        <<  "(6) dCdt"  ;

    for(size_t ix = 0; ix < nSpecies; ix++)
    {
        tfile_lem << std::setw(19)
                  << "(" << ix + 7 << ") Y_"
                  << lem.gas->speciesNames()[ix];
    }

    for(size_t ix = 0; ix < nSpecies; ix++)
    {
        tfile_lem << std::setw(19)
                  << "(" << ix + 7 + lem.ns << ") prod_"
                  << lem.gas->speciesNames()[ix];
    }

    tfile_lem << "\n";

    for(size_t ix = 0;ix < CBinSize; ix++){
        tfile_lem <<    Cstep *(ix) << std::setw(19)
                  << scalarC[ix][IRho] << std::setw(19)
                  << scalarC[ix][IT] << std::setw(19)
                  << scalarC[ix][Ih] << std::setw(19)
                  << scalarC[ix][IC2] << std::setw(19)
                  << scalarC[ix][IC];

        for(size_t jx = 0; jx < nSpecies; jx++)
        {
            tfile_lem << std::setw(19) << scalarC[ix][IY+jx] ;
        }

        for(size_t jx = 0; jx < nSpecies; jx++)
        {
            tfile_lem << std::setw(19) << speciesProdC[ix][jx] ;
        }
        tfile_lem << "\n";
    }

}

/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       writes scalars conditioned on C from scalarZC matrix to file
 * for given bin
 * @author      A.menon 23th Sep, 2021
 *---------------------------------------------------------------------------*/
void LEMLINE::writeZCProps()
{

    tfile_lem << "#--- \t  \t LEM Line data  --- "<< "\n"
        << "#(1) Z"<< std::setw(19)
        <<  "(2) C"<< std::setw(19)
        <<  "(3) T"<< std::setw(19)
        <<  "(4) h"<< std::setw(19)
        <<  "(5) dCdT"<< std::setw(19)
        <<  "(6) dedt"<< std::setw(19)
        <<  "(7) rho"  << std::setw(19)
        <<  "(8) FOUND"<<std::setw(19)
        <<  "(9) Found_dcdt";

    for(size_t ix = 0; ix < nSpecies; ix++)
    {
        tfile_lem << std::setw(19)
                  << "(" << ix + 10 << ") Y_"
                  << lem.gas->speciesNames()[ix];
    }
    tfile_lem << "\n";

    for(int ix = 0;ix < ZBinSize; ix++){
        double zNow = ix * (1./ZBinSize) + 0.5/ZBinSize;
        zNow = ix * Zstep;

        for(int jx = 0; jx < CBinSize; jx++){
            double cNow = jx * (1./CBinSize) + 0.5/CBinSize;
            cNow = jx * (Cstep);

            tfile_lem << zNow << std::setw(19)
                      << cNow << std::setw(19)
                      << scalarZC[ix][jx][IT] << std::setw(19)
                      << scalarZC[ix][jx][Ih] << std::setw(19)
                      << scalarZC[ix][jx][IC] << std::setw(19)
                      << scalarZC[ix][jx][ITheta] << std::setw(19)
                      << scalarZC[ix][jx][IRho] << std::setw(19)
                      << scalarZC_flag[ix][jx]<<std::setw(19)
                      << dCdt_flag[ix][jx];

            for(size_t kx = 0; kx < nSpecies; kx++)
            {
                tfile_lem << std::setw(19) << scalarZC[ix][jx][IY+kx] ;
            }
            tfile_lem << "\n";
        }

        tfile_lem << "\n";
    }

    tfile_lem << "#-- end\n";

}

/**----------------------------------------------------------------------------
 * @author      A.menon 2/1/2024
 *---------------------------------------------------------------------------*/
void LEMLINE::writeBetaPDFS(double m, double v)
{

    tfile_lem << "#--- \t  \t LEM Line data  --- "<< "\n"
        << "#(1) x"<< std::setw(19)
        << "(2) I(x)"<< std::setw(19)
        <<  "(3) p(x)"<< "\n";
    std::vector<double> fx = interpBeta(m,v);
    size_t bins = lem.BetaShapes[0][0].size();

    for(size_t kx = 0; kx < bins-1; kx++){
        double x = kx/(bins-1.);
        tfile_lem << x <<std::setw(19);
        tfile_lem << fx[kx]<<std::setw(19);
        tfile_lem << (fx[kx+1]-fx[kx])*(bins-1);
        tfile_lem << "\n";
    }


    tfile_lem << "#-- end\n";

}

/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       writes scalarZC and scalarZC_flag to disk
 * @author      A.menon 24 Oct, 2023
 *---------------------------------------------------------------------------*/
void LEMLINE::writeZCDat(){
    int RSize = scalarZC[0][0].size();
    tfile_lem << ZBinSize << ' ' << CBinSize << ' ' << RSize << '\n';
    for(int ix = 0;ix < ZBinSize; ix++){
        for(int jx = 0; jx < CBinSize; jx++){
            for(size_t kx = 0; kx < RSize; kx++){
                tfile_lem << scalarZC[ix][jx][kx]<<' ';
            }
        }
    }
}

void LEMLINE::writeZCFlagDat(){
    for(int ix = 0;ix < ZBinSize; ix++){
        for(int jx = 0; jx < CBinSize; jx++){
            tfile_lem << scalarZC_flag[ix][jx]<<' ';
        }
    }
}

/**----------------------------------------------------------------------------
 * @param       /none
 * @brief       writes scalarZC and scalarZC_flag to disk
 * @author      A.menon 24 Oct, 2023
 *---------------------------------------------------------------------------*/
void LEMLINE::writeCDat(){
    int RSize = scalarC[0].size();
    tfile_lem <<  CBinSize << ' ' << RSize << '\n';
    for(int jx = 0; jx < CBinSize; jx++){
        for(size_t kx = 0; kx < RSize; kx++){
            tfile_lem << scalarC[jx][kx]<<' ';
        }
    }
}



/**-----------------------------------------------------------------------------
 * @brief       Sets the dtStirr based on A.R Kerstein 1988
 * @author      M.Oevermann
 *
 * EDIT:        Eddy timing is now offloaded to sampleEddyTime()
 * EDIT:        returs false if insufficient resolution, true otherwise
 *---------------------------------------------------------------------------*/
bool LEMLINE::setdtStirr(double velocityScale)
{

    ut= velocityScale;
    double nu= getMeanKinematicViscosity(); // Kinematic viscosity (m^2)/sec.

    double Re_t = ut*lt/nu;
    eta         = N_ETA*lt*pow(Re_t, -3./4.);

    double lambda  = 54./5. * nu*Re_t/(lt*lt*lt * C_LAMBDA)
                            *(pow(lt/eta, 5./3.)-1.)/(1.-pow(eta/lt,4./3.));

    double eddyFreq  = lambda*lengthOfLine();

    dtStirr = 1./eddyFreq;

    // -- Check line for resolution
    std::list<CellData>::iterator it;
    for (it = cells.begin(); it != cells.end(); ++it) {
        if (eta < (*it).dx * 6.2){
            return false;
        }
    }

    return true;
}
/**-----------------------------------------------------------------------------
 *
 *---------------------------------------------------------------------------*/
void LEMLINE::setZMax(double ZM){

    ZMax = ZM;
    Zstep = ZM/(ZBinSize - 1.);

}


/**-----------------------------------------------------------------------------
 * @brief       Sets the dtStirr based on A.R Kerstein 1988
 * @author      M.Oevermann
 *
 * EDIT:        Eddy timing is now offloaded to sampleEddyTime()
 * EDIT:        returs false if insufficient resolution, true otherwise
 *---------------------------------------------------------------------------*/
bool LEMLINE::setdtStirrRe(double Re_t)
{
    if(Re_t < 1e-12){
        dtStirr = 1e9;
        eta = lt;
        return true;
    }


    double nu= getMeanKinematicViscosity(); // Kinematic viscosity (m^2)/sec.
    eta = N_ETA*lt*pow(Re_t, -3./4.);

    double lambda  = 54./5. * nu*Re_t/(lt*lt*lt * C_LAMBDA) *(pow(lt/eta, 5./3.)-1.)/(1.-pow(eta/lt,4./3.));

    double eddyFreq  = lambda*lengthOfLine();
    dtStirr = 1./eddyFreq;
    return true;
}

/**-----------------------------------------------------------------------------
 * @brief       Samples eddy time from Poisson process with a mean of dtStirr
 * @author      A. Menon
 *--------------------------------------------------------------------------*/
double LEMLINE::sampleEddyTime(){

    // -- generate random number [0,1]
    double F = (double) (rand()+1e-12) / RAND_MAX;
    // -- inverse CDF to find inter-eddy space
    double eddyGap = -log(F)*dtStirr;
    return eddyGap;
}

/**
 * @brief       advances chemistry using CVODE at 1/10th of chemStep,
 * conditions mass fractions based on C (O2 mass fractions)
 * OBS - species prod is not done here for now, might need a re-write to
 * include that
 *
 * @param       dt      \input          time step
 *
 * @author      A.Menon 25 Apr 2022
 *
 * EDIT:: species prod as well as C source term conditioning were added (check
 * github)
 */
void LEMLINE::advanceSourceTerm2(const double &dt)
{
    realtype tret;
    std::list<CellData>::iterator it;

    // -- advance chemistry and condition at the same time
    int indO2 = lem.gas->speciesIndex("O2");
    double Yu = Yu_p[indO2];
    double Yb = Yb_p[indO2];
    double hitCount[CBinSize]{};

    auto gas=lem.gas;
    double* Ms   = lem.Ms;
    double* RR   = lem.tmp_s1;

    // -- initialize temp to Zero
    size_t Ssize = tempC[0].size();
    for(size_t ix = 0; ix < CBinSize; ix++)
      for(size_t jx = 0; jx < Ssize ;jx++)
          tempC[ix][jx] = 0.0;

    // -- scan the line
    for (it = cells.begin(); it != cells.end(); ++it)
    {
        realtype *y = NV_DATA_S(lem.y);
        y[IT]       = (*it).T;
        y[IP]       = (*it).p;
        y[IRHO]     = (*it).rho;

        for (int s = 0; s < lem.ns; ++s) {
          y[IY0+s] = ((*it).Y)[s];
        }

        // -- init solver for this LEM cell
    //  CVodeReInit(lem.cvode_mem, time, lem.y);
        for(size_t tx = 1; tx < 21;tx++){
     //     CVode(lem.cvode_mem, time+dt/20*tx, lem.y, &tret, CV_NORMAL);
            // -- compute immediate C based on O2 mass fractions
            double C = (y[IY0+indO2]-Yu)/(Yb - Yu);
         // C = C < 0.0 ? 0.0:C;
         // C = C > 1.0 ? 1.0:C;
            //-- find production rates and dCdT
            gas->setState_TPY(y[IT], y[IP], &y[IY0]);
            lem.kin->getNetProductionRates(RR);// units Kmol/m^3/s
            double  rho = gas->density();
            double dCdt = RR[indO2]*Ms[indO2]/rho/(Yb - Yu);

            // -- condition
            int pos  = C/Cstep;
            pos = pos <      0      ?   0      :pos;
            pos = pos > (CBinSize-1)?CBinSize-1:pos;

            hitCount[pos]++;

            // -- temperature, Ys and species prod
            tempC[pos][IC]     += dCdt;
            tempC[pos][IT]     += y[IT];
            for(size_t jx = 0; jx < nSpecies; jx++){
                tempC[pos][IY+jx] += y[IY0+jx];
            }

            for(size_t jx = 0; jx < nSpecies; jx++){
                tempC[pos][IY+jx+nSpecies] += RR[jx];
            }

        }
        (*it).T   = y[IT];
        (*it).p   = y[IP];
        (*it).rho = y[IRHO];

        for (int s = 0; s < lem.ns; ++s){

            ((*it).Y)[s] = y[IY0+s];
        }
        (*it).dx = (*it).m / (*it).rho;
    }

    // -- exit scan and condition
    for(size_t ix = 0; ix < CBinSize; ix++){
      if(hitCount[ix] > 0){
          scalarC_flag[ix] = 1.0; // update the global flag matrix
          scalarC[ix][IC]  = bFac*(tempC[ix][IC]/hitCount[ix])
                             +(1.0-bFac)* scalarC[ix][IC];
          scalarC[ix][IT]  = bFac*(tempC[ix][IT]/hitCount[ix])
                             +(1.0-bFac)* scalarC[ix][IT];
          for(size_t jx = 0; jx < nSpecies; jx++){
              scalarC[ix][IY+jx] = bFac*(tempC[ix][IY+jx]/hitCount[ix])
                                 +(1.0-bFac)* scalarC[ix][IY+jx];
          }

          for(size_t jx = 0; jx < nSpecies; jx++){
              speciesProdC[ix][jx] = bFac*(tempC[ix][IY+jx+nSpecies]/hitCount[ix])
                                 +(1.0-bFac)* speciesProdC[ix][jx];
          }
      }
    }
}

