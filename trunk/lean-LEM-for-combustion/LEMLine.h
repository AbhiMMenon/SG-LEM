#ifndef LEMLINE_H
#define LEMLINE_H

//-- positions for conditioning matrix
#define IT      0
#define Ih      1
#define ITheta  2
#define IRho    3
#define IC      4
#define IC2     5
#define IY      6


//-- scalars other than mass fractions
#define NSCALAR  6

#define NETA_DEF    4
#define CLAMBDA_DEF 15

#include "CellData.h"
#include "LEM.h"

#include <iomanip>
#include <ctime>
#include <limits>
#include <boost/math/distributions/beta.hpp>
#include "omp.h"


//using namespace std;
//using namespace Cantera;


/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
class LEMLINE
{
  public:

    // Note: The constructor only creates an LEM line with
    // _nlem elements WITHOUT initializing, i.e. without any
    // meaningfull physical data. Initialization of the lem line
    // with physical data must be done after construction, e.g.
    // via member function initializeDataOnLine()
    LEMLINE(
        const int       _nlem,
        const double    _length,
        const double    _lt,
        std::ofstream&       tfile,
        LEM             &_lem,

        // -- A.Menon ↓
        int     _ZBinSize       = 200,
        int     _CBinSize       = 50
        );

//  LEMLINE(const LEMLINE &lemline);
   ~LEMLINE();

    void initializeDataOnLine();
    void initializeDataOnLineForMixingLayerTest();
    void initializeDataOnLineForPremixedFlameTest(long unsigned int nlem_div);


    void advanceLEM(const double &dt, double threshold=1e-5); // threshold added 19th jul 2021, A.Menon
    void advanceSourceTerm(const double &dt);
    void advanceSourceTerm2(const double &dt); // fractional step conditioning, A.Menon
    void advanceDiffusion(const double &dt);
    void advanceDiffDiffusion(const double &dt); // differential diffusion A.Menon
    void eddyEvent();
    void correctExpansion();
    double internalEnergyOnLine();
    double massOnLine();
    double lengthOfLine();
    int eddyLocation();
    int eddySize();
    void tripletMap( unsigned int iStart, const unsigned int nEddy);
    void solveWithThomasAlgorithm();
    void splice(double m, LEMLINE &toLine);
    void spliceToList(double m, LEMLINE &toLine, std::list<CellData>& SplicedMassList);
    void EraseFromOutlet(double m);
    void AddOneLEMCell(double m_in, double T_in, double p_in, std::vector<double> Y_in);

    bool setdtStirr(double velocityScale);
    bool setdtStirrRe(double Re_t);
    double sampleEddyTime();

    void initializeDataOnLineFromLES(double T_in, double P_in,std::vector<double> Y_in);
    void setLEMT(double T_in);

    std::vector<double> getSpeciesMeanMixtureFractions();
    double getMeanKinematicViscosity();
    double getMeanDynamicViscosity();
    double getMeanTemperature();
    double getMeanTemperatureCpWeighted();
    double getMinTemperature();
    double getMaxTemperature();
    double getMaxYfuel(int fuelindex);
    double getMeanDensity();
    double getMeanPressure();
    double getMeanGasConstant();
    double getMeanEnthalpy();
    double getMeanFormationEnthalpy();
    std::vector<double> getSpeciesMeanNetProductionRate();
    std::vector<double> getSpeciesAndHeatMeanNetProductionRate();
    std::vector<double> getSpeciesAndHeatMedianNetProductionRate();
    std::vector<double> getLEMcellsTemperature();
    double getMeanHeatProductionRate();
    double getMedianHeatProductionRate();
    double calculateMedian(std::vector<double> V);
    int getTotalNoofEddies();
    double getAverageEddiesLength();
    double getLEMDomainLength();
    int getLEMDomainCells();
    void setdtStirrRANS(double DiffTurbLEM, double L_kalmLEM , double
            timeStep_LES);


    friend std::ostream& operator<< (std::ostream& os, const LEMLINE &lem);
    friend std::istream& operator>> (std::istream& is, LEMLINE &lem);
    void binary_write(std::ofstream& o);
    void LEMoutput_write(std::ofstream& o, double timeNow,int noLES);
    void binary_read(std::ifstream& i);
    void binary_read_wo_ini(std::ifstream& i);


/*-----------------------------------------------------------------------------
 * All the RILEM stuff for SG-LEM, conditioning, integration, 0D reactor etc..
 * --------------------------------------------------------------------------*/
    void writeProps();
    void writeZProps();
    void writeCProps();
    void writeZCProps();
    void writeBetaPDFS(double m, double v);
    void writeZCDat();
    void writeCDat();
    void writeZCFlagDat();
    void setStreams();
    void setPremixedData(double phi);
    void setEquilTY_Z();
    void copyMats(LEMLINE &ref);
    void copyZCMats(LEMLINE &ref);
    void clearCmats();
    void setBilgerMixFrac();
    void setTracerMixFrac(std::string inertSpecie);
    void setC_select(std::vector<std::string> C_list);

    // OBS flaglinear does lead to arifacts in the C calculation
    void setProgressVariable_T(bool flagLinear = false, double dt = 1.0); // NOI
    void setProgressVariable_simple(double dt); // NOI
    void setProgVble_premix(double dt);
    void setProgVble_premixO2();
    void setProgVble_premixSp();
    void setProgVble_nonpremixSp();
    std::vector<double> getProgVble_nonpremixSp(double Z, double *Y, double T);
    void setHeatRelease();
    void setEnergy();
    void setdCdt(); // empty
    void setDivU(); //NOI
    void setZMax(double ZM); // can be used to set max Re_t in flamelet code

    // memory initialization
    void initZMemory();
    void initZCMemory();
    void initCMemory();


    // matrix initial value methods
    void initScalarZ();
    void initScalarC();
    void initScalarC_linear();
    void initScalarZC_zero();
    void initScalarZC_linear();
    void initScalarZC_unburnt();
    void initScalarZC_read();
    void initScalarC_read();
    void initScalarZC_fromCFDSum(double Z, double C, double T, double *Y, double dcDt);
    void initScalarZC_fromCFDNorm();
//  void initScalarZC_0D(double flamLimit = 1.0);
//  void initScalarC_0D();



    // -- Tabulated BetaPDFs
    void initBetaShapes(size_t bins,double mMax=1.0);
    std::vector<double> interpBeta(double mean, double var,double mMax=1.0);


    // conditioning methods
    void calcScalarZ();
    void calcScalarC();
    void calcSpeciesC();
    int  condProgVble();
    void calcScalarZC();
    void calcScalarReC(double Re_t);
    void interpCRe(double Re_t);

    double returnHsource(double *dYdt);  //NOI
    double getBilgerMixFrac(double *Y);
    double getMeanMixFrac(); //
    double explicitCsource(double *dYdt, double Z); //YES
    double getCFullness();
    double getProgVble_premixO2(double Y_O2);
    std::vector<double> getMinMaxMixFrac(); //NOI
    std::vector<double> getElementMassFrac(double *Y);


    // PDF weighted integration
    std::vector<double> getPDFMeanSZ(double zMean, double zVar, double press);
    std::vector<double> getPDFMeanSC(double cMean, double press,std::vector<int> pList);
    std::vector<double> getPDFMeanSC_Beta(double cMean, double cVar,std::vector<int> pList);
    std::vector<double> getPDFMeanSC_TH(double cMean, double cVar,std::vector<int> pList);
    std::vector<double> getPDFMeanS2(double zMean, double zVar, double cMean,double press);
    std::vector<double> getPDFMeanS2_TH(double zMean, double zVar, double cMean);
    std::vector<double> getPDFMeanS4(double zMean, double zVar, double cMean,double press); // PNOI
    std::vector<double> getPDFMeanS5(double zMean, double zVar, double cMean,double press); // PNOI
    std::vector<double> getPDFMeanSBV(double zMean, double zVar, double cMean,double cVar);

    // Interpolation
    std::vector<double> getInterpS(double zMean,double cMean,double press,bool uFlag); // YES
    std::vector<double> getInterpS_fillC(double zMean,double cMean); // YES
    double getdCdt_fillC(double zMean,double cMean); // YES
    std::vector<double> getInterpSZ(double zMean,double press);
    std::vector<double> getInterpSC(double cMean,double CFDpress, std::vector<int> pList);
    std::vector<double> getInterpSC_incomplete(double cMean,std::vector<int> pList); //YES

    void  fillZC_biased();

    // Averaging
    std::vector<double> getMeanS(std::vector<int> pList);
    double getCvar();
    double getCMean();


    // domain initialization
    void initializeJetTest();
    void initializeMixTest(); //NOI
    void boundaryMixTest(); //NOI
    void initializeOx();
    void initializeWithZ(double Z);
    void initializeWithPhi(double phi);
    void initializeWithZ_staggered(double Z); // redundant
    void sparkLine(double startFrac=0.5);
    void sparkPremixed();
    void sparkNonPremixed();
    void injectRadicals();
    void initCellByCell(double T_in, double P_in, std::vector<double> Y_in, int incr);

    // splicing
    void spliceFromList(std::list<CellData>& InSplicedMassList);
    void spliceFromList(std::list<CellData>& InSplicedMassList, double expRatio);
    void spliceToList(double m, std::list<CellData>& SplicedMassList); // modified from  previous version
    void regrid(double LEMres, double tolFac=10);

    // -- --  length fraction based
    void spliceToListLB(double Lfrac, std::list<CellData>& SplicedMassList);
    void spliceToListLBR(double Lfrac, double Cface, std::list<CellData>& SplicedMassList);


    // -- advancement options Or operator splitting methods
    void advanceLEM0(const double &dt); // no conditioning,no combustion, diffuse-eddy-diffuse,
    void advanceLEM1(const double &dt, bool chemFlag=true); // no conditioning, diffuse-eddy-diffuse
    void advanceLEM2(const double &dt, double threshold=1e-5, bool fastOn=true, bool chemFlag=true); // simplified algorithm
    void advanceLEM3(const double &dt, bool chemFlag=true); // only fast version
    void advanceLEM3DD(const double &dt, bool chemFlag=true);// fast version with mix-averaged diffusive coeffs
    void advanceLEM4(const double &dt, bool chemFlag=true); // only diffuse-eddy-diffuse
    void advanceLEM5(const double &dt, bool chemFlag=true); // fast(ish), for heavy mechanisms
    void advanceLEMPremixed(const double &dt, bool chemFlag=true);
    void advanceLEMFastPremixed(const double &dt, bool chemFlag=true);
    void advanceLEMFast(const double &dt, double &threshold);
    void advanceLEMFast2(const double &dt, bool chemFlag);



    // miscellaneous functions
    // inline double binToC(size_t ix,double step) {return std::pow(step*ix-1,3) +1;} // half sigmoid
    inline double binToCSigmoid(size_t ix) {return 1./(1. + std::exp(-(Cstep*10*ix-5)));} // full sigmoid
    inline double binToCLinear(size_t ix) {return Cstep*ix;} // linear
  //void linearPTcorrect(double pNew);

/*-----------------------------------------------------------------------------
 * Functions added by S.Arshad?
 * --------------------------------------------------------------------------*/

    void CorrectionLocalPressure(double Pnew);
    void CorrectionGlobalPressure(double Pnew);
    void CorrectionGlobalPressureAssumingAdiabaticProcess(double Pnew);
    double CorrectLEMTemperatures(double T, double Tlem, double Tmin, double
            Tmax);
    double CorrectLEMTemperaturesByFactor(double T, double Tlem, double Tmin,
            double Tmax);
    double CorrectLEMSingleSpecies(int s, double Yles, double Ylem);
    void CorrectLEMDensities();

    void initializeRestartData(double length_restart, double Energy_restart,
            double time_restart,
            std::vector<double> dx_restart,
            std::vector<double> T_restart, std::vector<double> m_restart,
            std::vector<double> E_restart, std::vector<double> rho_restart,
            std::vector<double> p_restart,
            std::vector<std::vector<double>> Y_restart);

    bool noLEMcells();
    void operator=(const LEMLINE &lemline);
    void makeEquidistantGridafterSplicing(int noLEMrequired, double minCellSize);
    void fuelMapping(double EvapMassCFDcell);


/*-----------------------------------------------------------------------------
 *                              Data for LEMLine
 * --------------------------------------------------------------------------*/


    std::list<CellData> cells;

    std::ofstream&   tfile_lem;
    LEM&        lem;

    int         nlem;                // number of grid points on line
    double      length;              // length of line
    double      Energy;              // internal energy of the line
    int         counter;

    double      time;                 // simulation time on LEM line
    double      timeNextTripletMap;   // time of next triple map event
    double      dtLEM;                // dt used on LEM line
    double      dtStirr;              // dt between eddy events / triplet maps
    double      ut;                   // characteristic velocity fluctuation
    double      lt;                   // integral length scale
    double      eta;                  // Kolmogorov length
                                      // Smith, Menon: PCI, 1996, 299-306
    int         noofEddies;
    double      Total_EddiesLength;
    double      Average_EddiesLength;

//-----------------------------------------------------------------------------
//                      Data additions by A.Menon
//-----------------------------------------------------------------------------

    double      N_ETA, C_LAMBDA;
    int         nSpecies;
    int         ZBinSize, CBinSize;
    double      Zstep, Cstep, Tstoich,bFac,ZMax;

    // burnt and un-burnt data for premixed cases
    double      TUnb_p, Tb_p;
    std::vector<double> Yb_p, Yu_p;

    // wall information
    bool        wallPresent;
    double      wallT;  // iso-thermal walls for now.

    std::vector<double> ySol; // solution vector for zeroD reactor
    double phi;               // equivalence ration for premixed


    // when scalarZC is being re-purposed for pretabulated flamelet code, Z becomes
    // trubulent Reynolds number
    matrix3D scalarZC;
    matrix2Dbool scalarZC_flag;
    matrix2Dbool dCdt_flag;
    matrix2D  scalarZ;
    matrix2D    tempC;
    matrix2D  scalarC;

    std::vector<double> scalarC_flag;
    matrix2D  speciesProdC;



};


//-----------------------------------------------------------------------------


#endif
