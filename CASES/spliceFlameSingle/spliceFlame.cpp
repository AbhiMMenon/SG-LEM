#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include "LEMLine.h"
#include "CellData.h"
#include "yaml-cpp/yaml.h"
#include "cantera/core.h"

#define NETA    0.4
#define CLAMBDA 6.7

/* 1-D flamelet code using LEM triplet-map stirring. The parameters are given
 * in input(dot)yaml. Note that the current version of Cantera uses (dot)yaml mechanism
 * files and so the previous (dot)cti files will have to be converted using the
 * `ctml2yaml' utility distributed with Cantera.
 *
 * The flame is stabilized by splicing out a length `spliceLengthOut' by
 * comparing the current (expanded due to combustion) domain size with the
 * original domain size, followed by splicing in `spliceLengthIn' domain
 * fragments based on the density ratio.
 *
 */


int main(int argc, char** argv){

    int      NLEM;
    int      zBinSize,cBinSize;
    double   velScale,tStep, tLES, domLength, phi,res;
    bool     combFlag, solveFast ;
    std::ofstream file;

    // {{{-- YAML file parser

    YAML::Node config   = YAML::LoadFile("input.yaml");
    std::string fname   = config["file"].as<std::string>();
   
//  fname   = "sk17/chem.yaml";
    auto sol = Cantera::newSolution(fname, "gas");
    auto gas = sol->thermo();
    std::cout << "N elements = " << gas->nElements() <<std::endl;
   
    // -- gas properties
    LEM oj(fname, "gas"); // LEM object contains a cantera gas object
                           //
                           //
    oj.fuelStream  = config["fuelStream"].as<std::string>();
    oj.oxStream    = config["oxStream"].as<std::string>();
    oj.fuelStreamT = config["fuelStreamT"].as<double>();
    oj.oxStreamT   = config["oxStreamT"].as<double>();
    oj.gammasCoeff = config["gammas"].as<std::vector<double>>();
    oj.pressure    = config["pressure"].as<double>();


        //-- Simulaton parameters
    velScale       = config["velScale"] .as<double>();
    tStep          = config["tStep"]    .as<double>();
    tLES           = config["tLES"]     .as<double>();
    domLength      = config["domLength"].as<double>();
    phi           = config["phi"]     .as<double>();
    zBinSize       = config["zBinSize"] .as<int>();
    cBinSize       = config["cBinSize"] .as<int>();
    res          = config["res"].as<double>();

    combFlag       = config["combFlag"].as<bool>();
    solveFast      = config["solveFast"].as<bool>();

    std::cout << "CombFlag from Yaml = " << combFlag;
    // }}}

    srand(time(NULL));          // seed for RNG

    NLEM = domLength/res;
    double lt = 2e-3;

    LEMLINE line1(NLEM, domLength ,lt ,file, oj, zBinSize, cBinSize);
    line1.makeEquidistantGridafterSplicing(NLEM,1e-9);
    line1.N_ETA = NETA;
    line1.C_LAMBDA = CLAMBDA;
    std::vector<std::string> C_list {"CO2", "H2O"};
    line1.setC_select(C_list);
    line1.initZCMemory();
    line1.initCMemory();
    line1.setZMax(1.0);
    line1.bFac = 1.0;
    line1.initializeWithPhi(phi);
    line1.initScalarC();
    line1.sparkLine(0.5);

    
    domLength = line1.lengthOfLine();
    std::cout << "\nsetdtStirr (" << velScale <<")"<< std::endl;
    line1.setdtStirrRe(velScale);
    

    double tx = 0;
    int fx = 0;

    for(; (tx < tLES) ; tx += tStep){
        std::list<CellData> dummyList;
        LEMLINE inLine(NLEM, domLength,lt ,file, oj, 0, 1);
        inLine.initCMemory();
        inLine.initializeWithPhi(phi);

        if(solveFast)
           line1.advanceLEMFastPremixed(tStep,  combFlag);
        else
           line1.advanceLEMPremixed(tStep,  combFlag);

        double lNow = line1.lengthOfLine();
        double spliceLengthOut = lNow-domLength;
        double lenghtFrac = spliceLengthOut;
        
       line1.spliceToListLB(lenghtFrac,dummyList);
        
        std::list<CellData> dummyList2;
        double tb = oj.TEquil_z[0.14/line1.Zstep];
        double tu = oj.TUnbZ[0.14/line1.Zstep];
        double spliceLengthIn = spliceLengthOut/(line1.Tb_p/line1.TUnb_p);
        lenghtFrac = spliceLengthIn;
        
        inLine.spliceToListLB(lenghtFrac,dummyList2);
        line1.spliceFromList(dummyList2);

        line1.makeEquidistantGridafterSplicing(domLength/res,1e-9);

     
        //mass based splicing, intention is still to maintain a domain
        //size
        
      //double lNow = line1.lengthOfLine();
      //double spliceLengthOut = lNow-domLength;
      //double massOut = spliceLengthOut*rhoB;
      //
      //line1.spliceToList(massOut,dummyList);
      //
      //std::list<CellData> dummyList2;
      //inLine.spliceToList(massOut,dummyList2);
      //line1.spliceFromList(dummyList2);

      
        std::string outFileName = "./results/" + std::to_string(fx) + ".dat";
        file.open(outFileName);
        file << "# Time "    << tx           << std::endl;
        file << "# eta = "   << line1.eta    << std::endl;
        file << "# length = "<< line1.length << std::endl;
        line1.writeProps();
        file.close();
        

     // outFileName = "./resultsZ/" + std::to_string(fx) + ".dat";
     // file.open(outFileName);
     // file << "# Time "    << tx           << std::endl;
     // file << "# eta = "   << line1.eta    << std::endl;
     // file << "# length = "<< line1.length << std::endl;
     // line1.writeZCProps();
     // file.close();
       
        fx++;
    }

    file.open("scalarCInit.dat");
    line1.writeCProps();
    file.close();
    
    file.open("scalarC.dat");
    line1.writeCDat();
    file.close();

    return 0;
}
