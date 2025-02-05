#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include <mpi.h>
#include "LEMLine.h"
#include "CellData.h"
#include "yaml-cpp/yaml.h"


#define NETA    4
#define CLAMBDA 15
#define Z_st 0.055


int main(int argc, char** argv){
    MPI_Init(NULL,NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    int      NLEM;
    int      zBinSize,cBinSize;
    int      spliceFreq;
    double   ZMax,tStep, tLES, domLength, res;
    bool     combFlag, solveFast ;
    std::ofstream file;

    // {{{-- YAML file parser

    YAML::Node config   = YAML::LoadFile("input.yaml");
    std::string fname   = config["file"].as<std::string>();


        // -- gas properties
    LEM oj(fname); // LEM object contains a cantera gas object
    oj.fuelStream  = config["fuelStream"].as<std::string>();
    oj.oxStream    = config["oxStream"].as<std::string>();
    oj.fuelStreamT = config["fuelStreamT"].as<double>();
    oj.oxStreamT   = config["oxStreamT"].as<double>();
    oj.gammasCoeff = config["gammas"].as<std::vector<double>>();
    oj.pressure    = config["pressure"].as<double>();


        //-- Simulaton parameters
    NLEM           = config["NLEM"]     .as<int>();
    ZMax       = config["ZMax"] .as<double>();
    tStep          = config["tStep"]    .as<double>();
    tLES           = config["tLES"]     .as<double>();
    domLength      = config["domLength"].as<double>();
    zBinSize       = config["zBinSize"] .as<int>();
    cBinSize       = config["cBinSize"] .as<int>();
    res          = config["res"].as<double>();
    spliceFreq     = config["spliceFreq"].as<int>();

    combFlag       = config["combFlag"].as<bool>();
    solveFast      = config["solveFast"].as<bool>();
    // }}}

    srand(time(NULL));          // seed for RNG

    if(world_rank==0) std::cout  << config << std::endl;
    NLEM = domLength/res;
    double lt = 1e-3;
    
    LEMLINE lineRef(NLEM, domLength ,lt ,file, oj, zBinSize, cBinSize);
    lineRef.initZCMemory();
    lineRef.initCMemory();
    lineRef.setZMax(ZMax);


    if(world_rank==0)std::cout << "\nCalculating..." << std::endl;
    
    for(int ix= world_rank ; ix < zBinSize; ix+=world_size){
        double Z = ix*lineRef.Zstep;
     // double Z = 0.1;
        
        LEMLINE line1(NLEM, domLength ,lt ,file, oj, zBinSize, cBinSize);
        line1.N_ETA = NETA;
        line1.C_LAMBDA = CLAMBDA;
        std::vector<std::string> C_list {"CO2", "H2O"};
        line1.setC_select(C_list);
        line1.setStreams();
        line1.initZCMemory();
        line1.initCMemory();
        line1.setZMax(ZMax);
        line1.setEquilTY_Z();
        
        double phi = Z*(1.-Z_st)/Z_st/(1.-Z);
        if(phi!=0.0){
            line1.setPremixedData(phi);
            line1.initScalarC_linear();
        }
        else {
            line1.setPremixedData(1e-5);
            line1.initScalarC_linear();
            lineRef.copyMats(line1);
            lineRef.calcScalarReC(Z);
            continue;
        }


        line1.initializeWithZ(Z);
        line1.sparkLine(0.25);
        
        double tb = oj.TEquil_z[ix];
        double tu = oj.TUnbZ[ix];

        double lPrev = line1.lengthOfLine();
        line1.setdtStirrRe(0.0);
        double tx = 0;
        int fx = 0;
        for(double fill = 0.0; (tx < tLES) ; tx += tStep){
            std::list<CellData> dummyList;
            LEMLINE inLine(NLEM, domLength,lt ,file, oj, 0, 1);
            inLine.initZCMemory();
            inLine.initializeWithZ(Z);

            if(solveFast)
               line1.advanceLEM3(tStep,  combFlag);
            else
               line1.advanceLEMPremixed(tStep,  combFlag);

            double lNow = line1.lengthOfLine();
            double spliceLength = lNow-lPrev;
            double lenghtFrac = spliceLength;
            
            line1.spliceToListLB(lenghtFrac,dummyList);
            
            std::list<CellData> dummyList2;

            spliceLength *= tu/tb;
            
            lenghtFrac = spliceLength;
            inLine.spliceToListLB(lenghtFrac,dummyList2);
            line1.spliceFromList(dummyList2);
            line1.makeEquidistantGridafterSplicing(lNow/res, 1e-9);
          
          //std::string outFileName = "./results/" + std::to_string(fx) + ".dat";
          //file.open(outFileName);
          //file << "# Time "    << tx           << std::endl;
          //file << "# eta = "   << line1.eta    << std::endl;
          //file << "# length = "<< line1.length << std::endl;
          //line1.writeProps();
          //file.close();
          //fx++;
        }
        lineRef.copyMats(line1);
        lineRef.calcScalarReC(Z);
        std::cout <<"Done Zbin = " << ix << "Z = " << Z<<std::endl;
    }
        
//  line1.fillZC_biased();
//  line1.initScalarZC_linear();

    int RSIZE = lineRef.scalarZC[0][0].size();
    std::vector<double> BUFF(zBinSize*cBinSize*RSIZE,0.0);
    std::vector<double> RECBUFF(zBinSize*cBinSize*RSIZE,0.0);
    if(world_size>0){
        for(int ix = 0; ix <zBinSize; ix++){
            for(int jx = 0; jx <cBinSize; jx++){
                for(size_t kx = 0; kx < RSIZE; kx++){
                    BUFF[ix*(cBinSize*RSIZE) + jx*RSIZE + kx] = lineRef.scalarZC[ix][jx][kx];
                }
            }
        }
        if(world_rank==0)std::cout << "All reduce (SUM) for all procs\n!";
        MPI_Allreduce(BUFF.data(), RECBUFF.data(), zBinSize*cBinSize*RSIZE, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

        for(int ix = 0; ix <zBinSize; ix++){
            for(int jx = 0; jx <cBinSize; jx++){
                for(int kx = 0; kx < RSIZE; kx++){
                  lineRef.scalarZC[ix][jx][kx]= RECBUFF[ix*(cBinSize*RSIZE) + jx*RSIZE + kx];
                }
              lineRef.scalarZC_flag[ix][jx] = true;
            }

        }
    }
        
    if(world_rank == 0){
        file.open("initZC.dat");
        lineRef.writeZCProps();
        file.close();
        file.open("scalarZC.dat");
        lineRef.writeZCDat();
        file.close();
        file.open("scalarZC_flag.dat");
        lineRef.writeZCFlagDat();
        file.close();
    }

    MPI_Finalize();
    return 0;
}
