#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "LEMLine.h"

/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
int main ()
{
    std::list<CellData>::iterator it;
    std::string ctiFile = "gri30.yaml";
    LEM oj(ctiFile, "");
    std::ofstream file;

    int nLEM    = 40;
    int nRealz  = 5;
    double tEnd = 1e-5;
    int zBinSize = 10, cBinSize = 10;

    LEMLINE lemLine1(nLEM, 0.04 ,0.04 ,file, oj, zBinSize, cBinSize);
    lemLine1.initCMemory();

    oj.fuelStream=  "CH4:1";
    oj.oxStream=  "O2:1,N2:3.76";
    oj.fuelStreamT=  oj.oxStreamT=  300;
    std::vector<std::string> C_list {"CO2", "H2O"};

    lemLine1.setC_select(C_list);
    lemLine1.initZCMemory();
    lemLine1.initCMemory();
    lemLine1.initializeWithPhi(0.5);
    file.close();
    return 0;
    lemLine1.initScalarC();
    lemLine1.initializeDataOnLineForPremixedFlameTest(160);
    
    file.open("lemdata.dat");
    file.close();
    return 0;

    for (int i=0; i<nRealz; i++)
    {
        file << "# Time "    << nRealz           << std::endl;
        file << "# eta = "   << lemLine1.eta    << std::endl;
        file << "# length = "<< lemLine1.length << std::endl;
        lemLine1.advanceLEMPremixed(tEnd,  false);
        std::cout << "Done realization # " << i << std::endl;
        lemLine1.writeProps();

    }

    file.close();

    return 0;
}
