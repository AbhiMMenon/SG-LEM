ofstream    file;       // for writing out state of any given line

// separate output file for each processor

if(Pstream::myProcNo()==0) {
    std::string fileName = "./LineProps/";
    std::string makeDir = "mkdir " + fileName;
    Info << "FILENAME = " << fileName << Foam::endl;
    system(makeDir.c_str());
}

int fileNo = 0;


// LEM gas object
LEM oj(fname);
int nSpecies = oj.ns;

// -- gas properties
  oj.fuelStream       = fuelStream;
  oj.oxStream         = oxStream;
  oj.fuelStreamT      = fuelStreamT;
  oj.oxStreamT        = oxStreamT;
  oj.pressure         = pressure;

srand(time(NULL));          // seed for RNG


// -- indices for prodSpecies
std::vector<int> pList(prodSpecies.size());
Info << " The positions are: ";
forAll(Yprod, i){
    pList[i] = oj.gas.speciesIndex(prodSpecies[i]);
    Info << pList[i] << " ";
}

const compressible::LESModel& lesModel =
    Y[0].db().lookupObject<compressible::LESModel>
    (
        momentumTransportModel::typeName
    );

const volScalarField& delta = lesModel.delta();
