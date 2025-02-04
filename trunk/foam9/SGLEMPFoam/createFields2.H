// Additional fields for superGrid LESLEM


Info<< "Creating field SuperCell\n" << endl;
volScalarField superCell
(
    IOobject
    (
        "superCell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);

Info<< "Creating field nLEMCells\n" << endl;
volScalarField nLEMCells
(
    IOobject
    (
        "nLEMCells",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);

Info<< "Creating field nSplicings\n" << endl;
volScalarField nSplicings
(
    IOobject
    (
        "nSplicings",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);

Info<< "Creating field nRatio\n" << endl;
volScalarField nRatio
(
    IOobject
    (
        "nRatio",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);

Info<< "Creating field nEddies\n" << endl;
volScalarField nEddies
(
    IOobject
    (
        "nEddies",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);


volScalarField intFAIL
(
    IOobject
    (
        "intFAIL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);

volScalarField intHoles
(
    IOobject
    (
        "intHoles",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);


Info<< "Creating field C_LEM\n" << endl;
volScalarField C_LEM
(
    IOobject
    (
        "C_LEM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)

);

Info<< "Creating field Csource_LEM\n" << endl;
volScalarField Csource_LEM
(
    IOobject
    (
        "Csource_LEM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0)

);



Info<< "Creating field massLEM\n" << endl;
volScalarField massLEM
(
    IOobject
    (
        "massLEM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMass, 0.0)

);

Info<< "Creating field phiDiff\n" << endl;
volScalarField phiDiff
(
    IOobject
    (
        "phiDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimVolume/dimTime, 0.0)

);

Info<< "Creating field phiIn\n" << endl;
volScalarField phiIn
(
    IOobject
    (
        "phiIn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimVolume/dimTime, 0.0)

);

// --

Info<< "Creating field T_mapped\n" << endl;
volScalarField T_mapped
(
    IOobject
    (
        "T_mapped",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimTemperature, 0.0)

);

T_mapped = T;

Info<< "Creating field C_mapped\n" << endl;
volScalarField C_mapped
(
    IOobject
    (
        "C_mapped",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0)

);


Info<< "Creating field T_LEM\n" << endl;
volScalarField T_LEM
(
    IOobject
    (
        "T_LEM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimTemperature, 0.0)

);

//T_LEM.boundaryField() = T.boundaryField(); // for visualization


Info<< "Creating field H_mapped\n" << endl;
volScalarField H_mapped
(
    IOobject
    (
        "H_mapped",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimMass, 0.0)

);

H_mapped = thermo.hs();

Info<< "Creating field RHO_mapped\n" << endl;

volScalarField RHO_mapped
(
    IOobject
    (
        "RHO_mapped",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMass/dimVolume, 0.0)

);
RHO_mapped = thermo.rho();

Info<< "Creating field dQ_cantera\n" << endl;

volScalarField dQ_cantera
(
    IOobject
    (
        "dQ_cantera",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)

);

Info<< "Creating field dCdT_mapped\n" << endl;

volScalarField dCdT_mapped
(
    IOobject
    (
        "dCdT_mapped",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless/dimTime, 0.0)

);

Info<< "Creating field H_thermo\n" << endl;
volScalarField H_thermo
(
    IOobject
    (
        "H_thermo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimMass, 0.0)

);

H_thermo = thermo.hs();

Info<< "Creating field cSource\n" << endl;
volScalarField cSource
(
    IOobject
    (
        "cSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("overTime", dimMass/dimVolume/dimTime, 0.0)

);

Info<< "Reading field Cmean\n" << endl;
volScalarField Cmean
(
    IOobject
    (
        "Cmean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Cmean.dimensions().reset(dimless);

fields.add(Cmean);

volScalarField Cvar
(
    IOobject
    (
        "Cvar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// RILEM stuff
Info<< "Reading RILEMProperties\n" << endl;
IOdictionary RILEMProperties
(
     IOobject
         (
                  "RILEMProperties",
                   runTime.constant(),
                   mesh,
                   IOobject::MUST_READ,
                   IOobject::NO_WRITE
         )
);


// sampling stuff
IOdictionary sampleFaceDict
(
     IOobject
         (
                  "sampleFaceDict",
                   runTime.constant(),
                   mesh,
                   IOobject::MUST_READ,
                   IOobject::NO_WRITE
         )
);


scalar Sc
(
    RILEMProperties.lookupOrDefault("Sc", 0.7)
);

scalar LEMres
(
    RILEMProperties.lookupOrDefault("LEMres", 24e-6)
);


Foam::string fname
(
    RILEMProperties.lookup("file")
);

scalar pressure
(

    RILEMProperties.lookupOrDefault("pressure", 101325.0)
);


scalar turbThresh
(

    RILEMProperties.lookupOrDefault("turbThresh", 0.0)
);



scalar eqRatio
(

    RILEMProperties.lookupOrDefault("eqRatio",0.0)
);

scalar betaThresh
(

    RILEMProperties.lookupOrDefault("betaThresh",1e-6)
);

scalar betaVarThresh
(

    RILEMProperties.lookupOrDefault("betaVarThresh",1e-5)
);


label NLEM_max
(

    RILEMProperties.lookupOrDefault("NLEM_max", 200)
);

label NLEM_min
(

    RILEMProperties.lookupOrDefault("NLEM_min", 20)
);


//label NLEM_LESLEM
//(
//
//    RILEMProperties.lookupOrDefault("NLEM_LESLEM", 5)
//);


label cBinSize
(

    RILEMProperties.lookupOrDefault("cBinSize", 300)
);

label clusterMax
(
    RILEMProperties.lookupOrDefault("clusterMax", 27)
);

label clusterMin
(
    RILEMProperties.lookupOrDefault("clusterMin", 27)
);


// Has to be non zero or else it seems to not read the dict, possible bug
label nSpark
(
    RILEMProperties.lookupOrDefault("nSpark",1)
);

Info << "nSpark = " << nSpark;


scalar N_ETA
(

    RILEMProperties.lookupOrDefault("N_ETA", 4.0)
);

scalar C_LAMBDA
(

    RILEMProperties.lookupOrDefault("C_LAMBDA", 15.0)
);

scalar bFac
(

    RILEMProperties.lookupOrDefault("bFac", 0.5)
);

scalar bFacCFD
(

    RILEMProperties.lookupOrDefault("bFacCFD", 0.05)
);


bool parallelChem
(

    RILEMProperties.lookupOrDefault("parallelChem", false)
);

bool adjustSuperGrid
(

    RILEMProperties.lookupOrDefault("adjustSuperGrid", false)
);

bool combFlag
(

    RILEMProperties.lookupOrDefault("combFlag", true)
);

bool expRatioFlag
(

    RILEMProperties.lookupOrDefault("expRatioFlag", true)
);

//bool correctExpansion
//(
//
//    RILEMProperties.lookupOrDefault("correctExpansion", false)
//);

bool enhancedMixing
(

    RILEMProperties.lookupOrDefault("enhancedMixing", false)
);


bool advanceLine
(

    RILEMProperties.lookupOrDefault("advanceLine", true)
);


bool integrateScalars
(

    RILEMProperties.lookupOrDefault("integrateScalars", true)
);


//bool writePlanes
//(
//
//    RILEMProperties.lookupOrDefault("writePlanes", false)
//);
//
bool wallDiffusionFlag
(

    RILEMProperties.lookupOrDefault("wallDiffusionFlag", false)
);


bool persistFlag
(

    RILEMProperties.lookupOrDefault("persistFlag", true)
);

bool initdCdT
(

    RILEMProperties.lookupOrDefault("initdCdT", true)
);

word initMatrices
(

    RILEMProperties.lookupOrDefault<word>("initMatrices", "zero")
);

word initEnthalpy
(

    RILEMProperties.lookupOrDefault<word>("initEnthalpy", "total")
);



word stirrInput (

    RILEMProperties.lookup<word>("stirrInput")
);

word initLEMlines (

    RILEMProperties.lookupOrDefault<word>("initLEMlines", "fromY")
);

Info << "mass fraction initialization method:" << initMatrices << endl;


word progSourceMethod (

    RILEMProperties.lookupOrDefault<word>("progSourceMethod", "mapped")
);

Info << "Progress variable source method :" << progSourceMethod << endl;

word PDFmethod (

    RILEMProperties.lookupOrDefault<word>("PDFmethod", "Beta")
);

wordList prodSpecies  //works fine
    (
         RILEMProperties.lookup("prodSpecies")
    );

wordList Cspecies
    (
         RILEMProperties.lookup("Cspecies")
    );
std::vector<std::string> Cspecies_vec(Cspecies.begin(), Cspecies.end());

PtrList<volScalarField> Yprod(prodSpecies.size());
PtrList<volScalarField> Y_LEM(prodSpecies.size()); //LEM averaged mass fractions
PtrList<volScalarField> Yprod_LEM(prodSpecies.size()); //LEM averaged production rates


forAll(Yprod, i){
    Info << "Setting up LEM_Y_ and LEM_Yprod_" << prodSpecies[i]<< Foam::endl;
    Yprod.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "prod_"+prodSpecies[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimless/dimTime, 0.0)
            )
        );

    Y_LEM.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "LEM_Y_"+prodSpecies[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimless, 0.0)
            )
        );

    Yprod_LEM.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "LEM_Yprod_"+prodSpecies[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimless/dimTime, 0.0)
            )
        );
}

Foam::string fuelStream
(
    RILEMProperties.lookup("fuelStream")
);

Foam::string oxStream
(
    RILEMProperties.lookup("oxStream")
);

scalar fuelStreamT
(

    RILEMProperties.lookupOrDefault("fuelStreamT", 500.0)
);

scalar oxStreamT
(

    RILEMProperties.lookupOrDefault("oxStreamT", 300.0)
);

volScalarField& he = thermo.he();
