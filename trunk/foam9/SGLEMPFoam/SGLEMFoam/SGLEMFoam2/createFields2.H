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

Info<< "Creating field lengthLEM\n" << endl;
volScalarField lengthLEM
(
    IOobject
    (
        "lengthLEM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);

volScalarField Re_Vis
(
    IOobject
    (
        "Re_vis",
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


label NLEM_max
(

    RILEMProperties.lookupOrDefault("NLEM_max", 200)
);

label NLEM_min
(

    RILEMProperties.lookupOrDefault("NLEM_min", 20)
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


label LEMsubSteps
(
    RILEMProperties.lookupOrDefault("LEMsubSteps",1)
);

scalar N_ETA
(

    RILEMProperties.lookupOrDefault("NETA", 4.0)
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

bool enhancedMixing
(

    RILEMProperties.lookupOrDefault("enhancedMixing", true)
);


bool parallelChem
(

    RILEMProperties.lookupOrDefault("parallelChem", false)
);

bool combFlag
(

    RILEMProperties.lookupOrDefault("combFlag", true)
);

bool expRatioFlag
(

    RILEMProperties.lookupOrDefault("expRatioFlag", true)
);

bool advanceLine
(

    RILEMProperties.lookupOrDefault("advanceLine", true)
);



bool wallDiffusionFlag
(

    RILEMProperties.lookupOrDefault("wallDiffusionFlag", false)
);

bool initdCdT
(

    RILEMProperties.lookupOrDefault("initdCdT", true)
);

word initMatrices
(

    RILEMProperties.lookupOrDefault<word>("initMatrices", "zero")
);

word stirrInput (

    RILEMProperties.lookup<word>("stirrInput")
);


word initLEMlines (

    RILEMProperties.lookupOrDefault<word>("initLEMlines", "fromY")
);

Info << "mass fraction initialization method:" << initMatrices << endl;

wordList prodSpecies  //works fine
    (
         RILEMProperties.lookup("prodSpecies")
    );

wordList Cspecies
    (
         RILEMProperties.lookup("Cspecies")
    );
std::vector<std::string> Cspecies_vec(Cspecies.begin(), Cspecies.end());

PtrList<volScalarField> Yprod(prodSpecies.size()); //LEM averaged production rates


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

scalar flamLimitLean
(

    RILEMProperties.lookupOrDefault("flamLimitLean", 0.0)
);



bool adjustSuperGrid = false;
volScalarField& he = thermo.he();
volScalarField k = turbulence->k();
volScalarField nut = turbulence->nut();
volScalarField epsilon = turbulence->epsilon();
