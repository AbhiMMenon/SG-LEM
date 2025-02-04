Info<< "Creating field Z_LEM\n" << endl;
volScalarField Z_LEM
(
    IOobject
    (
        "Z_LEM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)

);

// --
Info<< "Reading field Zmean\n" << endl;
volScalarField Zmean
(
    IOobject
    (
        "Zmean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Zmean.dimensions().reset(dimless);
//fields.add(Zmean);

Z_LEM = Zmean;

Info<< "Reading field Zvar\n" << endl;
volScalarField Zvar
(
    IOobject
    (
        "Zvar",
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


Foam::List<scalar> gammasCoeff
(

    RILEMProperties.lookup("gammasCoeff")
);

label zBinSize
(

    RILEMProperties.lookupOrDefault("zBinSize", 300)
);



scalar ZMax
(

    RILEMProperties.lookupOrDefault("ZMax", 1.0)
);

scalar flamLimitRich
(

    RILEMProperties.lookupOrDefault("flamLimitRich", 1.0)
);

scalar flamLimitLean
(

    RILEMProperties.lookupOrDefault("flamLimitLean", 0.0)
);

scalar thickeningFactor
(

    RILEMProperties.lookupOrDefault("thickeningFactor", 1.0)
);


bool initZC_LES
(

    RILEMProperties.lookupOrDefault("initZC_LES", false)
);

int spark = 0;
