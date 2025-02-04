Info<< "Creating field Re_t\n" << endl;
volScalarField Re_t
(
    IOobject
    (
        "Re_t",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)

);

Info<< "Creating field Karlovitz number\n" << endl;
volScalarField Ka
(
    IOobject
    (
        "Ka",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)

);

Info<< "Creating field Damkölher number\n" << endl;
volScalarField Da
(
    IOobject
    (
        "Da",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)

);

dimensionedScalar Su // Laminar flame speed
(

    RILEMProperties.lookupOrDefault("Su",0.15)
);

Su.dimensions().reset(dimLength/dimTime);
