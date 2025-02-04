Info<< "Creating field ltVis\n" << endl;
volScalarField ltVis
(
    IOobject
    (
        "ltVis",
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

