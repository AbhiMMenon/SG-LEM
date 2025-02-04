// based on Pierce and Moin (1998)

Info << "Calculating Cvar (PierceMoin)\n";
//scalar delta = 0.0;
//scalarField magSS = Foam::magSqr(fvc::grad(Cmean));
//forAll(magSS, ix){
//    delta = std::cbrt(mesh.V()[ix]);
//    Cvar[ix] = 1./12*delta*delta* magSS[ix];
//}

{
const compressible::LESModel& lesModel =
    Y[0].db().lookupObject<compressible::LESModel>
    (
        momentumTransportModel::typeName
    );

const volScalarField& delta = lesModel.delta();

volScalarField magSS = Foam::magSqr(fvc::grad(Cmean));
 Cvar = (1./12)*Foam::sqr(delta)*(magSS);
}
//Cvar.correctBoundaryConditions();
