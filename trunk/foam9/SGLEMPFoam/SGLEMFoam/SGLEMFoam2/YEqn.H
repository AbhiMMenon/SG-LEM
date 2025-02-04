
forAll(Y, i)
{

    
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];
        volScalarField& Y_LEMi = Y_LEM[i];
        volScalarField& Y_LEM_oldi = Y_LEM_old[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          + thermophysicalTransport->divj(Yi)
         ==
            (Y_LEMi - Y_LEM_oldi)/runTime.deltaT()*rho
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);

        Y_LEM_oldi = Y_LEMi;
    }
}

composition.normalise();
