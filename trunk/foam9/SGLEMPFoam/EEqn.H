{
    // -- update production rates explicit
    Info << "Computing heat release rates"<<Foam::endl;
    forAll(dYdt,i){
        volScalarField& dYdt_i = dYdt[i];
        volScalarField& Yi =Y[i];
        volScalarField& Yoldi =Yold[i];
  //    dYdt_i = fvc::ddt(Yi);
        dYdt_i = (Yi -Yoldi)/runTime.deltaT();
        Yoldi = Yi;

    }

    forAll(dQ_cantera,cellx){
      std::vector<double> dYdt_loc(oj.ns,0);

      forAll(Y,i){
          volScalarField& dYdt_i = dYdt[i];
          dYdt_loc[i] = dYdt_i[cellx];
      }

      dQ_cantera[cellx] = lemLines[0].returnHsource(dYdt_loc.data())*rho[cellx];

    }

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
    //  reaction->Qdot() //kg/m/s^3
        dQ_cantera
      + fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "min/max(T) = " << min(T).value()  << ", " << max(T_mapped).value() << endl;
    Info<< "min/max(T_mapped) = " << min(T_mapped).value() << ", " << max(T_mapped).value() << endl;
}
