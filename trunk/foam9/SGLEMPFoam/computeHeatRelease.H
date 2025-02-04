Info << "Computing heat release rates"<<Foam::endl;
forAll(dYdt,i){
    volScalarField& dYdt_i = dYdt[i];
    volScalarField& Yi =Y[i];
    volScalarField& Yoldi =Yold[i];
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

