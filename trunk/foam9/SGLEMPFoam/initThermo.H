Info << "Initializing H_Tot field"<< Foam::endl;


{
forAll(H_Tot, cellx){

  std::vector<double> y2(oj.ns,-100);

  forAll(Y,i){
      volScalarField& Yi = Y[i];
      y2[i] = Yi[cellx];
  }

  // initialiation
  oj.gas.setState_TPY(298.15,101325, y2.data());
  H_Tot[cellx] = he[cellx] + oj.gas.enthalpy_mass();

}

// -- initialize boundaries for H_Tot
Info << "Initializing H_tot boundaries"<< Foam::endl;


// for he boundaries
std::vector<double> YsPatch(oj.ns,-1);
forAll(H_Tot.boundaryField(), patchx) {
  fvPatchScalarField& pH = const_cast<fvPatchScalarField&> (H_Tot.boundaryField()[patchx]);
  const fvPatchScalarField& ph = he.boundaryField()[patchx];

  forAll(ph, facex){
      forAll(Y, i){
          volScalarField& Yi = Y[i];
          YsPatch[i] = Yi.boundaryField()[patchx][facex];
      }
      oj.gas.setState_TPY(298.15,101325, YsPatch.data());
      pH[facex] = ph[facex] + oj.gas.enthalpy_mass();
  }
}

// for T boundaries
forAll(T.boundaryField(), patchx) {
      fvPatchScalarField& pH = const_cast<fvPatchScalarField&> (H_Tot.boundaryField()[patchx]);
      const fvPatchScalarField& pT = thermo.T().boundaryField()[patchx];
      const fvPatchScalarField& pP = p.boundaryField()[patchx];
      const fvPatchScalarField& ph = he.boundaryField()[patchx];

      forAll(pT, facex){
          forAll(Y, i){
              volScalarField& Yi = Y[i];
              YsPatch[i] = Yi.boundaryField()[patchx][facex];
          }
          oj.gas.setState_TPY(pT[facex],pP[facex], YsPatch.data());
          pH[facex] = oj.gas.enthalpy_mass();
      }
}
H_Tot.correctBoundaryConditions();

}
