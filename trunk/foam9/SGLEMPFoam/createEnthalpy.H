wordList hBoundaryTypes(he.boundaryField().size());

const volScalarField::Boundary& tbf = T.boundaryField();

hBoundaryTypes = tbf.types();

forAll(tbf, patchi)
{
    if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
    {
        hBoundaryTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
    else if
    (
        isA<zeroGradientFvPatchScalarField>(tbf[patchi])
  // || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
    )
    {
        hBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
  //else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
  //{
  //    hBoundaryTypes[patchi] = mixedFvPatchScalarField::typeName;
  //}
}

volScalarField H_Tot
(
    IOobject
    (
        "H_Tot",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimEnergy/dimMass,
    hBoundaryTypes
 // thermo.he().boundaryField().types()
);

H_Tot = thermo.he();
fields.add(H_Tot);

Info << "Initializing H_Tot field"<< Foam::endl;
forAll(H_Tot, cellx){

  std::vector<double> y2(oj.ns,-100);

  forAll(Y,i){
      volScalarField& Yi = Y[i];
      y2[i] = Yi[cellx];
  }

  // initialiation
  if(initEnthalpy=="total"){
      oj.gas->setState_TPY(T_mapped[cellx],oj.pressure, y2.data());
      rho[cellx] = RHO_mapped[cellx];
      H_Tot[cellx] = oj.gas->enthalpy_mass();

      oj.gas->setState_TPY(298.15,101325, y2.data());
      he[cellx] = H_Tot[cellx] - oj.gas->enthalpy_mass();
  }
  else if (initEnthalpy == "sensible"){
      oj.gas->setState_TPY(298.15,101325, y2.data());
      H_Tot[cellx] = he[cellx] + oj.gas->enthalpy_mass();
  }
  else{Info<< "Incorrect option" << initEnthalpy<< ".. exitinig\n"; exit(0);}

}

thermo.correct();
he.correctBoundaryConditions();

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
      oj.gas->setState_TPY(298.15,101325, YsPatch.data());
      pH[facex] = ph[facex] + oj.gas->enthalpy_mass();
  }
}

// for T boundaries
forAll(T.boundaryField(), patchx) {
      fvPatchScalarField& pH = const_cast<fvPatchScalarField&> (H_Tot.boundaryField()[patchx]);
      fvPatchScalarField& ph = const_cast<fvPatchScalarField&> (he.boundaryField()[patchx]);
      const fvPatchScalarField& pT = thermo.T().boundaryField()[patchx];
      const fvPatchScalarField& pP = p.boundaryField()[patchx];

      if (pT.fixesValue())
      forAll(pT, facex){
          forAll(Y, i){
              volScalarField& Yi = Y[i];
              YsPatch[i] = Yi.boundaryField()[patchx][facex];
          }
          oj.gas->setState_TPY(pT[facex],pP[facex], YsPatch.data());
          pH[facex] = oj.gas->enthalpy_mass();

          oj.gas->setState_TPY(298.15,101325, YsPatch.data());
          ph[facex] = pH[facex] - oj.gas->enthalpy_mass();
      }
}

H_Tot.correctBoundaryConditions();
