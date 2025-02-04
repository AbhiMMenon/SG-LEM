/* 1. Solve H_Eqn
   2. correct interior he
   3. correct fixed T boundaries for H_Tot
   4. correct rest of boundaries for H_Tot
   5. correct he boundaries

 */

{

        // no kinetic energy, low-Mach-number assumptions
        fvScalarMatrix HEqn
        (
            fvm::ddt(rho, H_Tot)
          + mvConvection->fvmDiv(phi, H_Tot)
      //  - fvm::laplacian(turbulence->alphaEff(), H_Tot)
          + thermophysicalTransport->divq(H_Tot)
        );

        HEqn.relax();
        HEqn.solve();

        forAll(H_Tot, celli)
        {
            std::vector<double> y2(oj.ns,-1.0);
            forAll(Y,i){
                volScalarField& Yi = Y[i];
                y2[i] = Yi[celli];
            }

            // -- formation enthalpy
            oj.gas.setState_TPY(298.15, 101325, y2.data());
            double hForm = oj.gas.enthalpy_mass();

            // -- sensible enthalpy for Foam
            he[celli] = H_Tot[celli] - hForm;
        }

      std::vector<double> YsPatch(oj.ns,-1);

      // -- correct boundaries for fixed T boundaries
      forAll(H_Tot.boundaryField(), patchx) {
        fvPatchScalarField& pH = const_cast<fvPatchScalarField&> (H_Tot.boundaryField()[patchx]);
        const fvPatchScalarField& pT = T.boundaryField()[patchx];
        const fvPatchScalarField& pP = p.boundaryField()[patchx];
        const fvPatchScalarField& ph = he.boundaryField()[patchx];
        forAll(pH, facex){
            forAll(Y, i){
                volScalarField& Yi = Y[i];
                YsPatch[i] = Yi.boundaryField()[patchx][facex];
            }
            if(pT.fixesValue()){
                oj.gas.setState_TPY(pT[facex], pP[facex], YsPatch.data());
                pH[facex] = oj.gas.enthalpy_mass();
            }
        }
      }

    H_Tot.correctBoundaryConditions();

    // -- correct boundaries for hs
    forAll(H_Tot.boundaryField(), patchx) {
    const fvPatchScalarField& pH = H_Tot.boundaryField()[patchx];
    const fvPatchScalarField& pT = T.boundaryField()[patchx];
    const fvPatchScalarField& pP = p.boundaryField()[patchx];
    fvPatchScalarField& ph = const_cast<fvPatchScalarField&> (he.boundaryField()[patchx]);
        forAll(pH, facex){
            forAll(Y, i){
                volScalarField& Yi = Y[i];
                YsPatch[i] = Yi.boundaryField()[patchx][facex];
            }
            oj.gas.setState_TPY(298.15, 101325, YsPatch.data());
            ph[facex] = pH[facex] - oj.gas.enthalpy_mass();
        }
    }

    he.correctBoundaryConditions();

      //update Temperatures, viscocity and psi
    thermo.correct();
    Info<< "min/max(T) = " << min(T).value()  << ", " << max(T).value() << endl;
    Info<< "min/max(T_mapped) = " << min(T_mapped).value() << ", " << max(T_mapped).value() << endl;


}

