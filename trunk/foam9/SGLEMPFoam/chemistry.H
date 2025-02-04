{
    Info<< "Solving chemistry, source term for hs" << endl;

    chemistry.solve
    (
        (runTime.value() - runTime.deltaT().value()),
        runTime.deltaT().value()
    );

////forAll(H_Tot, cellx){
////    std::vector<double> y2(oj.ns,-1);
////    forAll(Y,i){
////        volScalarField& Yi = Y[i];
////        y2[i] = Yi[cellx];
////    }
////    // initialiation
////    chemistrySh_LEM[cellx] = lemLines[0].returnHsource(T[cellx],oj.pressure,y2.data());
////    chemistrySh[cellx] = (chemistrySh_LEM -chemistrySh)/tStep;
////}

    // turbulent time scale
    if (turbulentReaction)
    {
        volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
        volScalarField tc = chemistry.tc();

        // Chalmers PaSR model
        kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
    }
    else
    {
        kappa = 1.0;
    }

   chemistrySh = kappa*chemistry.Sh()();

}
