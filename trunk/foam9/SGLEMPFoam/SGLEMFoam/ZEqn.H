//  passive scalar

fvScalarMatrix ZEqn
(
    fvm::ddt(rho, Zmean)
//+ mvConvection->fvmDiv(phi, Zmean)
  + fvm::div(phi, Zmean)
- fvm::laplacian(turbulence->muEff()/Sc, Zmean)
);


ZEqn.relax();
ZEqn.solve();

Zmean.max(0.0);
Zmean.min(ZMax);
