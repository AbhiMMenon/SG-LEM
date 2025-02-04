if(progSourceMethod=="mapped") cSource = dCdT_mapped*rho;
//else if(progSourceMethod=="explicit") cSource = dCdT_exp*rho;
else {Info << "\n\nIncorrect progsource method" << progSourceMethod << "..exiting\n"; exit(0);}


fvScalarMatrix CEqn
(
    fvm::ddt(rho, Cmean)
  + mvConvection->fvmDiv(phi, Cmean)
//+ thermophysicalTransport->divj(Cmean)
  - fvm::laplacian(turbulence->muEff()/Sc, Cmean)
  ==
  cSource
);

CEqn.relax();
CEqn.solve();

Cmean.max(0.0);
Cmean.min(1.0);


