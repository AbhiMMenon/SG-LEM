Info << "Calculating Zvar\n";

{

volScalarField magSS = Foam::magSqr(fvc::grad(Zmean));
 Zvar = (1./12)*Foam::sqr(delta)*(magSS);
}

//forAll(magSS, ix){
//    delta = std::cbrt(mesh.V()[ix]);
//    Zvar[ix] = 1./12.*delta*delta* magSS[ix];
//}

//}

//else{
//    volScalarField Chi = C_Chi * turbulence->epsilon()/turbulence->k() * rho;
//
//    Info << "Calculating Zvar (transport Eqn)\n";
//    fvScalarMatrix ZvarEqn
//    (
//            (
//              fvm::ddt(rho, Zvar)
//            + fvm::div(phi, Zvar)
//            - fvm::laplacian(turbulence->muEff()/Sc, Zvar)
//            - 2*turbulence->muEff()/Sc*sqr(mag(fvc::grad(Zmean)))
//            + fvm::Sp(Chi, Zvar)
//    )
//
//
//    );
//
//    ZvarEqn.relax();
//    ZvarEqn.solve();
//}


