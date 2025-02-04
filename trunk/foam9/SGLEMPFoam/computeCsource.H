
//volScalarField k = turbulence->k();
//volScalarField nut = turbulence->nut();
//volScalarField epsilon = turbulence->epsilon();
//scalarField velScale = Foam::sqrt(2./3. * k);



















/*
 *
 * FLAME SPEED MODEL FOR TESTING
 *

 double Sl = 0.1763; // m per second, laminar flame speed

 double Pr = 1.0; // Prandtl number

 scalarField Re_delta = velScale*Delta*rho/turbulence->mu(); // Sub-grid Reynolds number

 scalarField tc = turbulence->alpha()*(1.0/Sl/Sl);

 scalarField Da_delta = Delta/velScale/tc; // sub-grid Damkölher numbea

 // Turbulent flame speed
 scalarField St = Sl*(scalar(1.0) + scalar(0.52)*Foam::sqrt(Re_delta*Pr) * Foam::pow(Da_delta,-0.25)); // Sub-grid Reynolds number

 scalar rho_unb = 1.1874;

 volScalarField CC = rho_unb * St * Foam::mag(fvc::grad(Cmean));
 cSource.internalField() =  Foam::mag(fvc::grad(Cmean))*rho_unb*St;

 volScalarField CC_Rr;
 CC_Rr.internalField() = reaction->R(Y[indO2]);
 /(lemLines[0].Yb_p[indO2] - lemLines[0].Yu_p[indO2]);

 */

