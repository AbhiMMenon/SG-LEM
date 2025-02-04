/* Compute turbulent Reynolds, Karlocitz and Damkölher numbers */

volScalarField k = turbulence->k();
volScalarField nut = turbulence->nut();
volScalarField epsilon = turbulence->epsilon();
scalarField velScale = Foam::sqrt(2./3. * k);


Re_t = Foam::sqrt(2/3.* k)*delta/turbulence->mu()*rho;

epsilon.max(1e-9); // prevent divide by zero, might skew results
volScalarField t_t = k/epsilon;
volScalarField t_k = Foam::sqrt((turbulence->mu()/rho)/epsilon);
volScalarField t_c = (turbulence->mu()/rho)/(Su*Su);


Ka = t_c/t_k;
Da = t_t/t_c;
