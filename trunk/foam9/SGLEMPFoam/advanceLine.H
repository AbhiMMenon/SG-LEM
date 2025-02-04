/*------------------------------------------------------------------------------
 * advanceLine.H
 *
 * extract turbulence information for all the superCells, set dtStrirr for the
 * LEM lines, advance in time and perform required calculations
 *-----------------------------------------------------------------------------*/

//volScalarField k = turbulence->k();

auto advStart = std::chrono::high_resolution_clock::now();

Info << "Advancing lines..\n";
volScalarField magU = mag(U);
if(advanceLine){
for(size_t ix = 0; ix < nCells;ix++){


    double Re_t_loc = 0.0;
    double velScale = 0.0;
    double C_SG = 0.0;

    if(stirrInput=="resolved"){
      double Vsum  = 0;
      double Umean = 0;
      for(int cellx:cellListSC[ix]){

          Vsum  += mesh.V()[cellx];
          Umean += magU[cellx]* mesh.V()[cellx];
      }

      Umean /= Vsum;

      double vS = 0;
      for(int cellx:cellListSC[ix]){
          vS += abs(magU[cellx] - Umean)* mesh.V()[cellx];
      }
      velScale = vS/Vsum;
      Re_t_loc = velScale*(domLengthSC[ix]*0.5)/lemLines[ix].getMeanKinematicViscosity();
    }
    else if (stirrInput == "SGS") {

        double Vsum  = 0.0;
        double kDotV = 0.0;
        for(int cellx:cellListSC[ix]){

            Vsum  += mesh.V()[cellx];
            kDotV += mesh.V()[cellx] * k[cellx];
            C_SG = Cvar[cellx] > C_SG? Cvar[cellx]: C_SG;
        }

        double kSC = kDotV / Vsum;
        velScale = std::sqrt(2./3. * kSC);
        Re_t_loc = velScale*(domLengthSC[ix]*0.5)/lemLines[ix].getMeanKinematicViscosity();
    }

 // uniLine[ix] = lemLines[ix].getCvar() < 1e-3? false: true;
    uniLine[ix] = C_SG > 0.005? true: false; // the SG cell contains a flame-front

    lemLines[ix].noofEddies = 0;
    if(enhancedMixing){
        lemLines[ix].lt  = domLengthSC[ix]/2.;
        lemLines[ix].eta = lemLines[ix].lt/2.;
        lemLines[ix].eddyEvent();
    }

    //---  threshold
    if(turbThresh > 0.0){
        Re_t_loc = Re_t_loc > turbThresh? turbThresh: Re_t_loc;
    }

    lemLines[ix].lt  = domLengthSC[ix]/2.;
    lemLines[ix].setdtStirrRe(Re_t_loc);


    if(uniLine[ix])
    {
        lemLines[ix].advanceLEMFastPremixed(tStep,  combFlag);
    }
    else{
        lemLines[ix].lengthOfLine();
    }

}
}
else
for(size_t ix = 0; ix < nCells;ix++){
    lemLines[ix].setProgVble_premixSp();
}


auto advStop = std::chrono::high_resolution_clock::now();
auto advTime   = std::chrono::duration_cast<std::chrono::microseconds>
(advStop - advStart);
double advTime_sec = advTime.count()/1e6;
double maxAdvTime;

if(procPatch.size()>0)
    MPI_Reduce(&advTime_sec,&maxAdvTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
else
    maxAdvTime = advTime_sec;


Info <<"Advanced all lines, time (s)  [max all procs] = " << maxAdvTime <<endl;

