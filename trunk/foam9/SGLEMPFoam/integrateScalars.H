/*------------------------------------------------------------------------------
 * integrateScalars.H
 *
 * Evaluates the integral of conditional PDFs based on Zmean, Zvar and Cmean.
 * Temperature can either be updated using the energy equation or directly from
 * integration, the same for dQ.
 *
 *
 * In case the dQ is 0 from the line (due to poor statistical fidelity), it is
 * evaluated using the standard option (PaSR or WSR) using the equilibrium
 * solution which was used to initialize the conditional matrix on the line.
 *
 * A.Menon
 * *-----------------------------------------------------------------------------*/


Info << "Mapping scalars..\n";

auto chemDiffStart = std::chrono::high_resolution_clock::now();
if(integrateScalars ){
//#pragma omp parallel for schedule(dynamic) default(shared)
    for(size_t ix = 0; ix < nCells;ix++){

        for(int cellx:cellListSC[ix]){
            double C  = Cmean[cellx];
            double Cv  = Cvar[cellx];
            double CvMax  =C*(1.-C);
            double CvMin  =betaVarThresh;
            // solution vector
            std::vector<double> SOL(NSCALAR + nSpecies+pList.size(),0.0);

            if(Cv < CvMin || C < betaThresh || (C > (1. -betaThresh))) {
                SOL = lemLines[ix].getInterpSC(C,oj.pressure, pList);
                intFAIL[cellx] = 0;
            }
            else if (PDFmethod == "Beta"){
                Cv = Cv >= CvMax? 0.999*CvMax:Cv;
                SOL = lemLines[ix].getPDFMeanSC_Beta(C,Cv,pList);
                intFAIL[cellx] = 2;
            }
            else if (PDFmethod == "TopHat"){

                C = std::max(C, betaThresh);
                C = std::min(C,1.-betaThresh);
                SOL = lemLines[ix].getPDFMeanSC_TH(C,Cv,pList);
                intFAIL[cellx] = 2;

            }


            else {Info << "PDFmethod error!\n"; exit(0);}

            if (SOL[IT]  < 100){
                intFAIL[cellx] = -10;
                continue;
            }

            // -- Assign values
            H_mapped[cellx]    = SOL[Ih];
            RHO_mapped[cellx]  = SOL[IRho];
            dQ_cantera[cellx]  = SOL[ITheta];
            T_mapped[cellx]    = SOL[IT];
            dCdT_mapped[cellx]       = SOL[IC];

            // -- update production rates
            forAll(Yprod,i){
                volScalarField& Ypi = Yprod[i];
                Ypi[cellx] = SOL[IY+i+oj.ns];

            }

            // -- update species
            forAll(Y,i){
                volScalarField& Yi = Y[i];
                Yi[cellx] = SOL[IY + i]*bFacCFD + (1.-bFacCFD)*Yi[cellx];
            }

        }
        if(!persistFlag && PDFmethod != "Beta")
                lemLines[ix].clearCmats();

    }

}
forAll(Y,i){
    volScalarField& Yi = Y[i];
    Yi.correctBoundaryConditions();
}

auto chemDiffStop = std::chrono::high_resolution_clock::now();
auto mapTime   = std::chrono::duration_cast<std::chrono::microseconds>
(chemDiffStop - chemDiffStart);
double mapTime_sec = mapTime.count()/1e6;
double maxMapTime;

if(procPatch.size()>0)
    MPI_Reduce(&mapTime_sec,&maxMapTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
else
    maxMapTime = mapTime_sec;

composition.normalise();

Info <<"Done Mapping, time (s) [max all procs] = " << maxMapTime <<endl;
