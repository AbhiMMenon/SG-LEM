/*------------------------------------------------------------------------------
 * integrateScalars.H
 *
 *
 * intFAIL codes, can visualize with paraview
 *
 * 0 -> delta PDF used for skip
 * 1 -> succesfull beta PDF integration
 * 2 -> beta PDF failure, boost library
 * 3 -> incorrect integrated value, based on temperature
 *
 * A.Menon
 * *-----------------------------------------------------------------------------*/


Info << "Mapping scalars..\n";

auto chemDiffStart = std::chrono::high_resolution_clock::now();
if(integrateScalars){
//#pragma omp parallel for schedule(dynamic) default(shared)
    for(size_t ix = 0; ix < nCells;ix++){
        lemLines[ix].fillZC_biased();

        for(int cellx:cellListSC[ix]){
            double Z  = Zmean[cellx];
            double Zv = Zvar[cellx];
            double Cm = Cmean[cellx];
            double Cv = Cvar[cellx];
            double ZvMin = betaVarThresh;

            // solution vector
            std::vector<double> SOL(NSCALAR + oj.ns,-2);

            if (Zv < ZvMin) //&& Cv < ZvMin )
         // if (true) //&& Cv < ZvMin )
            {
                // -- assume delta PDF
                SOL=lemLines[ix].getInterpS(Z,Cm,101325,false);
                intFAIL[cellx] = 4;
            }
            else if (PDFmethod == "Beta"){

                Z = std::max(Z,betaThresh);
                Z = std::min(Z,1.-betaThresh);
                Cm = std::max(Cm, betaThresh);
                Cm = std::min(Cm,1.-betaThresh);

                SOL = lemLines[ix].getPDFMeanSBV(Z, Zv,Cm,Cv);
                intFAIL[cellx] = 1;

            }
            else if (PDFmethod == "TopHat"){

                Z = std::max(Z,betaThresh);
                Z = std::min(Z,1.-betaThresh);
                Cm = std::max(Cm, betaThresh);
                Cm = std::min(Cm,1.-betaThresh);

                SOL = lemLines[ix].getPDFMeanS2_TH(Z, Zv,Cm);
                intFAIL[cellx] = 2;

            }
            else{ Info << "Unkown option : " << PDFmethod<<", exiting";  exit(0);}

            if (SOL[IT]  < 100){
                   intFAIL[cellx] = 3;
                   continue;
            }

            // -- Assign values
            H_mapped[cellx]   = SOL[Ih];
            RHO_mapped[cellx] = SOL[IRho];
            dCdT_mapped[cellx] = SOL[IC];
            T_mapped[cellx]       = SOL[IT];
            C_mapped[cellx]       = SOL[IC2];
      //    dCdT_exp[cellx]    = SOL[IC2];

            // -- update production rates
            forAll(Yprod,i){
                volScalarField& Ypi = Yprod[i];
                Ypi[cellx] = SOL[IY+i+oj.ns];

            }

            // -- update species
            forAll(Y,i){
                if(composition.solve(i)){
                volScalarField& Yi = Y[i];
                Yi[cellx] = SOL[IY + i]*bFacCFD + (1.-bFacCFD)*Yi[cellx];
                }
            }

        }
    }
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

forAll(Y,i){
    volScalarField& Yi = Y[i];
    Yi.correctBoundaryConditions();
}

forAll(Yprod,i){
    // just for visualisation
    Y_LEM[i].correctBoundaryConditions();
}
composition.normalise();

Info <<"Done Mapping, time (s) [max all procs] = " << maxMapTime <<endl;
