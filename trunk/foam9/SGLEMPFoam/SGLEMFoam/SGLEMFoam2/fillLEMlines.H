if(initLEMlines == "fromY"){
    for(size_t ix = 0; ix < nCells; ix++){
        double avgT = 0, avgP = 0;
        double maxT = 0;
        std::vector<double> avgYi(oj.ns, -0.0);

        for(int cellx:cellListSC[ix]){

            // ignition option ?
            maxT = maxT < T[cellx] ? T[cellx]:maxT;
            avgT       += T[cellx];
            avgP       += p[cellx];
        //  T_LEM[cellx] = T[cellx];
            std::vector<double> y2(oj.ns,-1);

            forAll(Y,i){
                volScalarField& Yi = Y[i];
                y2[i] = Yi[cellx];
                avgYi[i] += Yi[cellx];
            }

        }

        avgT /= cellListSC[ix].size();
        avgP /= cellListSC[ix].size();

        for(size_t sx = 0; sx < oj.ns; sx++)
            avgYi[sx] /= cellListSC[ix].size();
        lemLines[ix].initializeDataOnLineFromLES(avgT, oj.pressure, avgYi);
        lemLines[ix].setBilgerMixFrac();
    }
}


else{
    Info << "\n Unrecognized option : " << initLEMlines << endl;
    Info << "\n exiting.. "  << endl;
    exit(0);
}

Info << "Done!" << endl;
