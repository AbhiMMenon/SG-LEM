if(initLEMlines == "fromY"){
    for(size_t ix = 0; ix < nCells; ix++){
        double avgT = 0, avgP = 0;
        double maxT = 0;
        std::vector<double> avgYi(Y.size(), 0.0);

        lemLines[ix].setPremixedData(eqRatio);
        for(int cellx:cellListSC[ix]){

            // ignition option ?
            maxT = maxT < T[cellx] ? T[cellx]:maxT;
            avgT       += T[cellx];
            avgP       += p[cellx];
            T_LEM[cellx] = T[cellx];
            std::vector<double> y2(oj.ns,-1);

            forAll(Y,i){
                volScalarField& Yi = Y[i];
                y2[i] = Yi[cellx];
                avgYi[i] += Yi[cellx];
            }

        }

        avgT /= cellListSC[ix].size();
        avgP /= cellListSC[ix].size();

        for(int sx = 0; sx < oj.ns; sx++)
            avgYi[sx] /= cellListSC[ix].size();
        lemLines[ix].initializeDataOnLineFromLES(avgT, oj.pressure, avgYi); // ugly name, fix later
        lemLines[ix].setProgVble_premixSp();
        lemLines[ix].setEnergy();
        lemLines[ix].calcScalarC();

    }
}
else if (initLEMlines == "fromPhi"){
    Info << "Phi = " << eqRatio<< Foam::endl;
    for(size_t ix = 0; ix < nCells; ix++){
        lemLines[ix].setPremixedData(eqRatio);
        lemLines[ix].initializeWithPhi(eqRatio);
        double avgC = 0;
        for(int cellx:cellListSC[ix]){
            avgC += Cmean[cellx];
        }
        avgC /= cellListSC[ix].size();
        if(avgC > 0.5){
            lemLines[ix].sparkLine(0);
        }
    }

}
else {
    Info << "\n Unrecognized option : " << initLEMlines << endl;
    Info << "\n exiting.. "  << endl;
    exit(0);
}

// -- ingnite those  LEM lines where c > 0.5


int spark = 0;
