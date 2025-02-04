 /*
 * Depopulate lists to in-splice using the ordered face list inFacesSC_loc
 *
 ------------------------------------------------------------------------------*/

// de-populate splice list INSPLICE

for(size_t ix = 0; ix < nCells; ix++){
    int nSplice = 0;

    if(expRatioFlag){
        for(int facex:inFacesSC_loc[ix]){

            double expRatio = GammaList[facex] / Gamma[ix];
            lemLines[ix].spliceFromList(spliceList[facex], expRatio);
            nSplice++;
        }
    }
    else{
        for(int facex:inFacesSC_loc[ix]){
            lemLines[ix].spliceFromList(spliceList[facex]);
            nSplice++;
        }
    }

    lemLines[ix].setProgVble_premixSp();
    lemLines[ix].makeEquidistantGridafterSplicing(resMax[ix] , 1e-9);
  //if(!uniLine[ix])
  //    lemLines[ix].makeEquidistantGridafterSplicing(2 , 1e-9);
  //else
  //    lemLines[ix].makeEquidistantGridafterSplicing(resMax[ix] , 1e-9);

    // for visualization
    int    nCell   = lemLines[ix].cells.size();
    double lemMass = uniLine[ix]? 1.0: 0.0;
    double nEd     = lemLines[ix].noofEddies;

    std::vector<double> vizLine = lemLines[ix].getMeanS(pList);
    double fullness = lemLines[ix].getCFullness();

    for(int cellx:cellListSC[ix]){
        nLEMCells[cellx]= nCell;
        massLEM[cellx]  = lemMass;
        T_LEM[cellx]    = vizLine[IT];
        C_LEM[cellx]    = vizLine[IC];
        Csource_LEM[cellx]  = vizLine[ITheta]; // note that ITheta and IC switch places here
        nEddies[cellx]  = nEd;
        nSplicings[cellx]  = nSplice;
        intHoles[cellx] = fullness;

        // -- LEM averaged mass fractions, pList for both production rates and
        // LEM averaged mass fractions
        forAll(Yprod,i){
            int s = pList[i];
            volScalarField& Y_LEM_loc = Y_LEM[i];
            volScalarField& Yprod_LEM_loc = Yprod_LEM[i];
            Y_LEM_loc[cellx] = vizLine[IY+s];
            Yprod_LEM_loc[cellx] = vizLine[IY+i+oj.ns];
        }


    }

    // -- check for accidental data write due to pointer issues
    if(lemLines[ix].N_ETA != N_ETA || lemLines[ix].C_LAMBDA != C_LAMBDA){
        Pout << "N_ETA after Splice2 changed from = "<< N_ETA << "to "  << lemLines[ix].N_ETA << "\n";
        Pout << "C_LAMBDA after Splice2 changed from = "<< C_LAMBDA << "to "  << lemLines[ix].C_LAMBDA << "\n";
    }
}



for(int patchx:outletPatch){
    for(int facex:boundaryFaceSC[patchx]){
        spliceList[facex].clear();
    }
}


Info << "All splicing completed\n\n\n";
