 /*
 * Depopulate lists to in-splice using the ordered face list inFacesSC_loc
 *
 * 7th Oct
 * Reorganized 13th Oct
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

    lemLines[ix].makeEquidistantGridafterSplicing(resMax[ix], 1e-9);

    // for visualization
    int    nCell   = lemLines[ix].cells.size();
    double lemMass = lemLines[ix].massOnLine();
    double lenLEM  = lemLines[ix].lengthOfLine();
//  double Zline   = lemLines[ix].getMeanMixFrac();
    double nEd     = lemLines[ix].dtStirr>0.0? 1./lemLines[ix].dtStirr:1.0;
    double rat     = nSplice/nEd;

    std::vector<double> vizLine = lemLines[ix].getMeanS(pList);
    for(int cellx:cellListSC[ix]){
        nLEMCells[cellx]= nCell;
        massLEM[cellx]  = lemMass;
        lengthLEM[cellx]= lenLEM;
        nEddies[cellx]  = nEd;
        nSplicings[cellx]  = nSplice;
        nRatio[cellx]   = rat;

        // -- LEM averaged mass fractions, pList for both production rates and
        // LEM averaged mass fractions
        forAll(Yprod,i){
            volScalarField& Yprod_LEM_loc = Yprod[i];
            Yprod_LEM_loc[cellx] = vizLine[IY+i+oj.ns];
        }


    }

}



// clear outlets
for(int patchx:outletPatch){
    for(int facex:boundaryFaceSC[patchx]){
        spliceList[facex].clear();
    }
}

Info << "All splicing completed\n\n\n";
