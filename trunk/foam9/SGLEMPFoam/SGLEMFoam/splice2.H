 /*
 * Depopulate lists to in-splice using the ordered face list inFacesSC_loc
 *
 * 7th Oct
 * Reorganized 13th Oct
 ------------------------------------------------------------------------------*/

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
    double lemMass = lemLines[ix].lengthOfLine()/domLengthSC[ix];
    double Zline   = lemLines[ix].getMeanMixFrac();
    double nEd     = lemLines[ix].noofEddies;

    std::vector<double> vizLine = lemLines[ix].getMeanS(pList);
    for(int cellx:cellListSC[ix]){
        nLEMCells[cellx]= nCell;
        massLEM[cellx]  = lemMass;
        Z_LEM[cellx]     = Zline;
        T_LEM[cellx]    = vizLine[IT];
        C_LEM[cellx]    = vizLine[IC];
        Csource_LEM[cellx]  = vizLine[ITheta]; // note that ITheta and IC switch places here
        nEddies[cellx]  = nEd;
        nSplicings[cellx]  = nSplice;

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

}



// clear outlets
for(int patchx:outletPatch){
    for(int facex:boundaryFaceSC[patchx]){
        spliceList[facex].clear();
    }
}

Info << "All splicing completed\n\n\n";
