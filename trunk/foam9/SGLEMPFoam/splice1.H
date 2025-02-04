/*------------------------------------------------------------------------------
 * internalSplice.H
 *
 * Populate face lists out-spliced LEM segments following the face number order
 * sorted by Phi, outFacesSC_loc.  
 *
 * 7th Oct,
 * 
 * Reorganized 13th Oct
 *
 * MODIFIED 30th Nov, length fraction based splicing
 ------------------------------------------------------------------------------*/


// populate splice list OUTSPLICE
// *combine with previous loop*
for(size_t ix = 0; ix < nCells; ix++){
    double lLen = lemLines[ix].lengthOfLine();
    size_t cell1 = cellListSC[ix][0];
    for(int facex:outFacesSC_loc[ix]){
        double lenSp  = volFrac[facex]; // lengthFraction
        if((lenSp>1e-9)&& (lenSp < 1.0-1e-9)){
            lemLines[ix].spliceToListLB(lenSp*lLen,  spliceList[facex]);
            // send Gamma 16th Nov
            GammaList[facex] = Gamma[ix];
        }
    }
}

Info << "Splicing list populated" << endl;
