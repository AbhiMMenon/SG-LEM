/*------------------------------------------------------------------------------
 * initLEMlines.H
 *
 * needs new description for pre-mixed solver
 *
 *-----------------------------------------------------------------------------*/


// --- make an LEM line for each of the superCells, domLength
// --- computed in createSuperGrid.H

Info << "Initializing LEM lines in super-cells ;initializing memory"<< endl;

std::vector<LEMLINE> lemLines;
std::vector<LEMLINE> lemLinesCopy;

std::vector<int> resMax(nCells);
std::vector<bool> uniLine(nCells);
int indO2 = oj.gas.speciesIndex("O2");

// create lemLines, initialize memory and push into the vector. Good to check
// memory leaks etc. here
for(size_t ix = 0; ix < nCells; ix++){

    int localN = domLengthSC[ix]/LEMres;
    localN = std::max(NLEM_min, localN);
    localN = std::min(localN, NLEM_max);
    resMax[ix] = localN;
    
    LEMLINE le = LEMLINE(localN, domLengthSC[ix],domLengthSC[ix]/2., file, oj, 0, cBinSize);
    le.N_ETA = N_ETA;
    le.C_LAMBDA = C_LAMBDA;
    le.initCMemory();
    le.bFac = bFac;

    // this is done here and not in the constructor for better control
    lemLines.push_back(le);
    for(int cellx:cellListSC[ix]){
        nLEMCells[cellx]= resMax[ix];
    }

}
lemLines[0].setC_select(Cspecies_vec);

Info << "All memory initialized\n";
Info << "Initialize data method: "<< initLEMlines << endl;


