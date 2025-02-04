/*------------------------------------------------------------------------------
 * initLEMlines.H
 *
 *
 *
 *-----------------------------------------------------------------------------*/


// --- make an LEM line for each of the superCells, domLength
// --- computed in createSuperGrid.H

Info << "Initializing LEM lines in super-cells ;initializing memory"<< endl;

std::vector<LEMLINE> lemLines;

std::vector<bool> uniLine(nCells,false);
std::vector<int> resMax(nCells,0);

for(size_t ix = 0; ix < nCells; ix++){

    int localN = domLengthSC[ix]/LEMres;
    localN = std::max(NLEM_min, localN);
    localN = std::min(localN, NLEM_max);
    resMax[ix] = localN;

    LEMLINE le = LEMLINE(resMax[ix], domLengthSC[ix],domLengthSC[ix]/2., file, oj, zBinSize, cBinSize);
    le.bFac = bFac;
    le.initZCMemory();
    le.setZMax(ZMax);
    le.C_LAMBDA = C_LAMBDA;
    le.N_ETA = N_ETA;
    lemLines.push_back(le);

    for(int cellx:cellListSC[ix]){
        nLEMCells[cellx]= resMax[ix];
    }
}
lemLines[0].setC_select(Cspecies_vec);
lemLines[0].setStreams();
lemLines[0].setEquilTY_Z();

oj.F = thickeningFactor;

Info << "All memory initialized\n";
Info << "Initialize data method: "<< initLEMlines << endl;

