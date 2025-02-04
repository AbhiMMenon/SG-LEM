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
std::vector<LEMLINE> lemLinesCopy;

std::vector<int> resMax(nCells);
std::vector<bool> uniLine(nCells);

// create lemLines, initialize memory and push into the vector. Good to check
// memory leaks etc. here
for(size_t ix = 0; ix < nCells; ix++){

    int localN = domLengthSC[ix]/LEMres;
    localN = std::max(NLEM_min, localN);
    localN = std::min(localN, NLEM_max);
    resMax[ix] = localN;
    LEMLINE le = LEMLINE(resMax[ix], domLengthSC[ix],domLengthSC[ix]/2., file, oj, 200,0);
    le.initZMemory();
    le.bFac = bFac;
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
Info << "All memory initialized\n";
Info << "Initialize data method: "<< initLEMlines << endl;

int fuelInd = std::distance(oj.fuelY.begin(), std::max_element(oj.fuelY.begin(), oj.fuelY.end()));

Info << "Done!" << endl;
int spark = 0;
