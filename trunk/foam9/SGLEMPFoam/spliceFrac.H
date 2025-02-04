// Determine volume fractions for splicing
std::vector<double>    volFrac(globalFluxSC.size(), 0.0);
for(size_t ix = 0; ix < nCells;ix++){
    double cellVol = cellVolSC[ix];
    for(int facex:outFacesSC_loc[ix]){
        double spliceVol = fabs(globalFluxSC[facex])*tStep;
        double frac =  spliceVol/cellVol;
        volFrac[facex] = frac;
    }
}

// volFrac from a virtual inlet superCell since inlet faces don't have a donor Cell
for(int patchx:inletPatch){
    int owner      = boundaryOwnerSC[patchx][0];
    double cellVol = cellVolSC[owner]; // some reasonable size

    for(int facex:boundaryFaceSC[patchx]){
        double spliceVol = fabs(globalFluxSC[facex])*tStep;
        double frac = spliceVol/ cellVol;
        volFrac[facex] = frac;
    }
}

Info << "Splicing fractions determined" << endl;
