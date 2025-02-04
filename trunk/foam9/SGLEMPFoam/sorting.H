/**-----------------------------------------------------------------------------
 * sorting.H
 *
 * @author A.Menon 29th sep 2020
 *
 * Sort the lists inFacesSC_loc and outFacesSC_loc based on globalFluxSC, in
 * line with S.Arshad splicing method i.e larger flux corresponds to larger
 * displacement of spliced LEM segment
 * -----------------------------------------------------------------------------*/

// sort influx list
for(size_t ix = 0; ix < nCells; ix++){
    std::vector<int> &s = inFacesSC_loc[ix];
    // lambda expression with capture clause
    std::sort(s.begin(), s.end(), [=](int a, int b){
            // sort influxes in descending order (-ve flux)
            return fabs(globalFluxSC[a]) > fabs(globalFluxSC[b]);
    });
}

Info << "inFlux list sorted"<<endl;

// sort outflux list
for(std::vector<int> &s:outFacesSC_loc){
    // lambda expression with capture clause
    std::sort(s.begin(), s.end(), [=](int a, int b){
            // sort out-fluxes in ascending order (+ve flux)
            return fabs(globalFluxSC[a]) < fabs(globalFluxSC[b]);
    });

}

Info << "outFlux list sorted"<<endl;


surfaceScalarField rhoFace(linearInterpolate(rho));
// additional check for flux buildup; (owner, neighbour) test
for(size_t ix = 0; ix < nCells; ix++){
    double phiBuildup = 0.0;


    for(int facex:inFacesSC_loc[ix]){
        phiBuildup += fabs(globalFluxSC[facex]);
    }

    // record phi IN
    for(int cellx:cellListSC[ix])
        phiIn[cellx] = phiBuildup;

    for(int facex:outFacesSC_loc[ix]){
        phiBuildup -= fabs(globalFluxSC[facex]);
    }

    // record total Phi
    for(int cellx:cellListSC[ix])
        phiDiff[cellx] = phiBuildup;

}

