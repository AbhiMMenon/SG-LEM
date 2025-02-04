/*------------------------------------------------------------------------------
 * inletSplice.H
 *
 * Per patch we find the average Z (usually the Z is constant throughout the
 * patch. Each inlet patch feeds multiple super-grid cells, this is found using
 * boundaryFaceSC[][]. Create one line with enough length and cell volume, and
 * LEM numbers, initialize this to either fuel or oxidiser, spark if needed
 * (heat the entire line to a high temperature) and splice in to the necessary
 * super-cell. This should maintain the non-premixed nature of the fuel unless
 * the patch is given a value other than 1 or 0.
 *-----------------------------------------------------------------------------*/

// -- total mass flow for processor
double mDotFuel = 0.0, mDotOx = 0.0;

// -- for timing

for(size_t ix = 0; ix < inletPatch.size(); ix++){
    int patchx = inletPatch[ix];
    int nFaces = boundaryFaceSC[patchx].size();
    
    double  TFace =std::accumulate(T.boundaryField()[patchx].begin(), T.boundaryField()[patchx].end(),0.0 );
    TFace/=T.boundaryField()[patchx].size();

    std::vector<double> YsPatch(oj.ns,-1);
    forAll(Y, i){
        volScalarField& Yi = Y[i];
        YsPatch[i] = std::accumulate(Yi.boundaryField()[patchx].begin(),Yi.boundaryField()[patchx].end(),0.0);
        YsPatch[i]/=T.boundaryField()[patchx].size();;
    }
   
    
    //-- here we make a large inlet superCell with a large line, refreshed
    //every time step.

    int owner = boundaryOwnerSC[patchx][0]; // typical owner
    int mult  = boundaryFaceListSC[patchx].size();

    // make an LEM line per inlet patch, some typical size
    LEMLINE inletLine = LEMLINE(NLEM_max*mult, domLengthSC[owner]*mult, 1.0, file, oj, 0, 0);

    double typGamma = Gamma[owner];
    inletLine.initializeDataOnLineFromLES(TFace, oj.pressure, YsPatch); // ugly name, fix later
//  inletLine.initializeWithZ(Zline);

    // -- inlet 'pilot' flame
    if(spark < nSpark){
        inletLine.sparkNonPremixed();
        spark++;
        Pout << "\t\t Sparking inlet line " << owner << Foam::endl;}

    for(int jx =0 ; jx < nFaces; jx++){
        std::vector<int> faceList = boundaryFaceListSC[patchx][jx];
        int facex = boundaryFaceSC[patchx][jx];

        GammaList[facex] = typGamma;

        // inlet splicing
        double lenSp = volFrac[facex];
        inletLine.spliceToListLB(lenSp, spliceList[facex]);
    }
}

if(!inletPatch.empty() ){
  //Pout << "Calculated fuel    mass rate (kg/s):" << mDotFuel << endl;
  //Pout << "Calculated Ox      mass rate (kg/s):" << mDotOx   << endl;
  //Pout << "Calculated inlet   mass rate (kg/s):" << mDotOx + mDotFuel << endl;
    Info << "Inlet splicing done.."<< endl;
}
