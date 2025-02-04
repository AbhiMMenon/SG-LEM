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


for(size_t ix = 0; ix < inletPatch.size(); ix++){
    int patchx = inletPatch[ix];
    int nFaces = boundaryFaceSC[patchx].size();


    // -- find the average T, P and Y for the patch
    std::vector<double> TFace(T.boundaryField()[patchx].begin(),
                                T.boundaryField()[patchx].end());
    std::vector<double> PFace(p.boundaryField()[patchx].begin(),
                                p.boundaryField()[patchx].end() );

    double T_avg = std::accumulate(TFace.begin(), TFace.end(),0.0);
    double P_avg = std::accumulate(PFace.begin(), PFace.end(),0.0);

    int cfdFaces = TFace.size();

    T_avg /= cfdFaces;
    P_avg /= cfdFaces;

    std::vector<double> avgYi(Y.size(), 0.0);
    forAll(Y,i){
        volScalarField& Yi = Y[i];
        std::vector<double> YiFace(Yi.boundaryField()[patchx].begin(),
                Yi.boundaryField()[patchx].end());
        double y_avg = std::accumulate(YiFace.begin(), YiFace.end(),0.0);
        avgYi[i] = y_avg/cfdFaces;
    }


    //-- here we make a large inlet superCell with a large line, refreshed
    //every time step.

    int owner = boundaryOwnerSC[patchx][0]; // typical owner
    int mult  = boundaryFaceListSC[patchx].size();
    mult = 1;

    // make an LEM line per inlet patch, some typical size
    LEMLINE inletLine = LEMLINE(NLEM_max*mult, domLengthSC[owner]*mult, 1.0, file, oj, 0, 0);

    inletLine.initCMemory();
  //inletLine.setC_select(Cspecies_vec);
  //inletLine.bFac = bFac;

    double typGamma = Gamma[owner];
    if(initLEMlines == "fromY")
        inletLine.initializeDataOnLineFromLES(T_avg, P_avg, avgYi);
    else if(initLEMlines == "fromPhi")
        inletLine.initializeWithPhi(eqRatio);
    else exit(0);

    if(spark < nSpark){
        Pout << "Sparking Inlet line\n";
        inletLine.sparkPremixed();
        spark++;

    }

    double lLen = inletLine.lengthOfLine();

    for(int jx =0 ; jx < nFaces; jx++){
        std::vector<int> faceList = boundaryFaceListSC[patchx][jx];
        int facex = boundaryFaceSC[patchx][jx];

        GammaList[facex] = typGamma;

        // inlet splicing
        double lenSp = volFrac[facex];
        inletLine.spliceToListLB(lenSp*lLen, spliceList[facex]);
    }
}

if(!inletPatch.empty() ){
    Info << "Inlet splicing done.."<< endl;
}
