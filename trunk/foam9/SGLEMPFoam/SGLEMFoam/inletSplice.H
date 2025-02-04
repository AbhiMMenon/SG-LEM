/*------------------------------------------------------------------------------
 * inletSplice.H
 *-----------------------------------------------------------------------------*/

// -- total mass flow for processor
double mDotFuel = 0.0, mDotOx = 0.0;

// -- for timing

for(size_t ix = 0; ix < inletPatch.size(); ix++){
    int patchx = inletPatch[ix];
    int nFaces = boundaryFaceSC[patchx].size();
    std::vector<double> phiFace(phiVol.boundaryField()[patchx].begin(),
                                phiVol.boundaryField()[patchx].end() );
    std::vector<double>    ZFace(Zmean.boundaryField()[patchx].begin(),
                                 Zmean.boundaryField()[patchx].end() );
    std::vector<double>    rhoFace(rho.boundaryField()[patchx].begin(),
                                   rho.boundaryField()[patchx].end() );
    // find average ZFace
    double Zline = std::accumulate(ZFace.begin(),ZFace.end(),0.0);
    Zline = Zline/ZFace.size();


    int owner = boundaryOwnerSC[patchx][0]; // typical owner
    int mult  = boundaryFaceListSC[patchx].size();
        mult  = 1.0;


    double typGamma = Gamma[owner];


    for(int jx =0 ; jx < nFaces; jx++){
        // make an LEM line per inlet patch, some typical size
        // -- inlet 'pilot' flame
        LEMLINE inletLine = LEMLINE(resMax[owner]*mult, domLengthSC[owner]*mult, 1.0, file, oj, 0, 0);
        inletLine.initializeWithZ(Zline);
        double lLen = inletLine.lengthOfLine();
        if(spark < nSpark && Zline > 0.0){
            inletLine.injectRadicals();
            spark++;
            Pout << "\t\t Injecting radicals on line " << owner << Foam::endl;}
        
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
