/*------------------------------------------------------------------------------
 * internalFlux.H
 *
 * Calculate internal fluxes for all superCell faces for this processor
 * Also, identify which are the outFlux faces and inFlux faces for splicing.
 *
 * @author A.Menon 21th sep 2020
 *
 *
 * LAST EDIT 30th Sep
 *
 * MODIFIED 23rd October: direction modifier for CFD fluxes using dirMod in
 * internalFaces.H 
 *-----------------------------------------------------------------------------*/

// volume flux
// possible that dirMod is not needed if the dot product wit mesh.Sf is removed
surfaceScalarField phiVol(linearInterpolate(U) & mesh.Sf());
surfaceScalarField CSF(linearInterpolate(Cmean));

// clear the lists
for(size_t ix = 0; ix < outFacesSC_loc.size();ix++){
    outFacesSC_loc[ix].clear();
     inFacesSC_loc[ix].clear();
}

for(int ix = 0; ix < totFaces; ix++){

    double flux = 0.0;
    double cLoc = 0.0;

    for(size_t jx = 0; jx < faceList[ix].size(); jx++){

        // faceList contains all the fineGrid faces that make up this superGrid
        // face
        int facex = faceList[ix][jx];
        flux += phiVol[facex] * dirMod[facex];
        cLoc += CSF[facex];
    }
    cLoc /= faceList[ix].size();

    fluxSC[ix] = flux;
    CmeanSC_face[ix] = cLoc;

    // phiVol > 0 indicates flux from owner to neighbour
    if(flux > 0){

        outFacesSC_loc[owner_SC[ix]].    push_back(ix);
         inFacesSC_loc[neighbour_SC[ix]].push_back(ix);
    }
    
    // phiVol < 0 indicates flux from neighbour to owner
    else{

         inFacesSC_loc[owner_SC[ix]].    push_back(ix);
        outFacesSC_loc[neighbour_SC[ix]].push_back(ix);
    }
}

Info << "Supergrid internal fluxes determined" << endl;
