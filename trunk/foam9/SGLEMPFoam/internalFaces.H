/*------------------------------------------------------------------------------
 * internalFaces.H
 *
 * Create lists of owners and neighbours for each (super) internal face
 * (between supercells)
 *
 * Create lists of fine-mesh faces for each (super) internal face (for flux
 * calculation)
 *
 * @author A.Menon 18th sep 2020
 *
 * MODIFIED 7th oct
 *
 * owner_SC and neighbour_SC had to be modified since the original library
 * notation is not consistent with OpenFoam
 *
 * MODIFIED 23rd Oct: owner_SC and neighbour_SC are again fixed by the library,
 * internal flux directions are now checked for consistency at the CFD level by
 * comparing it's superGrid owner with the owner of the face. Based on a
 * discussion with Michael
 *
 *-----------------------------------------------------------------------------*/


std::vector<int> owner_SC(superGrid.owner().begin(),superGrid.owner().end());
std::vector<int> neighbour_SC(superGrid.neighbour().begin(), superGrid.neighbour().end());

// list of all fine-mesh faces in internal face, new internal faces depending upon
std::vector<std::vector<int>> faceList(superGrid.nInternalFaces());
std::vector<int> dirMod(mesh.nInternalFaces(),INT_MAX); // direction modifier

for(int ix =0; ix < mesh.nInternalFaces();ix++){
    int pos = superGrid.faceRestrictAddr()[ix];

    // faceRestrictAddr returns > 0 if the internal face is part of the
    // super-internal Face and < 0 if it's in the interior

    if(pos >= 0){
        faceList[pos].push_back(ix);

        int CFDOwner = mesh.owner()[ix];
        int CFDNeigh = mesh.neighbour()[ix];

        int owner    = fineMesh.child()[CFDOwner];
        int neigh    = fineMesh.child()[CFDNeigh];


        if (owner == owner_SC[pos] && neigh == neighbour_SC[pos])
            dirMod[ix] = 1;
        else if (owner == neighbour_SC[pos] && neigh == owner_SC[pos])
            // reverse flux direction if superGrid owner of the CFD face is not the
            // superGrid owner of the superGrid face
            dirMod[ix] = -1;
        else {
            Pout << "More than two owners for superGrid face " << pos << ".. exiting!!";
            exit(0);
        }
    }


}

//Pout << "Total internal faces created      : " << totFaces << endl;


// internal flux
std::vector<double> fluxSC(totFaces,0.0);
// averge Cmean on SC face
std::vector<double> CmeanSC_face(totFaces,0.0);
