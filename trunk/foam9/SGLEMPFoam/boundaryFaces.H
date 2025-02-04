/*------------------------------------------------------------------------------
 * boundaryFaces.H
 *
 * Create lists of super-cells, fine-cells and faces correspoding to boundary
 * patches for this processor

 * coarseMgMeshLevel does not coarsen boundaries
 *
 * @author A.Menon 14th sep 2020
 *
 *-----------------------------------------------------------------------------*/

std::vector<int> inletPatch;    // list of 'inlet' patches etc...
std::vector<int> fuelInletPatch; // added 3rd Dec
std::vector<int> outletPatch;
std::vector<int> procPatch;
std::vector<int> cyclicPatch;
std::vector<int> wallPatch;

// proc number on the other side of proc patch ; added 13th Oct
std::vector<int> procPatchNeighNo;
std::vector<int> cyclicNeighPatch;

std::vector<std::vector<int>> inletSuperCells; // superCell owner of the (CFD) face
std::vector<std::vector<int>> outletSuperCells;
std::vector<std::vector<int>> procSuperCells;
std::vector<std::vector<int>> procNeighSuperCells; // superCell owner across proc boundary
std::vector<std::vector<int>> wallSuperCells;

std::vector<std::vector<int>> cyclicSuperCells;
std::vector<std::vector<int>> cyclicNeighSuperCells; // superCell owner across cyclic patch

// local (superGrid) face number for influx and out-flux
std::vector<std::vector<int>> outFacesSC_loc(nCells);
std::vector<std::vector<int>>  inFacesSC_loc(nCells);



Info << "Detecting patches, creating super-grid boundaries" <<Foam::endl;

//{{{ Patch number identification
// cycle through all patches and find inlets, outlets and procs and walls
forAll(mesh.boundary(), patch){
    std::string patchName = mesh.boundary()[patch].name();

    if(mesh.boundary()[patch].size() == 0)
        continue; // only the patches found in this processor

    // --- Valid names t.ex "inlet1" or "InletAir" or "FuelInlet" etc..
    if(patchName.find("inlet")!=std::string::npos){
        inletPatch.push_back(patch);
        //Pout << "Inlet     patch found: " << patchName << endl;
    }

    if(patchName.find("Inlet")!=std::string::npos){
        inletPatch.push_back(patch);
        //Pout << "Inlet     patch found: " << patchName << endl;
    }
    
    if(patchName.find("INLET")!=std::string::npos){
        inletPatch.push_back(patch);
        //Pout << "Inlet     patch found: " << patchName << endl;
    }

    if(patchName.find("fuel")!=std::string::npos ||
            patchName.find("Fuel")!=std::string::npos){
        fuelInletPatch.push_back(patch);
        //Pout << "fuelInlet patch found: " << patchName << endl;
    }

    // ---

    // --- Valid names t.ex "outlet" or "2Outlet" or "OutletX" etc..
    if(patchName.find("outlet")!=std::string::npos){
        outletPatch.push_back(patch);
        //Pout << "Outlet    patch found: " << patchName << endl;
    }

    if(patchName.find("Outlet")!=std::string::npos){
        outletPatch.push_back(patch);
        //Pout << "Outlet    patch found: " << patchName << endl;
    }
    
    if(patchName.find("OUTLET")!=std::string::npos){
        outletPatch.push_back(patch);
        //Pout << "Outlet    patch found: " << patchName << endl;
    }
    // ---

    if(mesh.boundary()[patch].type()=="wall"){
        wallPatch.push_back(patch);
        //Pout << "Wall patch found: " << patchName << endl;
    }

    if(isType<processorFvPatch>(mesh.boundary()[patch])){
        auto polyPatch = refCast<const processorPolyPatch>
                            (mesh.boundaryMesh()[patch]);
        int neighProc = polyPatch.neighbProcNo();

        procPatch.push_back(patch);
        procPatchNeighNo.push_back(neighProc);
        //Pout << "Processor patch found: " << patchName << endl;
    }

    if(isType<processorCyclicFvPatch>(mesh.boundary()[patch])){
        auto polyPatch = refCast<const processorCyclicPolyPatch>
                            (mesh.boundaryMesh()[patch]);
        int neighProc = polyPatch.neighbProcNo();

        procPatch.push_back(patch);
        procPatchNeighNo.push_back(neighProc);
        //Pout << "Processor cyclic patch found: " << patchName << endl;
    }

 // if(mesh.boundary()[patch].type()=="cyclic"){  //(isType<cyclic>(mesh.boundary()[patch])){
 //     auto polyPatch = refCast<const cyclicPolyPatch>
 //                         (mesh.boundaryMesh()[patch]);
 //     int neighPatch = polyPatch.neighbPatchID();

 //     cyclicPatch.push_back(patch);
 //     cyclicNeighPatch.push_back(neighPatch);
 //     Pout << "Single proc cyclic patch found: " << patchName << endl;
 // }

}


// concatenate patch lists
std::vector<int> patches = inletPatch;
patches.insert(patches.end(), outletPatch.begin(), outletPatch.end());
patches.insert(patches.end(),   procPatch.begin(),   procPatch.end());

//}}}

//{{{ Patch super-cell identification
inletSuperCells      .resize(inletPatch .size());
outletSuperCells     .resize(outletPatch.size());
procSuperCells       .resize(procPatch  .size());
procNeighSuperCells  .resize(procPatch  .size());
wallSuperCells       .resize(wallPatch  .size());
cyclicSuperCells     .resize(cyclicPatch.size());
cyclicNeighSuperCells.resize(cyclicNeighPatch.size());

// inlet patches
for(size_t ix = 0; ix < inletPatch.size();ix++){
    int patchx = inletPatch[ix];
    // super-cell numbers
    inletSuperCells[ix] = std::vector<int>(superGrid.faceCells(patchx).begin(),
            superGrid.faceCells(patchx).end());
}

// outlet patches
for(size_t ix = 0; ix < outletPatch.size();ix++){
    int patchx = outletPatch[ix];
    outletSuperCells[ix] =
        std::vector<int>(superGrid.faceCells(patchx).begin(),
                superGrid.faceCells(patchx).end());
}

// wall patches
for(size_t ix = 0; ix < wallPatch.size();ix++){
    int patchx = wallPatch[ix];
    wallSuperCells[ix] =
        std::vector<int>(superGrid.faceCells(patchx).begin(),
                superGrid.faceCells(patchx).end());
}

// proc patches
for(size_t ix = 0; ix < procPatch.size();ix++){
    int patchx    = procPatch[ix];
    int neighProc = procPatchNeighNo[ix];


    /*
     * communicate with the neighbouring processor: the super cell owner list
     */

    int nFaces = mesh.boundaryMesh()[patchx].size(); // CFD faces

    std::vector<int> recvOwner(nFaces,0);
    std::vector<int> sendOwner (superGrid.faceCells(patchx).begin(),
            superGrid.faceCells(patchx).end());



    // MPI send recv to neighbouring proc , tag = nFaces
    MPI_Sendrecv( sendOwner.data(), nFaces, MPI_INT, neighProc, nFaces,
                  recvOwner.data(), nFaces, MPI_INT, neighProc, nFaces,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE );


    // update the owner/owner2 lists
         procSuperCells[ix] = sendOwner;
    procNeighSuperCells[ix] = recvOwner;

}

// cyclic single proc patches
for(size_t ix = 0; ix < cyclicPatch.size();ix++){
    int patchx     = cyclicPatch[ix];
    int neighPatch = cyclicNeighPatch[ix];


    // the superCells for this cyclic patch
    cyclicSuperCells[ix] = std::vector<int>(superGrid.faceCells(patchx).begin(),
                superGrid.faceCells(patchx).end());

    // the neighbour superCells for this cyclic patch
    cyclicNeighSuperCells[ix] = std::vector<int>(superGrid.faceCells(neighPatch).begin(),
                superGrid.faceCells(neighPatch).end());

}

//}}}

/*--------------------------------------------------------------------------------/
 * Creating superGrid boundary faces 25th sep
 * Rules:
 *
 *      1. A face must have a unique patch number AND a unique superGrid owner
 *      2. Create a new superGrid face when either superGrid number or
 *      procNeighSuperCells changes,as you scan through CFD faces
 *      3. CFD faces that make up a superGrid boundary face can be identified
 *      using patch number and patch local face number
 *
 *      This would ensure contiguous face agglomeration across processor
 *      boundaries
 *
 * MODIFIED 14th Oct: bug in face proc code produced empty face lists. Fixed
 * now

 * MODIFIED 22th Oct: simplified method for inlet and outlet boundaries, cell
 * by cell approach for proc boundaries should reduce number of proc boundaries

 * MODIFIED 02 Aug 2022: single processor cyclic patches added
 *
 *------------------------------------------------------------------------------*/

// superGrid lists
int totPatch = mesh.boundary().size();
std::vector<std::vector<std::vector<int>>>boundaryFaceListSC(totPatch); // foam numbering
std::vector<std::vector<int>> boundaryOwnerSC (totPatch);               // foam numbering
std::vector<std::vector<int>> boundaryFaceSC  (totPatch);               // super-grid numbering
std::vector<std::vector<int>> boundaryOwner2SC(procPatch.size());
std::vector<std::vector<int>>   cyclicOwner2SC(cyclicPatch.size());

/* Note that superGrid numbering starts with internal faces then inlet faces,
 * outlet faces and then processor faces, wall faces are not here since walls
 * are not involved in splicing
 */

int bFSCcount = 0;

// inlets
for(size_t ix = 0; ix < inletPatch.size(); ix++){

    int patchx = inletPatch[ix];

    // all the boundary faces for each superCell
    std::vector<std::vector<int>> boundaryFaceList(nCells);

    // find the inlet faces for all superCells
    for(size_t jx = 0; jx < inletSuperCells[ix].size();jx++){
        int cellx = inletSuperCells[ix][jx];
        boundaryFaceList[cellx].push_back(jx);
    }

    for(size_t jx = 0; jx < nCells; jx++){
        if(!boundaryFaceList[jx].empty()){
            boundaryFaceListSC[patchx].push_back(boundaryFaceList[jx]);
            boundaryOwnerSC[patchx].push_back(jx);
            boundaryFaceSC [patchx].push_back(bFSCcount + totFaces);
            bFSCcount++;
        }
    }
}

// outlets
for(size_t ix = 0; ix < outletPatch.size(); ix++){

    int patchx = outletPatch[ix];

    // all the boundary faces for each superCell
    std::vector<std::vector<int>> boundaryFaceList(nCells);

    // find the outlet faces for all superCells
    for(size_t jx = 0; jx < outletSuperCells[ix].size();jx++){
        int cellx = outletSuperCells[ix][jx];
        boundaryFaceList[cellx].push_back(jx);
    }

    for(size_t jx = 0; jx < nCells; jx++){
        if(!boundaryFaceList[jx].empty()){
            boundaryFaceListSC[patchx].push_back(boundaryFaceList[jx]);
            boundaryOwnerSC[patchx].push_back(jx);
            boundaryFaceSC [patchx].push_back(bFSCcount + totFaces);
            bFSCcount++;
        }
    }
}

// Walls, does not update boundary face boundaryFaceSC or bFSCcount since it is
// not involved in splicing
for(size_t ix = 0; ix < wallPatch.size(); ix++){

    int patchx = wallPatch[ix];

    // all the boundary faces for each superCell
    std::vector<std::vector<int>> boundaryFaceList(nCells);

    // find the wall faces for all superCells
    for(size_t jx = 0; jx < wallSuperCells[ix].size();jx++){
        int cellx = wallSuperCells[ix][jx];
        boundaryFaceList[cellx].push_back(jx);
    }

    for(size_t jx = 0; jx < nCells; jx++){
        if(!boundaryFaceList[jx].empty()){
            boundaryFaceListSC[patchx].push_back(boundaryFaceList[jx]);
            boundaryOwnerSC[patchx].push_back(jx);
         // boundaryFaceSC [patchx].push_back(bFSCcount + totFaces);
         // bFSCcount++;
        }
    }
}


// Proc and proc cyclic patches {{{
for(size_t ix = 0; ix < procPatch.size(); ix++){

    std::vector<int> tempFace;

    int patchx = procPatch[ix];
    int size   = mesh.boundary()[patchx].size();

    int owner;
    int owner2;

    for(int jx = 0; jx < size-1; jx++){
        // fill up a boundary face
        owner  =      procSuperCells[ix][jx];
        owner2 = procNeighSuperCells[ix][jx];
        tempFace.push_back(jx);

        // store face and start a new face when the owner OR owner-2
        // changes
        if(
              (owner  !=      procSuperCells[ix][jx+1])
                 ||
              (owner2 != procNeighSuperCells[ix][jx+1])
          ){

            // detect overhangs, this is needed due to OpenFoam's face
            // numbering being a bit jumpy, location wise


            // -- make a vector pair for searching
            std::vector<std::pair<int,int>> ownerPair;
            for(size_t kx = 0; kx < boundaryOwnerSC[patchx].size(); kx++){

                std::pair<int, int> p = std::make_pair(boundaryOwnerSC [patchx][kx],
                                                       boundaryOwner2SC[  ix  ][kx]);
                ownerPair.push_back(p);

            }

           // -- search the vector pair for the current combo
           std::pair<int, int> p = std::make_pair(owner, owner2);
           std::vector<std::pair<int, int>>::iterator pos1i
               = std::find(ownerPair.begin(), ownerPair.end(),p);


           // -- if an overhang is detected simply merge it with the appropriate
           // face
           if(pos1i != ownerPair.end()){

               int pos1 = std::distance(ownerPair.begin(), pos1i);

             //Pout << "Boundary Overhang found, owner: " << owner
             //     << " owner2 : " << owner2
             //     << "   CFD faces :";
             //for(int f:tempFace)
             //    Pout << f << " ";
             //Pout << endl;

               std::vector<int>::iterator end =
                   boundaryFaceListSC[patchx][pos1].end();


               // -- merging with older face
               boundaryFaceListSC[patchx][pos1].insert(end,
                       tempFace.begin(), tempFace.end());
               tempFace.clear();
               continue;
           }


            // -- create new face if not an overhang
            boundaryOwnerSC   [patchx].push_back(owner);
            boundaryOwner2SC  [ix]    .push_back(owner2);
            boundaryFaceSC    [patchx].push_back(bFSCcount + totFaces);
            boundaryFaceListSC[patchx].push_back(tempFace);
            bFSCcount++;

            tempFace.clear();
        }
    }


    // final CFD face
    tempFace.push_back(size-1);

    // repeat everything fort the final face since neither owner nor owner 2
    // will change

    // -- make a vector pair for searching
    std::vector<std::pair<int,int>> ownerPair;
    for(size_t kx = 0; kx < boundaryOwnerSC[patchx].size(); kx++){

        std::pair<int, int> p = std::make_pair(boundaryOwnerSC [patchx][kx],
                                               boundaryOwner2SC[  ix  ][kx]);
        ownerPair.push_back(p);

    }

    // -- search the vector pair for the current combo

    std::pair<int, int> p = std::make_pair(owner, owner2);
    std::vector<std::pair<int, int>>::iterator pos1i
        = std::find(ownerPair.begin(), ownerPair.end(),p);

    if(pos1i != ownerPair.end()){

        int pos1 = std::distance(ownerPair.begin(), pos1i);
      //Pout << "Boundary Overhang found, owner: " << owner
      //     << " owner2 : " << owner2
      //     << "   CFD faces :";
      //for(int f:tempFace)
      //    Pout << f << " ";
      //Pout << endl;

        std::vector<int>::iterator end =
            boundaryFaceListSC[patchx][pos1].end();


        // -- merging with older face
        boundaryFaceListSC[patchx][pos1].insert(end,
                tempFace.begin(), tempFace.end());
        tempFace.clear();
    }
    else{

        // -- create new face
        boundaryOwnerSC   [patchx].push_back(owner);
        boundaryOwner2SC  [ix]    .push_back(owner2);
        boundaryFaceListSC[patchx].push_back(tempFace);
        boundaryFaceSC    [patchx].push_back(bFSCcount + totFaces);
        bFSCcount++;
    }


    // some outputs for debugging
  //Pout << " faces created for [" << procPatchNeighNo[ix] << "]"
  //     << ": "<< nProcFaces
  //     << endl;

  //if(Pstream::myProcNo() == 5 && procPatchNeighNo[ix] == 2){

  //    Pout << "BOUNDARY FACE TEST----"<< endl << endl;
  //    for(int jx = 0; jx < nProcFaces; jx++){
  //        std::vector<int> faces = boundaryFaceListSC[patchx][jx];
  //        Pout << "superFace no:" << jx
  //             << " |owner    :" << boundaryOwnerSC[patchx][jx]
  //             << " |owner2   :" << boundaryOwner2SC[ix][jx]
  //             << " |       ----> ";
  //        for(int i:faces)
  //            Pout << i <<"  ";
  //        Pout << endl;
  //    }

  //    for(int jx = 0; jx < size; jx++){
  //        Pout << "face combo "   << jx  << ":  "
  //             << procSuperCells[ix][jx] << "   "
  //                     << procNeighSuperCells[ix][jx] << endl;
  //    }

  //}
}
// }}}

//{{{ single processor cyclic
for(size_t ix = 0; ix < cyclicPatch.size(); ix++){

    std::vector<int> tempFace;

    int patchx = cyclicPatch[ix];
    int size   = mesh.boundary()[patchx].size();

    int owner;
    int owner2;

    for(int jx = 0; jx < size-1; jx++){
        // fill up a boundary face
        owner  =      cyclicSuperCells[ix][jx];
        owner2 = cyclicNeighSuperCells[ix][jx];
        tempFace.push_back(jx);

        // store face and start a new face when the owner OR owner-2
        // changes
        if(
              (owner  !=      cyclicSuperCells[ix][jx+1])
                 ||
              (owner2 != cyclicNeighSuperCells[ix][jx+1])
          ){

            // detect overhangs, this is needed due to OpenFoam's face
            // numbering being a bit jumpy, location wise


            // -- make a vector pair for searching
            std::vector<std::pair<int,int>> ownerPair;
            for(size_t kx = 0; kx < boundaryOwnerSC[patchx].size(); kx++){

                std::pair<int, int> p = std::make_pair(boundaryOwnerSC [patchx][kx],
                                                       cyclicOwner2SC  [  ix  ][kx]);
                ownerPair.push_back(p);

            }

           // -- search the vector pair for the current combo
           std::pair<int, int> p = std::make_pair(owner, owner2);
           std::vector<std::pair<int, int>>::iterator pos1i
               = std::find(ownerPair.begin(), ownerPair.end(),p);


           // -- if an overhang is detected simply merge it with the appropriate
           // face
           if(pos1i != ownerPair.end()){

               int pos1 = std::distance(ownerPair.begin(), pos1i);

             //Pout << "Boundary Overhang found, owner: " << owner
             //     << " owner2 : " << owner2
             //     << "   CFD faces :";
             //for(int f:tempFace)
             //    Pout << f << " ";
             //Pout << endl;

               std::vector<int>::iterator end =
                   boundaryFaceListSC[patchx][pos1].end();


               // -- merging with older face
               boundaryFaceListSC[patchx][pos1].insert(end,
                       tempFace.begin(), tempFace.end());
               tempFace.clear();
               continue;
           }


            // -- create new face if not an overhang
            boundaryOwnerSC   [patchx].push_back(owner);
            cyclicOwner2SC    [ix]    .push_back(owner2);
            boundaryFaceSC    [patchx].push_back(bFSCcount + totFaces);
            boundaryFaceListSC[patchx].push_back(tempFace);
            bFSCcount++;

            tempFace.clear();
        }
    }


    // final CFD face
    tempFace.push_back(size-1);

    // repeat everything fort the final face since neither owner nor owner 2
    // will change

    // -- make a vector pair for searching
    std::vector<std::pair<int,int>> ownerPair;
    for(size_t kx = 0; kx < boundaryOwnerSC[patchx].size(); kx++){

        std::pair<int, int> p = std::make_pair(boundaryOwnerSC [patchx][kx],
                                                 cyclicOwner2SC[  ix  ][kx]);
        ownerPair.push_back(p);

    }

    // -- search the vector pair for the current combo

    std::pair<int, int> p = std::make_pair(owner, owner2);
    std::vector<std::pair<int, int>>::iterator pos1i
        = std::find(ownerPair.begin(), ownerPair.end(),p);

    if(pos1i != ownerPair.end()){

        int pos1 = std::distance(ownerPair.begin(), pos1i);
      //Pout << "Boundary Overhang found, owner: " << owner
      //     << " owner2 : " << owner2
      //     << "   CFD faces :";
      //for(int f:tempFace)
      //    Pout << f << " ";
      //Pout << endl;

        std::vector<int>::iterator end =
            boundaryFaceListSC[patchx][pos1].end();


        // -- merging with older face
        boundaryFaceListSC[patchx][pos1].insert(end,
                tempFace.begin(), tempFace.end());
        tempFace.clear();
    }
    else{

        // -- create new face
        boundaryOwnerSC   [patchx].push_back(owner);
        cyclicOwner2SC    [ix]    .push_back(owner2);
        boundaryFaceListSC[patchx].push_back(tempFace);
        boundaryFaceSC    [patchx].push_back(bFSCcount + totFaces);
        bFSCcount++;
    }


    // some outputs for debugging
    int nProcFaces = boundaryFaceSC[patchx].size();
  //Pout << " faces created for [" << cyclicNeighPatch[ix] << "]"
  //     << ": "<< nProcFaces
  //     << endl;
}//}}}

if(procPatch.size() >0) MPI_Barrier(MPI_COMM_WORLD);

Info << "-----------------------DONE---------------------------" << endl;
Info << endl;

//Pout << "SuperGrid boundary faces created: " << bFSCcount << endl;

// -- initialize fluxes
std::vector<std::vector<double>> mDotFuelSC(inletPatch.size());
std::vector<std::vector<double>> mDotOxSC  (inletPatch.size());
std::vector<std::vector<double>> mDotOutSC (outletPatch.size());


// create splice list
std::vector<std::list<CellData>>  spliceList(totFaces + bFSCcount);
std::vector<double> GammaList(totFaces + bFSCcount);
