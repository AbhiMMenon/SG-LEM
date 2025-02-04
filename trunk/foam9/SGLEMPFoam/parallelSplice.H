/**--------------------------------------------------------------------------------/
 * @file parallelSplice.H
 *
 * 1. Send vectors to the neighbour processors about how many LEM cells per face
 * one each processor will be sent
 * 2. Receive the same info from neighbour processors
 *
 * 3. Prepare vectors for sending and receiving cell data (dx, rho, T, p,M, Ys)
 * {5 + nSpecies}
 *
 * 4. Comm cell data vector with neighbour
 * 5. Parse received data and populate spliceList at the proc boundary faces
 *
 *
 * @author A.Menon 13th Oct
 *
 * MODIFIED 15th oct: Trying 1 MPI operation per patch
 *------------------------------------------------------------------------------*/



Info << "Parallel splicing..." << endl;
for(size_t ix = 0; ix < procPatch.size(); ix++){

    int patchx    = procPatch[ix];
    int neighProc = procPatchNeighNo[ix];
    int nFaces = boundaryFaceSC[patchx].size();

    std::vector<int> _sendSizes(nFaces,0);
    std::vector<int> _recvSizes(nFaces);

    // -- count no. of LEM cells
    for(int jx = 0; jx < nFaces; jx++){
        int facex = boundaryFaceSC[patchx][jx];
        _sendSizes[jx] = spliceList[facex].size();

    }

    // MPI send recv tag = nFaces
    // Comm with neighbour about LEM cell count on each proc face
    MPI_Sendrecv( _sendSizes.data(), nFaces, MPI_INT, neighProc, nFaces,
                  _recvSizes.data(), nFaces, MPI_INT, neighProc, nFaces,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE );

    // -- find how many LEM wafers in total for comm
    int sendLEMcount = std::accumulate(_sendSizes.begin(),_sendSizes.end(), 0);
    int recvLEMcount = std::accumulate(_recvSizes.begin(),_recvSizes.end(), 0);

    int nScalar  = 5 + nSpecies; // T, p, dx, rho, E Y1 Y2 Y3...

    std::vector<double> SEND(nScalar * sendLEMcount,-1);
    std::vector<double> RECV(nScalar * recvLEMcount,-1);

    // -- package for sending
    unsigned int kx = 0;
    for(int jx = 0; jx < nFaces; jx++){

        int facex = boundaryFaceSC[patchx][jx];
        std::list<CellData> &cellList = spliceList[facex];
        std::list<CellData>::iterator it = cellList.begin();

        for(;it!= cellList.end();it++){
            SEND[kx + 0] = (*it).T;
            SEND[kx + 1] = (*it).dx;
            SEND[kx + 2] = (*it).p;
            SEND[kx + 3] = (*it).rho;
            SEND[kx + 4] = (*it).E;
            std::memcpy(SEND.data() + kx + 5, (*it).Y, nSpecies * sizeof(double));
            kx += nScalar;
        }
        cellList.clear();
    }

    // MPI send recv tag = nFaces
    // -- Comm with neighbour about LEM cell data for patch
    MPI_Sendrecv( SEND.data(), SEND.size(), MPI_DOUBLE, neighProc, nFaces,
                  RECV.data(), RECV.size(), MPI_DOUBLE, neighProc, nFaces,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE );

    // -- unpack RECV and fill up the cells
    kx = 0;
    for(int jx = 0; jx < nFaces; jx++){

        int facex = boundaryFaceSC[patchx][jx];
        std::list<CellData> &cellList = spliceList[facex];
        cellList.clear(); // just for safety
        CellData CELL(nSpecies, -1); // some arbitrary size

        for(int nx = 0; nx < _recvSizes[jx]; nx++){
            CELL.T   = RECV[kx + 0];
            CELL.dx  = RECV[kx + 1];
            CELL.p   = RECV[kx + 2];
            CELL.rho = RECV[kx + 3];
            CELL.E   = RECV[kx + 4];
            std::memcpy(CELL.Y, RECV.data() + kx + 5,  nSpecies * sizeof(double));
            
            CELL.m   = CELL.rho * CELL.dx;

            // push into cellList
            cellList.push_back(CELL);
            kx += nScalar;
        }
    }

    // send and recieve the proc boundary Gammas: 16th Nov 2020
    std::vector<double> GammaSEND(nFaces, 0);
    std::vector<double> GammaRECV(nFaces, 0);

    // package
    for(int jx = 0; jx < nFaces; jx++){
        int facex = boundaryFaceSC[patchx][jx];
        GammaSEND[jx] = GammaList[facex];
    }
    
    MPI_Sendrecv( GammaSEND.data(), GammaSEND.size(), MPI_DOUBLE, neighProc, nFaces,
                  GammaRECV.data(), GammaRECV.size(), MPI_DOUBLE, neighProc, nFaces,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE );


    // un-package
    for(int jx = 0; jx < nFaces; jx++){
        int facex = boundaryFaceSC[patchx][jx];
        GammaList[facex] = GammaRECV[jx];
    }

}

Info << "Parallel splicing done" << endl;
