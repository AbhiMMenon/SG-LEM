/*
 * This is specifically hard-coded for the JChen case as a trial
 *
 * Brief: make a large LEM line with all the CELLs that lie on each plane,
 * condition on C, write this for each plane
 */

if(wProps && writePlanes){
    size_t nPlanes = dumpLines.size();
    size_t nScalar  = 5 + nSpecies; // T, p, dx, rho, E Y1 Y2 Y3...
    Info << "Writing condtitioned planes" << Foam::endl;

    for(size_t np = 0; np < nPlanes;np++){
        //-- clear line
        dummy.cells.clear();

        size_t plainCells=0;
        size_t totPlainCells=0;
        size_t nP = Pstream::nProcs();
        std::vector<double> SEND;
        std::vector<double> RECV;
        std::vector<int> recvSizes(nP,0);
        std::vector<int> disp(nP,0);

        for(size_t dx:dumpLines[np]){
            LEMLINE copy = lemLines[dx];
            dummy.spliceFromList(copy.cells);
        }

        //-- find send size
        plainCells = dummy.cells.size();
        SEND.resize(nScalar*plainCells);

        //-- find recv size
        if(procPatch.size()>0){
            MPI_Reduce(&plainCells,&totPlainCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Gather(&plainCells, 1 , MPI_INT, recvSizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
        else{
            totPlainCells = plainCells;
        }

        Info << "Recv sizes for plane " << np << ":" ;
        for(size_t px = 0; px < nP;px++) Info<< recvSizes[px]*nScalar << " ";
        Info << Foam::endl;

        // -- package for sending
        if(plainCells>0){
            std::list<CellData>::iterator it = dummy.cells.begin();
            size_t kx = 0;
            for(;it!= dummy.cells.end();it++){
                SEND[kx + 0] = (*it).T;
                SEND[kx + 1] = (*it).dx;
                SEND[kx + 2] = (*it).p;
                SEND[kx + 3] = (*it).rho;
                SEND[kx + 4] = (*it).E;
                std::memcpy(SEND.data() + kx + 5, (*it).Y, nSpecies * sizeof(double));
                kx += nScalar;
            }
        }

        // --- allot memory for MPI gather on proc 0
        RECV.resize(nScalar*totPlainCells,-1);
        // -- small operation to modify recvSizes

        recvSizes[0] *= nScalar;
        disp[0] = 0;
        for(size_t px = 1; px < nP;px++){
            disp[px] += recvSizes[px-1]+disp[px-1];
            recvSizes[px] *= nScalar;
        }

        Info << "Disp sizes for plane " << np << ":" ;
        for(size_t px = 0; px < nP;px++) Info<< disp[px] << " ";
        Info << Foam::endl;

       ///-- gather all LEM cells that lie on this plane in proc 0
       Info << "Line data gathered for plane " << np << Foam::endl;
       MPI_Gatherv(SEND.data(), SEND.size() , MPI_DOUBLE,
                   RECV.data() , recvSizes.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // -- unpackage recv data to dummy
        if (Pstream::myProcNo()==0){
            dummy.cells.clear();
            CellData CELL(nSpecies, 1e-4); // some arbitrary size

            for(size_t nx = 0; nx < RECV.size(); nx+=nScalar){
                CELL.T   = RECV[nx + 0];
                CELL.dx  = RECV[nx + 1];
                CELL.p   = RECV[nx + 2];
                CELL.rho = RECV[nx + 3];
                CELL.E   = RECV[nx + 4];
                std::memcpy(CELL.Y, RECV.data() + nx + 5,  nSpecies * sizeof(double));
                CELL.m   = CELL.rho * CELL.dx;


                // push into cellList
                dummy.cells.push_back(CELL);
            }


            // -- condition
            dummy.clearCmats();
            dummy.setProgVble_premixO2();
            dummy.setEnergy();
            dummy.calcScalarC();
            dummy.calcSpeciesC();

            //-- write to the file
            file.open("condPlanes/plane"+std::to_string(np)+"/"+ std::to_string(fileNo)+".dat");
            file << "## nCells in dummy line for plane " + std::to_string(np) +" :" << dummy.cells.size();
            file << endl;
            dummy.writeCProps();
            file.close();
        }
        Info << "Conditioned and writen for plane " << np << Foam::endl;
    }

}
