/* Initialize the scalarC matrix with either zeros or with the laminar
 * flamelet, only if the beta PDFs are being used.
 *
 */

//if(!LESLEMflag){
if(true){
    Info << "Initializing matrices..\n";
    lemLines[0].setEquilTY_Z();

    if(initMatrices == "simple"){
        lemLines[0].initScalarZC_linear(); // initialized to simple
    }
    else if (initMatrices == "zero"){
        lemLines[0].initScalarZC_zero(); // initialized to zeros
    }

    else if(initMatrices=="zeroD" || initMatrices == "read" ){
        size_t RSIZE = lemLines[0].scalarZC[0][0].size();

        std::vector<double> BUFF(zBinSize*cBinSize*RSIZE,-100);
        std::vector<int> BUFF_found(zBinSize*cBinSize,0);
        if(Pstream::myProcNo()==0){
          //if(initMatrices == "zeroD"){
          //    lemLines[0].initScalarZC_0D(flamLimitRich);
          //    file.open("init1.dat");
          //    lemLines[0].writeZCProps();
          //    file.close();
          //}
            if(initMatrices == "read"){
                Info << "Reading init files";
                lemLines[0].initScalarZC_read();
            }
        }

        if(procPatch.size()>0){
            // package for sending
            for(size_t ix = 0; ix <zBinSize; ix++){
                for(size_t jx = 0; jx <cBinSize; jx++){
                    for(size_t kx = 0; kx < RSIZE; kx++){
                        BUFF[ix*(cBinSize*RSIZE) + jx*RSIZE + kx] = lemLines[0].scalarZC[ix][jx][kx];
                    }
                    BUFF_found[ix*(cBinSize) + jx ] = (lemLines[0].scalarZC_flag[ix][jx]==true)? 1: 0;
                }
            }
            Info << "Broadcasting to other procs\n";
            MPI_Bcast(BUFF.data(), zBinSize*cBinSize*RSIZE, MPI_DOUBLE,0, MPI_COMM_WORLD);
            MPI_Bcast(BUFF_found.data(), zBinSize*cBinSize, MPI_INT,0,MPI_COMM_WORLD);


            for(size_t ix = 0; ix <zBinSize; ix++){
                for(size_t jx = 0; jx <cBinSize; jx++){
                    for(size_t kx = 0; kx < RSIZE; kx++){
                      lemLines[0].scalarZC[ix][jx][kx]= BUFF[ix*(cBinSize*RSIZE) + jx*RSIZE + kx];
                    }
                    lemLines[0].scalarZC_flag[ix][jx]= BUFF_found[ix*(cBinSize) + jx ] == 1? true:false;
                    lemLines[0].dCdt_flag[ix][jx]= lemLines[0].scalarZC_flag[ix][jx];
                }
            }

        }
    }
    else{
        Info << "Unkown option '" << initMatrices << "'..exiting\n";
        exit(0);
    }

    // -- erase C production if needed
    if(!initdCdT){
        cSource = 0.0*cSource;
        for(size_t ix = 0; ix < zBinSize; ix++)
            for(size_t jx = 0; jx < cBinSize; jx++){
                lemLines[0].scalarZC[ix][jx][ITheta] = 0.0;
                lemLines[0].dCdt_flag[ix][jx] = false;
            }


    }

    // write this initial matrix
    if(Pstream::myProcNo()==2){
        file.open("init.dat");
        lemLines[0].writeZCProps();
        file.close();
    }


    Info << "Updating matrices..\n";
    // copy all the matrices from lemLines[0]
    for(size_t ix = 1; ix < nCells; ix++){
        lemLines[ix].copyZCMats(lemLines[0]);
     // lemLines[ix].initScalarZ();
    }
}

