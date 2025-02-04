/* Initialize the scalarC matrix with either zeros or with the laminar
 * flamelet, only if the beta PDFs are being used.
 *
 */

if(true){
    Info << "Initializing matrices..\n";

    if(initMatrices=="unburnt"){
        lemLines[0].initScalarC();
    }
    
    else if(initMatrices=="zeroD" || initMatrices == "read" ){
        lemLines[0].initializeWithPhi(eqRatio); // just in case
        int RSIZE = lemLines[0].scalarC[0].size();

        std::vector<double> BUFF(cBinSize*RSIZE,-100);
        if(Pstream::myProcNo()==0){

         // if(initMatrices == "zeroD") 
         //     lemLines[0].initScalarC_0D();
         // NOTE: This function is no longer functional since the move to Cantera 3.0
         // It is needed for the first test-case (backward-facing step) but not for the Volvo case
         // else if(initMatrices == "read")
            
            if(initMatrices == "read")
                lemLines[0].initScalarC_read();
            file.open("init1.dat");
            lemLines[0].writeCProps();
            file.close();
        }

        if(procPatch.size()>0){
            // package for sending
            for(int ix = 0; ix <cBinSize; ix++){
                for(int jx = 0; jx < RSIZE; jx ++){
                    BUFF[jx + ix*RSIZE] = lemLines[0].scalarC[ix][jx];
                }
            }
            Info << "Broadcasting to other procs\n";
            MPI_Bcast(BUFF.data(), cBinSize*RSIZE, MPI_DOUBLE,0, MPI_COMM_WORLD);

            for(int ix = 0; ix <cBinSize; ix++){
                lemLines[0].scalarC_flag[ix]=1.0;
                for(int jx = 0; jx < RSIZE; jx ++){
                    lemLines[0].scalarC[ix][jx]=BUFF[jx + ix*RSIZE];
                }
            }
        }

    }
    else if(initMatrices=="zero"){
        lemLines[0].initializeWithPhi(eqRatio); // just in case
        lemLines[0].initScalarC();
    }

    else if(initMatrices=="linear"){
        lemLines[0].initializeWithPhi(eqRatio); // just in case
        lemLines[0].initScalarC_linear();
    }
    else{
        Info << "Unkown option '" << initMatrices << "'..exiting\n";
        exit(0);
    }

    // write this initial matrix
    if(Pstream::myProcNo()==2){
        file.open("init.dat");
        lemLines[0].writeCProps();
        file.close();
    }


    // -- erase C production if needed
    if(!initdCdT){
        for(int ix = 0; ix < cBinSize; ix++){
            lemLines[0].scalarC[ix][IC] = 0;
        }
    }


    Info << "Updating matrices..\n";
    // copy all the matrices from lemLines[0]
    for(size_t ix = 1; ix < nCells; ix++){
        lemLines[ix].copyMats(lemLines[0]);
    }

    Info << "Done! Rquilibrium Temp = " << lemLines[0].Tb_p << endl;
    
    Info << "Initializing BetaPDFs..\n";
    lemLines[0].initBetaShapes(cBinSize,1.0);
}
