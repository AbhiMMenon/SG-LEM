/* Identify all the super-cells that fall into this arbitrary plane, the
 * conditioned  values can then be output in timing.H
 *
 *
 * This is specifically hard-coded for the JChen case as a trial
 */


// dummy line for conditioning on the plane
std::vector<std::vector<int>> dumpLines(5);
LEMLINE dummy = LEMLINE(1, 1, 1.0, 1, 1.0, file, oj, 0, cBinSize);
dummy.initCMemory();
dummy.setPremixedData(eqRatio);
dummy.cells.clear();

// -- make directory for planes
if(writePlanes){

    if(Pstream::myProcNo()== 0){
    Info << "MAKING condPLANES" << Foam::endl;
    makeDir = "rm -rf condPlanes";
    system(makeDir.c_str());

    makeDir = "mkdir  condPlanes";
    system(makeDir.c_str());

    makeDir = "mkdir  condPlanes/plane0";
    system(makeDir.c_str());

    makeDir = "mkdir  condPlanes/plane1";
    system(makeDir.c_str());

    makeDir = "mkdir  condPlanes/plane2";
    system(makeDir.c_str());

    makeDir = "mkdir  condPlanes/plane3";
    system(makeDir.c_str());

    makeDir = "mkdir  condPlanes/plane4";
    system(makeDir.c_str());
    }


    //-- locations of planes is hard-coded for now, look to have a dictionary
    //-- instead for future cases.

    double H = 1.466e-2;
    Foam::point  cond0(0.5*H, 0 ,0.005);
    Foam::point  cond1(0.75*H, 0 ,0.005);
    Foam::point  cond2(1.00*H, 0 ,0.005);
    Foam::point  cond3(1.25*H, 0 ,0.005);
    Foam::point  cond4(1.5*H, 0 ,0.005);

    Foam::vector dir1(1,0,0);

    Foam::plane pl0(cond0, dir1);
    Foam::plane pl1(cond1, dir1);
    Foam::plane pl2(cond2, dir1);
    Foam::plane pl3(cond3, dir1);
    Foam::plane pl4(cond4, dir1);

    Foam::cuttingPlane cutPl0(pl0,mesh,true,labelUList::null());
    Foam::cuttingPlane cutPl1(pl1,mesh,true,labelUList::null());
    Foam::cuttingPlane cutPl2(pl2,mesh,true,labelUList::null());
    Foam::cuttingPlane cutPl3(pl3,mesh,true,labelUList::null());
    Foam::cuttingPlane cutPl4(pl4,mesh,true,labelUList::null());


    labelList cutList0 = cutPl0.cutCells();
    labelList cutList1 = cutPl1.cutCells();
    labelList cutList2 = cutPl2.cutCells();
    labelList cutList3 = cutPl3.cutCells();
    labelList cutList4 = cutPl4.cutCells();


    // -- identify which LEM lines correspond to which plane, for THIS
    // -- processor.

    for(size_t ix = 0; ix < nCells; ix++){
        for(int cellx:cellListSC[ix]){
            bool exists = std::find(std::begin(cutList0), std::end(cutList0), cellx)
                != std::end(cutList0);
            if (exists){dumpLines[0].push_back(ix);continue;}

            exists = std::find(std::begin(cutList1), std::end(cutList1), cellx)
                != std::end(cutList1);
            if (exists){dumpLines[1].push_back(ix);continue;}

            exists = std::find(std::begin(cutList2), std::end(cutList2), cellx)
                != std::end(cutList2);
            if (exists){dumpLines[2].push_back(ix);continue;}

            exists = std::find(std::begin(cutList3), std::end(cutList3), cellx)
                != std::end(cutList3);
            if (exists){dumpLines[3].push_back(ix);continue;}

            exists = std::find(std::begin(cutList4), std::end(cutList4), cellx)
                != std::end(cutList4);
            if (exists){dumpLines[4].push_back(ix);continue;}
        }
    }
}

//Pout << "Dump lines :";
//for(int dL:dumpLines) Pout << dL <<" ";
//Pout << Foam::endl;

