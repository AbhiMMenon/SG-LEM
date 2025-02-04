/*-----------------------------------------------------------------------------
 * createSuperGrid.H
 *
 * Create superGrid using dbns library, a scalar field for visualizing the
 * coarse mesh, list of fine cells in each superCell
 *
 * @author A.Menon 10th sep 2020
 *---------------------------------------------------------------------------*/


// MGridGen agglomeration - sizes, aspect ratios, limits controlled in dbns
// library, must recompile shared files for changes to take effect

if(adjustSuperGrid){
Info<< "Reading  field volMax, adjusting superGrid sizes\n" << endl;
    volScalarField volMax
    (
        IOobject
        (
            "volMax",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    double mmHere = *std::min_element(volMax.begin(), volMax.end());
    volMax = mmHere;
    clusterMin = std::round(mmHere);
    clusterMax = std::round(mmHere);
}


Info << "-----------------MGRIDGEN Agglomeration---------------" << endl;
Info << "Sizes are " << clusterMin << " " << clusterMax<<endl;
Info << endl;

//if(LESLEMflag){
//    clusterMin = 1;
//    clusterMax = 1;
//    NLEM = NLEM_LESLEM;
//    cBinSize = 1;
//    zBinSize = 1;
//}
Foam::fineMgMeshLevel fineMesh(mesh);
Foam::coarseMgMeshLevel superGrid(fineMesh, clusterMin, clusterMax);

int totFaces = superGrid.nInternalFaces();
size_t nCells   = superGrid.nCells();

std::vector<double> cellVol (mesh.V().begin(), mesh.V().end());


// list of all fine-mesh cells in superCell, superCell volume, integral length
// scale
std::vector<std::vector<int>> cellListSC(nCells);
std::vector<double> cellVolSC(nCells, 0.0);
std::vector<double> Gamma(nCells, 0.0);// cross section
std::vector<double> domLengthSC(nCells);
std::vector<double> ltAvg(nCells);

// -- superCell is used for visualizing
forAll(superCell, i){

    int sCell    = fineMesh.child()[i];
    superCell[i] = sCell;

    cellListSC[sCell].push_back(i);
    cellVolSC [sCell] += cellVol[i];

}

// find integral length scale for each cell
for(size_t ix = 0; ix < nCells; ix++){
    domLengthSC[ix] = std::cbrt(cellVolSC[ix]);
    ltAvg[ix] = std::cbrt(cellVolSC[ix]/cellListSC[ix].size());
    // -- initial Gamma values
    Gamma[ix] = cellVolSC[ix] / domLengthSC[ix];
}


