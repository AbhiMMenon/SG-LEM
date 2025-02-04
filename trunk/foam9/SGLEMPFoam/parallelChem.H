if(advanceLine && parallelChem && procPatch.size() > 0){
MPI_Barrier(MPI_COMM_WORLD);

Info << "\n----------Parallel chemistry-------------" << endl;
auto paraChemStart = std::chrono::high_resolution_clock::now();

int totLEMCells = 0;
int totLEMCellsDom = 0;
for(size_t ix = 0; ix <nCells; ix++){
    if(!uniLine[ix])
    totLEMCells += lemLines[ix].cells.size();
}
MPI_Allreduce(&totLEMCells,&totLEMCellsDom, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // so all procs can make the DISTSizes properly

Info << "nLEM cells in domain=" << totLEMCellsDom<<"\n";

int nP = Pstream::nProcs();
int thisProc = Pstream::myProcNo();

int cellsPerProc = totLEMCellsDom/nP;
int remainCell = totLEMCellsDom - cellsPerProc*nP;

// -- DISTsizes is needed by all procs while RECVsizes only by root
std::vector<int> DISTsizes(nP, cellsPerProc);
std::vector<int> DISTsizesActual(nP, 0);
std::vector<int> RECVsizes(nP, 0);
std::vector<int> RECVsizesActual(nP, 0);

MPI_Gather(&totLEMCells, 1 , MPI_INT, RECVsizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

std::vector<int> disp(nP,0);
DISTsizes[nP-1] += remainCell;

int nScalar  = 3 + nSpecies; // T,p,rho and species
std::vector<double> SEND(totLEMCells * nScalar,-1);
std::vector<double> MasterDATA(1,0);


if(thisProc ==0) MasterDATA.resize(totLEMCellsDom*nScalar, -1);


// package all LEM cells for sending
size_t kx = 0;
for(size_t ix = 0; ix <nCells; ix++){
    if(!uniLine[ix]){
    LEMLINE &lNow = lemLines[ix];
    std::list<CellData>::iterator it = lNow.cells.begin();
    for(;it!= lNow.cells.end();it++){
        SEND[kx + 0] = (*it).T;
        SEND[kx + 1] = (*it).p;
        SEND[kx + 2] = (*it).rho;
        std::memcpy(SEND.data() + kx + 3, (*it).Y, nSpecies * sizeof(double));
        kx += nScalar;
    }
    }
}

int sum = 0;
for(size_t px = 0; px < nP;px++){
    RECVsizesActual[px] = RECVsizes[px]*nScalar;
    DISTsizesActual[px] = DISTsizes[px]*nScalar;

    disp[px] = sum;
    sum += RECVsizesActual[px];

}

Info << "MPI gather LEM data on root" <<Foam::endl;
MPI_Gatherv(SEND.data(), SEND.size() , MPI_DOUBLE, MasterDATA.data() , RECVsizesActual.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
Info << "Distributing LEM data evenly on procs" <<Foam::endl;

sum = 0;
for(size_t px = 0; px < nP;px++){
    disp[px] = sum;
    sum += DISTsizesActual[px];
}

std::vector<double> RECV(DISTsizesActual[thisProc],-1);
MPI_Scatterv(MasterDATA.data(), DISTsizesActual.data(),disp.data() , MPI_DOUBLE, RECV.data() , DISTsizesActual[thisProc], MPI_DOUBLE, 0, MPI_COMM_WORLD);

for(size_t ix = 0; ix < RECV.size();ix++){
  if(RECV[ix] == -1){std::cerr<<"COPYING ERROR, exiting"; exit(0);}
}

Info << "Distributed chemistry advancement" <<Foam::endl;
realtype tret;
realtype *y = NV_DATA_S(oj.y);

for(kx = 0; kx < DISTsizes[thisProc];kx+=nScalar){
    // check if data is valid
    double time = 0;
    y[IT]       = RECV[kx +0];
    y[IP]       = RECV[kx +1];
    y[IRHO]     = RECV[kx +2];
    for (int s = 0; s <nSpecies; ++s) {
      y[IY0+s] = RECV[kx + s + 3];;
    }
    oj.gas->setState_TPY(y[IT], y[IP], y + IY0); // check if this gas->state is relalizable or we made an error somewhere
    CVodeReInit(oj.cvode_mem,time, oj.y);
    CVode(oj.cvode_mem, time+tStep, oj.y, &tret, CV_NORMAL);

    RECV[kx +0] =y[IT];
    RECV[kx +1] =y[IP];
    RECV[kx +2] =y[IRHO];
    for (int s = 0; s <nSpecies; ++s) {
        RECV[kx + s + 3] = y[IY0+s];
    }
}

Info << "Gathering CVODE advanced LEM cell data" <<Foam::endl;
std::fill(MasterDATA.begin(), MasterDATA.end(), -1);

// the send buffer here is RECV, disp was calculated using DISTsizesActual so
// it remains the same
MPI_Gatherv(RECV.data(), RECV.size() , MPI_DOUBLE, MasterDATA.data() , DISTsizesActual.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
Info << "Redistributing and local binned data" <<Foam::endl;

sum = 0; // maybe have two separate disp things instead of recalculating
for(size_t px = 0; px < nP;px++){
    disp[px] = sum;
    sum += RECVsizesActual[px];
}

std::fill(SEND.begin(), SEND.end(), -1);
MPI_Scatterv(MasterDATA.data(), RECVsizesActual.data(),disp.data() , MPI_DOUBLE, SEND.data() , totLEMCells*nScalar, MPI_DOUBLE, 0, MPI_COMM_WORLD);

for(size_t ix = 0; ix < SEND.size();ix++){
  if(SEND[ix] == -1){std::cerr<<"COPYING ERROR, exiting"; exit(0);}
}

kx = 0;
for(size_t ix = 0; ix <nCells; ix++){
    if(!uniLine[ix]){
    LEMLINE &lNow = lemLines[ix];
    std::list<CellData>::iterator it = lNow.cells.begin();
    for(;it!= lNow.cells.end();it++){
        (*it).T   = SEND[kx + 0];
        (*it).p   = SEND[kx + 1];
        (*it).rho = SEND[kx + 2];
        std::memcpy((*it).Y,SEND.data() + kx + 3, nSpecies * sizeof(double));
        kx += nScalar;
    }
    }
}

// condtition all lines
for(size_t ix = 0; ix <nCells; ix++){
    LEMLINE &lNow = lemLines[ix];
    lNow.setEnergy();
    lNow.setHeatRelease();
    if(lNow.ZBinSize >0){
    lNow.setBilgerMixFrac();
        if(lNow.CBinSize >0){
        lNow.setProgVble_nonpremixSp();
        lNow.calcScalarZC();
        }
    }
    else if(lNow.scalarC.size()>0){
        lNow.setProgVble_premixSp();
        lNow.calcScalarC();
        lNow.calcSpeciesC();
    }
}

auto paraChemStop = std::chrono::high_resolution_clock::now();
auto paraChemTime   = std::chrono::duration_cast<std::chrono::microseconds>
(paraChemStop - paraChemStart);
double advTime_sec = paraChemTime.count()/1e6;
double maxAdvTime;

if(procPatch.size()>0)
    MPI_Reduce(&advTime_sec,&maxAdvTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
else
    maxAdvTime = advTime_sec;

Info <<"Done! time [max all procs] = " << maxAdvTime <<endl;
Info << "-----------------------------------------" << endl <<endl;

}
