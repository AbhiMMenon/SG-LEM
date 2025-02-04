Info << "Mapping species \n";
for(size_t ix = 0; ix < nCells;ix++){

    std::vector<double> SOL = lemLines[ix].getMeanS(pList);
    for(int cellx:cellListSC[ix]){
        // -- averaged solution

        // -- assign scalars
        RHO_mapped[cellx]  = SOL[IRho];
        H_mapped[cellx]    = SOL[Ih];
        T_mapped[cellx]    = SOL[IT];

        // -- update species
        forAll(Y,i){
            if (composition.solve(i)){
            volScalarField& Yi = Y[i];
            Yi[cellx] = SOL[IY + i];
            }
        }

        // -- update production rates
        forAll(Yprod,i){
            volScalarField& Ypi = Yprod[i];
            Ypi[cellx] = SOL[IY+i+oj.ns];
        }
    }
}
