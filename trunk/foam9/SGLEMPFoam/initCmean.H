
if(true){
    volScalarField& Yi = Y[indO2];
    forAll(Cmean, cellx){
        Cmean[cellx] = lemLines[0].getProgVble_premixO2(Yi[cellx]);
    }
}
