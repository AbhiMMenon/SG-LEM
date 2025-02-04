/* Initializize Cmean using species and Zmean */

if(initZC_LES){
dimensionedScalar Tmax = Foam::max(T);
dimensionedScalar Tmin = Foam::min(T);

Cmean = (T-Tmin)/(Tmax- Tmin);

forAll(Cmean, cellx){
    std::vector<double> y2(oj.ns,-1);

    forAll(Y,i){
        volScalarField& Yi = Y[i];
        y2[i] = Yi[cellx];
    }

    Zmean[cellx]= lemLines[0].getBilgerMixFrac(y2.data());
  //std::vector<double> CC11 = lemLines[0].getProgVble_nonpremixSp(Zmean[cellx],y2.data(),T[cellx]) ;
  //Cmean[cellx]= CC11[0];
  //cSource[cellx]= CC11[1];

}
}
