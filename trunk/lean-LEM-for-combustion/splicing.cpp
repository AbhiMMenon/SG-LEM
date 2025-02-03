#include "LEM.h"
#include "LEMLine.h"
#include "CellData.h"


/*-----------------------------------------------------------------------------
 *				SPLICING
 *
 * Subroutines for splicing to and from the LEM line,
 *
 * @authors  - M.Oevermann, S. Arshad
 * --------------------------------------------------------------------------*/




/*-----------------------------------------------------------------------------
   LEMLINE::splice() removes elements with mass m from the right end of
   the current LEM line and adds them to the left side of toLine

-----------------------------------------------------------------------------*/
void LEMLINE::splice(double m,
                     LEMLINE &toLine)
{
  std::list<CellData>::iterator it;
  int count = cells.size();


  double perc_mass = 0.0;
  double mass_n= 0.0;
  mass_n=massOnLine();
  perc_mass= m/mass_n;
  perc_mass *= 100;


  it = cells.end(); --it;

  double splicedMass   = 0.0;
  double splicedLength = 0.0;
  double splicedE      = 0.0;

  double nSplice = 0;

  while (splicedMass < m)
  {
    splicedMass += (*it).m;
    splicedLength += (*it).dx;
    splicedE += (*it).E;
    ++nSplice;
    --it;
  }
  ++it;


  // Add a cell to compensate for residual mass
  double dm = splicedMass - m;
  if (dm / m > 1.e-5)
  {
    double alpha = dm / (*it).m;
    cells.insert(it, (*it));
    (*it).dx = (1. - alpha) * (*it).dx;
    (*it).m  = (1. - alpha) * (*it).m;
    (*it).E  = (1. - alpha) * (*it).E;
    --it;

    splicedLength -= alpha * (*it).dx;
    splicedMass -= alpha * (*it).m;
    splicedE -= alpha * (*it).E;

    (*it).dx = alpha * (*it).dx;
    (*it).m  = alpha * (*it).m;
    (*it).E  = alpha * (*it).E;
    ++it;
    ++nlem;
  }

  length -= splicedLength;
  nlem   -= nSplice;
  Energy -= splicedE;

  toLine.length += splicedLength;
  toLine.nlem   += nSplice;
  toLine.Energy += splicedE;

  toLine.cells.splice(toLine.cells.begin(), cells, it, cells.end());

}

/*------------------------------------------------------------------------------
   LEMLINE::splice() removes elements with mass m from the right end of
   the current LEM line and adds them to the the global list
------------------------------------------------------------------------------*/
void LEMLINE::spliceToList(double m,
    LEMLINE &toLine, std::list<CellData>& SplicedMassList)
{

  std::list<CellData>::iterator it;

  double perc_mass = 0.0;
  double mass_n= 0.0;
  mass_n=massOnLine();
  perc_mass= m/mass_n;
  perc_mass *= 100;


  it = cells.end(); --it;

  double splicedMass   = 0.0;
  double splicedLength = 0.0;
  double splicedE      = 0.0;

  double nSplice = 0;

  while (splicedMass < m)
  {
    splicedMass += (*it).m;
    splicedLength += (*it).dx;
    splicedE += (*it).E;
    ++nSplice;
    --it;
  }
  ++it;

  double dm = splicedMass - m;
  if (dm / m > 1.e-5)
  {
    double alpha = dm / (*it).m;
    cells.insert(it, (*it));
    (*it).dx = (1. - alpha) * (*it).dx;
    (*it).m  = (1. - alpha) * (*it).m;
    (*it).E  = (1. - alpha) * (*it).E;
    --it;

    splicedLength -= alpha * (*it).dx;
    splicedMass -= alpha * (*it).m;
    splicedE -= alpha * (*it).E;

    (*it).dx = alpha * (*it).dx;
    (*it).m  = alpha * (*it).m;
    (*it).E  = alpha * (*it).E;

    ++it;
    ++nlem;
  }

  length -= splicedLength;
  nlem   -= nSplice;
  Energy -= splicedE;

  toLine.length += splicedLength;
  toLine.nlem   += nSplice;
  toLine.Energy += splicedE;

  SplicedMassList.splice(SplicedMassList.begin(), cells, it, cells.end());

}

/*------------------------------------------------------------------------------
   LEMLINE::splice() removes elements with mass m from the right end of
   the current LEM line and adds them to the the global list

   Simpler version with no toLine

   @author A. Menon, 7th oct 2020
------------------------------------------------------------------------------*/
void LEMLINE::spliceToList(double m, std::list<CellData>& SplicedMassList)
{

  std::list<CellData>::iterator it;

  double perc_mass = 0.0;
  double mass_n= 0.0;
  mass_n=massOnLine();
  perc_mass= m/mass_n;
  perc_mass *= 100;


  it = cells.end(); --it;

  double splicedMass   = 0.0;
  double splicedLength = 0.0;
  double splicedE      = 0.0;

  double nSplice = 0;

  while (splicedMass < m)
  {
    splicedMass += (*it).m;
    splicedLength += (*it).dx;
    splicedE += (*it).E;
    ++nSplice;
    --it;
  }
  ++it;

  double dm = splicedMass - m;
  if (dm / m > 1.e-5)
  {
    double alpha = dm / (*it).m;
    cells.insert(it, (*it));
    (*it).dx = (1. - alpha) * (*it).dx;
    (*it).m  = (1. - alpha) * (*it).m;
    (*it).E  = (1. - alpha) * (*it).E;
    --it;

    splicedLength -= alpha * (*it).dx;
    splicedMass   -= alpha * (*it).m;
    splicedE      -= alpha * (*it).E;

    (*it).dx = alpha * (*it).dx;
    (*it).m  = alpha * (*it).m;
    (*it).E  = alpha * (*it).E;

    ++it;
    ++nlem;
  }

  length -= splicedLength;
  nlem   -= nSplice;
  Energy -= splicedE;
  SplicedMassList.splice(SplicedMassList.begin(), cells, it, cells.end());

}

/*------------------------------------------------------------------------------
   LEMLINE::splice() removes elements with length L from the right end of the
   current LEM line and adds them to the the global list

   Simpler version with no toLine

   @author A. Menon, 30th nov 2020
------------------------------------------------------------------------------*/
void LEMLINE::spliceToListLB(double Lfrac, std::list<CellData>& SplicedMassList)
{

    if(Lfrac <= 0.){
        std::cout << "\nSplice Length error .. skipping";
        std::cout << "\n L = " << Lfrac<< "; lLen = "<< lengthOfLine()<<  "\n";
        return;
    }

    Lfrac = Lfrac > lengthOfLine()/2? lengthOfLine()/2: Lfrac;

    std::list<CellData>::iterator it;


    double L = Lfrac;

    it = cells.end(); --it;
    double splicedMass   = 0.0;
    double splicedLength = 0.0;
    double splicedE      = 0.0;

    int nSplice = 0;

    while (splicedLength < L)
    {
        splicedMass   += (*it).m;
        splicedLength += (*it).dx;
        splicedE      += (*it).E;
        ++nSplice;
        --it;
    }
    ++it;

    double dl = splicedLength - L;
    if (dl / L > 1.e-5)
    {
        double alpha = dl / (*it).dx;
        cells.insert(it, (*it));
        (*it).dx = (1. - alpha) * (*it).dx;
        (*it).m  = (1. - alpha) * (*it).m;
        (*it).E  = (1. - alpha) * (*it).E;
        --it;

        splicedLength -= alpha * (*it).dx;
        splicedMass   -= alpha * (*it).m;
        splicedE      -= alpha * (*it).E;

        (*it).dx = alpha * (*it).dx;
        (*it).m  = alpha * (*it).m;
        (*it).E  = alpha * (*it).E;

        ++it;
        ++nlem;
    }

    length -= splicedLength;
    nlem   -= nSplice;
    Energy -= splicedE;
    SplicedMassList.splice(SplicedMassList.begin(), cells, it, cells.end());

}

/*------------------------------------------------------------------------------
   LEMLINE::splice() removes elements with length L from the left end of the
   current LEM line and adds them to the the global list


   @author A. Menon, 8th Feb 2024
------------------------------------------------------------------------------*/
void LEMLINE::spliceToListLBR(double Lfrac,double Cface, std::list<CellData>& SplicedMassList)
{

    if(Lfrac >= lengthOfLine() || Lfrac <= 0.){
        std::cout << "\nSplice Length error .. skipping";
        std::cout << "\n L = " << Lfrac<< "; lLen = "<< lengthOfLine()<<  "\n";
        return;
    }


    std::list<CellData>::iterator it ;
    std::list<CellData>::iterator itEnd = cells.end();
    double L = Lfrac;
    double CmeanLEM = getCMean();

    //-- find starting point of splicing based on Cface
    if(CmeanLEM > 1e-5){
     // int steps = (1.-Cface)*cells.size()*(1. - 2.*L/lengthOfLine());
        int steps = (1.-Cface)*cells.size()*(0.75);
        std::advance(itEnd,-1*steps);
    }

    it = itEnd; --it;
    double splicedMass   = 0.0;
    double splicedLength = 0.0;
    double splicedE      = 0.0;

    int nSplice = 0;

    while (splicedLength < L)
    {
        splicedMass   += (*it).m;
        splicedLength += (*it).dx;
        splicedE      += (*it).E;
        ++nSplice;
        --it;
    }
    ++it;

    double dl = splicedLength - L;
    if (dl / L > 1.e-5)
    {
        double alpha = dl / (*it).dx;
        cells.insert(it, (*it));
        (*it).dx = (1. - alpha) * (*it).dx;
        (*it).m  = (1. - alpha) * (*it).m;
        (*it).E  = (1. - alpha) * (*it).E;
        --it;

        splicedLength -= alpha * (*it).dx;
        splicedMass   -= alpha * (*it).m;
        splicedE      -= alpha * (*it).E;

        (*it).dx = alpha * (*it).dx;
        (*it).m  = alpha * (*it).m;
        (*it).E  = alpha * (*it).E;

        ++it;
        ++nlem;
    }

    length -= splicedLength;
    nlem   -= nSplice;
    Energy -= splicedE;
    SplicedMassList.splice(SplicedMassList.begin(), cells, it, itEnd);

    if(cells.size() != nlem)
    std::cout << "\n nLEM before " << nlem<< "; after splice = "<< cells.size()<<  "\n";
}


/*------------------------------------------------------------------------------
   LEMLINE::spliceFromList() removes elements from the global list and put in
   the
   left side of current LEM line
------------------------------------------------------------------------------*/
void LEMLINE::spliceFromList(std::list<CellData>& InSplicedMassList)
{
   std::list<CellData>::iterator it;

   // update line stats, A. Menon Oct 7th 2020
   it =  InSplicedMassList.begin();
   for(;it!=InSplicedMassList.end();it++){

       length +=  (*it).dx;
       Energy +=  (*it).E;
   }

   if(wallPresent){
       it =  cells.end();
       --it;
       cells.splice(it, InSplicedMassList);
   }
   else{

       cells.splice(cells.begin(), InSplicedMassList);
   }
   nlem = cells.size();
   if(nlem == 0){
       std::cout <<"\n\n WARNING NO LEM CELLS";
   }
}

/*------------------------------------------------------------------------------
   LEMLINE::spliceFromList() removes elements from the global list and put in
   the left side of current LEM line after contracting/expanding the LEM cells.

   16th Nov 2020
------------------------------------------------------------------------------*/
void LEMLINE::spliceFromList(std::list<CellData>& InSplicedMassList, double
        expRatio)
{
   std::list<CellData>::iterator it;

   // update line stats, A. Menon Oct 7th 2020
   it =  InSplicedMassList.begin();
   for(;it!=InSplicedMassList.end();it++){
       (*it).dx *= expRatio;
       (*it).m *= expRatio; // instead of  having a separate variable for cross-section, simply modify the mass
       length +=  (*it).dx;
       Energy +=  (*it).E;
   }

   if(wallPresent){
       it =  cells.end();
       --it;
       cells.splice(it, InSplicedMassList);
   }
   else
       cells.splice(cells.begin(), InSplicedMassList);
   nlem = cells.size();

   if(nlem == 0){
       std::cout <<"\n\n WARNING NO LEM CELLS";
   }
}


/*------------------------------------------------------------------------------
   LEMLINE::EraseFromOutlet() removes elements with mass m from the right end of
   the current LEM line and adds them to the the global list
------------------------------------------------------------------------------*/
void LEMLINE::EraseFromOutlet(double m)
{
  std::list<CellData>::iterator it;

  double perc_mass = 0.0;
  double mass_n= 0.0;
  mass_n=massOnLine();
  perc_mass= m/mass_n;
  perc_mass *= 100;


  it = cells.end(); --it;

  double splicedMass   = 0.0;
  double splicedLength = 0.0;
  double splicedE      = 0.0;

  double nSplice = 0;

  while (splicedMass < m)
  {
    splicedMass   += (*it).m;
    splicedLength += (*it).dx;
    splicedE      += (*it).E;
    ++nSplice;
    --it;
  }
  ++it;

  double dm = splicedMass - m;
  if (dm / m > 1.e-5)
  {
    double alpha = dm / (*it).m;
    cells.insert(it,  (*it));
    (*it).dx = (1. - alpha) * (*it).dx;
    (*it).m  = (1. - alpha) * (*it).m;
    (*it).E  = (1. - alpha) * (*it).E;
    --it;

    splicedLength -= alpha * (*it).dx;
    splicedMass   -= alpha * (*it).m;
    splicedE      -= alpha * (*it).E;

    (*it).dx = alpha * (*it).dx;
    (*it).m  = alpha * (*it).m;
    (*it).E  = alpha * (*it).E;

    ++it;
    ++nlem;
  }

  length -= splicedLength;
  nlem   -= nSplice;
  Energy -= splicedE;
  cells.erase(it, cells.end());

}

/*------------------------------------------------------------------------------
   LEMLINE::splice() removes elements with mass m from the right end of
   the current LEM line and adds them to the the global list
------------------------------------------------------------------------------*/
void LEMLINE::AddOneLEMCell(double m_in, double T_in, double p_in,
std::vector<double> Y_in)
{
  auto gas = lem.gas;
  std::list<CellData>::iterator it;
  it = cells.begin();

  cells.insert(it, (*it));
  --it;

  (*it).m  = m_in;
  (*it).T  = T_in;
  (*it).p  = p_in;
  for (int s = 0; s < lem.ns; s++)
  {
    (*it).Y[s] = Y_in[s];
  }

  gas->setState_TPY((*it).T, (*it).p, (*it).Y);
  (*it).rho = gas->density();
  (*it).dx = (*it).m/(*it).rho;
  double E_cell = (gas->enthalpy_mass() - (*it).p / (*it).rho) * (*it).m;
  (*it).E = E_cell;

  Energy += E_cell;
  length += (*it).dx;
  nlem++;


}

/*-----------------------------------------------------------------------------
 * @brief   new re-grid code, preserves length of line, perserves gradients.
 * More suitable for SG-LEM to preserve flame structures.
 *
 * @param \input LEMres
 *        \input tolFac
 *
 * tolFac is the trigger for cell merging, which is mass based, merging should
 * only be used for small spliced-in fragments.
 *
 * @author  A.Menon
 * --------------------------------------------------------------------------*/

void LEMLINE::regrid( double LEMres, double tolFac)
{
  auto gas = lem.gas;
  std::list<CellData>::iterator it, it2, it3;
  double *Y  = lem.tmp_s1;

  //LEMres = lengthOfLine()/LEMres;

  // Two step procedure
  // 1 - merge: cell too small? merge with next cell, delete next cell
  // 2 - split: cell too big? Add cells before cell  of correct resolution,
  // linear interpolate species and temperature to capture
  // gradients, adjust cell size and mass

  it = cells.begin();
  it2 = it;
  ++it2;
  for(;it2!=cells.end()&&it!=cells.end();){
      double DX = ((*it).dx + (*it2).dx);

      if(LEMres > DX){
        double M =  (*it).m + (*it2).m;
        double MR1 =  (*it).m/M;
        double MR2 =  (*it2).m/M;
        (*it).m = M;

        for (int s = 0; s < lem.ns; s++)
          Y[s] = MR1*(*it).Y[s] + MR2*(*it2).Y[s];

        double T = MR1*(*it).T + MR2*(*it2).T;
        double P = 0.5*(*it).p + 0.5*(*it2).p; // should not matter
        gas->setState_TPY(T,P,Y);
        gas->getMassFractions((*it).Y);

        (*it).T = T;
        (*it).p = P;
        (*it).rho = gas->density();
        (*it).dx = (*it).m/gas->density();
    //  (*it).E = (*it).m*(gas->enthalpy_mass() - (*it).p/(*it).rho);
        it2 = cells.erase(it2);
      }
      else{++it2; ++it;
      }
  }

//return;
 it2 = cells.begin();
 it = it2;

 for(;it2!=cells.end()&&it!=cells.end();++it2){

      double Lnow = (*it2).dx;
      if(Lnow > LEMres){
          int needCell = std::floor(Lnow/LEMres);
          double comp = (Lnow - LEMres*needCell);
          // -- insert cells
          for(int ix = 0; ix < needCell-1; ix++){
             it2=  cells.insert(it2, (*it2));
             ++it2;
          }

          (*it2).m=(LEMres+comp)* (*it2).rho;
          (*it2).dx=LEMres+comp;
          // -- adjust dx and m
          it = it2; --it;
          for(int ix = 0; ix < needCell-1; ix++){
              (*it).dx = LEMres;
              (*it).m  = (*it).rho*LEMres;
              --it;
          }
      }
}

  setEnergy();

}






/*-----------------------------------------------------------------------------
 *			Making EQUIDISTANT grid after SPLICING
 *
 *			S. Arshad and Esteban
 * --------------------------------------------------------------------------*/
void LEMLINE::makeEquidistantGridafterSplicing(int noLEMrequired,
        double minCellSize)
{
  //const double minCellSize = 1.e-5;
  auto gas = lem.gas;
  std::list<CellData>::iterator it, itNext, itlast;

  double L = 0.0;
  for(it = cells.begin(); it != cells.end(); ++it)
      L += (*it).dx;
  double nolemrequired=noLEMrequired;
  // constant grid spacing we want to achieve
  const double Dx = L / nolemrequired;
  double x, y, yNext;           // x = distance on equidistant grid
                                // y = distance on current LEM line
  int nInsertedCells = 0;

  it = itNext = cells.begin();
  ++itNext;
  y = (*it).dx;
  yNext = y + (*itNext).dx;
  x = Dx;

  do
  {
//     cout << "the number of current lem cell is" <<
//     std::distance(cells.begin(), (*it) <<  endl;
//     cout << "x is" << x << "and y is" << y << endl;
    if (x - y > minCellSize)
    {
      // cout << "the current lem cell" << std::distance(cells.begin(), (*it) <<
      // "is smaller than equidistant cell" << endl;
      // we might need to merge several cells
      // when we exit this while loop, we should still have x > y
      // but then we need to split the next cell exactly at position x
      while (   (x - yNext > - minCellSize)
             && (it != cells.end())
             && (itNext != cells.end()) )
      {
// 	cout << "merging cells"<< endl;
        //gas->setState_TPY((*(*it).T, (*(*it).p,(*(*(*(*it).Y);
        //double u1 = gas->enthalpy_mass();
        //double u1 = (*(*it).E / (*(*it).m + (*(*it).p / (*(*it).rho;
        //gas->setState_TPY((*itNext).T, (*itNext).p,(*itNext).Y);
        //double u2 = (gas->enthalpy_mass());// - (*itNext).p / (*itNext).rho);
        //double u2 = (*itNext).E / (*itNext).m + (*itNext).p / (*itNext).rho;
        //
        // merging cells it and itNext
        double L = (*it).dx + (*itNext).dx;
        double m = (*it).m +  (*itNext).m;
        double a = (*it).m / m;
        double oMa = 1. - a;
        double rho = m / L;
        double p   = 0.5 * ((*it).p + (*itNext).p);
        double *Y  = lem.tmp_s1;
        //double u = a * u1 + oMa * u2;
        for (int s = 0; s < lem.ns; s++)
          Y[s] = a *   (*it).Y[s] + oMa * (*itNext).Y[s];
        double T = a * (*it).T + oMa * (*itNext).T;

        //gas->setState_TPY(T,p,Y);
        //gas->setMassFractions(Y);
        //gas->setState_HP(u, p);
        //gas->setState_UV(u, 1. / rho);
        (*it).T = T; //gas->temperature();
        (*it).p = p;
        (*it).rho = rho;
        for (int s = 0; s < lem.ns; s++)
        (*it).Y[s] = Y[s];
        (*it).dx = L;
        (*it).m = m;
        gas->setState_TPY((*it).T, (*it).p, (*it).Y); // addition by Esteban
        (*it).E = (*it).m * (gas->enthalpy_mass() - (*it).p / (*it).rho);

        itNext = cells.erase(itNext);
        --nInsertedCells;

        y = yNext;
        yNext = y + (*itNext).dx;

      }
      if (    (x - y > minCellSize)
           && (it != cells.end())
           && (itNext != cells.end())  )
      {
// 	cout << "splitting cell" << endl;
        // splitting cell itNext: left part will be merged with it, right
        // part will be introduced as a new cell in list
        //gas->setState_TPY((*(*it).T, (*(*it).p,(*(*(*(*it).Y);
        // double u1 = (gas->enthalpy_mass());// - (*(*it).p / (*(*it).rho);
        //double u1 = (*(*it).E / (*(*it).m + (*(*it).p / (*(*it).rho;
        //gas->setState_TPY((*itNext).T, (*itNext).p,(*itNext).Y);
        //double u2 = (gas->enthalpy_mass());// - (*itNext).p / (*itNext).rho);
        //double u2 = (*itNext).E / (*itNext).m + (*itNext).p / (*itNext).rho;
        double lR = yNext - x;
        double lL = x - y;
        double alpha = lL/(*itNext).dx;

        double L = Dx;
        double m = (*it).m + alpha * (*itNext).m;
        double a = (*it).m / m;
        double oMa = 1. - a;
        double rho = m/L;
        double p   = 0.5 * ((*it).p + (*itNext).p);
        double *Y  = lem.tmp_s1;
        double T = a * (*it).T + oMa * (*itNext).T;

        //double u = a * u1 + oMa * u2;
        for (int s = 0; s < lem.ns; s++)
          Y[s] = a * (*it).Y[s] + oMa * (*itNext).Y[s];
        //gas->setMassFractions(Y);
        //gas->setState_HP(u, p);
        //gas->setState_UV(u, 1. / rho);

        (*it).T = T; //gas->temperature();
        (*it).p = p;
        (*it).rho = rho;
        for (int s = 0; s < lem.ns; s++)
          (*it).Y[s] = Y[s];
        (*it).dx = L;
        (*it).m = m;
        gas->setState_TPY((*it).T, (*it).p, (*it).Y); // addition by Esteban
        (*it).E = (*it).m * (gas->enthalpy_mass() - (*it).p / (*it).rho);

        (*itNext).dx = lR;
        (*itNext).m = lR * (*itNext).rho;
        (*itNext).E *= (1.-alpha);

        y = x;

        ++it;
        ++itNext;

        if (it != cells.end())
          y += (*it).dx;
        if (itNext != cells.end())
          yNext += (*itNext).dx;

        x += Dx;
      }
    }
    else if (y - x > minCellSize)
    {
//             cout << "the current lem cell" << std::distance(cells.begin(),
//             (*it) << "is bigger than equidistant cell" << endl;

      while (y > x && y - x > minCellSize)
      {
        // splitting cell by inserting a new cell with length dx before
        // position of the iterator
        double dxR = (*it).dx - Dx;
// 	cout << "dxR is" << dxR << "dx is "<< (*(*it).dx << "Dx is "<< Dx << endl;

//        if (fabs(dxR) > minCellSize)
        {
// 	  cout << "inside if loop"<< endl;

          double a = Dx / (*it).dx;
// 	  cout << "a is" << a << endl;
          //cells.insert(it,
          //    CellData(lem.ns, Dx, (*(*it).T, (*(*it).p, (*(*it).rho, (*(*(*(*it).Y));
// 	  cout << "before inserting cell" << endl;
          cells.insert(it, (*it));
          --it;
          (*it).m *= a;
          (*it).E *= a;
          ++it;
          ++nInsertedCells;
          // shrinking size of original cell by dx
          (*it).dx -= Dx;
// 	  cout << "dx now becomes" << (*(*it).dx << endl;
          (*it).m = (*it).rho * (*it).dx;
// 	  cout << "mass now becomes "<< (*(*it).m << endl;
          x += Dx;
// 	  cout << "x now becomes"<< x << endl;
          (*it).E *= (1.-a);
// 	  cout << "energy now becomes" <<  (*(*it).E << endl;
        }

      }
    }
    else
    {
//       cout << "the current lem cell" << std::distance(cells.begin(), (*it) <<
//       "is neither smaller nor bigger than the equidistant cell" << endl;

      double diff = x-y;
//       cout << "difference is" << diff << endl;

      ++it;
//       cout << "it now is" << std::distance(cells.begin(), (*it) << endl;
      ++itNext;
//       cout << "itnext  now is" << std::distance(cells.begin(), itNext) <<
//       endl;
//       cout << "cells size is" << cells.size() << endl;
//       cout << "length is" << L << endl;

      if(std::distance(cells.begin(), itNext) == (nolemrequired+1))
      {
// 	cout << "dx on last cell is" << (*itNext).dx << "only" << endl;
// 	cout << "itnext  now is" << std::distance(cells.begin(), itNext) <<
// 	endl;

	if (
              (*itNext).dx < minCellSize  ||
              (abs(L-x)    < minCellSize) ||
              (abs(L-y)    < minCellSize)
           ){
	  it=cells.erase(it);
	  itNext = cells.begin();
	  break;
	}
      }
      x += Dx;
//       cout << "x now becomes" << x << "after adding dx" << Dx << endl;

      if (it != cells.end())
      {
        y += (*it).dx;
//       cout << "y now becomes" << y << endl; // on last lem cell it reaches
//       cells.end() therfore nothing is added to y
      }
      if (itNext != cells.end())
      {
      yNext += (*itNext).dx;
//       cout << "ynext now becomes" << yNext << endl; //ynext becomes higher
//       after lastcall. it adds the dx of lem cell number 0 becuase itnext
//       goes to 0
      }

    }
  }

  while (it != cells.end());

//   cout << "end of makeEquidistantGrid" << endl;
//   cout << "nlem required are" << nolemrequired << "and nlem now is"<< nlem
//   << endl;
//   cout << "before nlem assigning to nolemrequired" << endl;
  nlem =nolemrequired;

  if(noLEMcells())
  {
//   cout << "Warning: nlem is wrong in makeEquidistantGridafterSplicing"<<
//   endl;
  }
}

