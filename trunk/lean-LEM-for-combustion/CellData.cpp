#include "CellData.h"
#include <iostream>

//using namespace std;

int CellData::ns = 0.0;


/**
 * @brief       constructor for Cell object
 * @param       \input  _nsp    number of species
 * @param       \input  _dx     cell width
 *
 * @author      M. Oevermann
 */
CellData::CellData
    (
        const int       _nsp,
        const double    _dx
    )
{
    dx = _dx;
    Y  = new double[_nsp];
    DD = new double[_nsp];

    // initialize Z and C A.Menon
    Z = 0;
    C = 0;
    dCdt = 0;

}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
CellData::CellData(
    const int _nsp,
    const double _dx,
    const double _T,
    const double _p,
    const double _rho,
    const double* _Y) : dx(_dx), T(_T), p(_p), rho(_rho)
{
    Y   = new double[_nsp];
    DD = new double[_nsp];
    for (int s = 0; s < _nsp; s++)
      Y[s] = _Y[s];

    m   = _rho * _dx;
    ns  = _nsp;

    C = 0;

}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
CellData::CellData(const CellData& c)
{
    T   = c.T;
    rho = c.rho;
    p   = c.p;
    m   = c.m;
    dx  = c.dx;
    Y   = new double[ns];
    DD  = new double[ns];

    for (int s = 0; s < ns; s++)
      Y[s] = c.Y[s];


    E   = c.E;
    // initialize Z and C A.Menon
    Z = c.Z;
    C = c.C;
    dCdt = c.dCdt;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void CellData::operator=(const CellData& c)
{

    T   = c.T;
    rho = c.rho;
    p   = c.p;
    m   = c.m;
    dx  = c.dx;

    delete[] Y;
    Y   = new double[ns];
    DD  = new double[ns];


    for (int s = 0; s < ns; s++)
      Y[s] = c.Y[s];



    E   = c.E;

    // initialize Z and C A.Menon
    Z = c.Z;
    C = c.C;
    dCdt = c.dCdt;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
CellData::~CellData()
{
  delete[] Y;
  delete[] DD;
}


