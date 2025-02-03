#ifndef CELLDATA_H
#define CELLDATA_H


/*-----------------------------------------------------------------------------
 *
 * --------------------------------------------------------------------------*/
class CellData
{

  public:

    static int ns; // number of species

    CellData(
        const int _nsp,
        const double _dx);


    CellData(
        const int _nsp,
        const double _dx,
        const double _T,
        const double _p,
        const double _rho,
        const double *_Y);
    CellData(const CellData& c);
    void operator=(const CellData& c);
    ~CellData();

    double dx;  // dx of cell
    double m;   // mass in cell: m = rho * dx (redundant)
    double E;   // internal energy of the cell E = e * m

    double  T, p, rho;
    double* Y;
    double* DD;
//  std::vector<double> Y;      // Menon 16th Oct

    // -- A.Menon
    double  Z;         // MixF
    double  C;         // progress variable
    double  dCdt;     // change in progress variable | 6 Jan 2021
    double  dedt;     // heat production rate



    double lambda, cp;  // thermal conductivity, cp
    double a, b, c, d, x; // coefficients for Thomas algorithm

};

#endif
