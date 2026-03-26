#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>

namespace Interpolation
{
namespace Chebyshev
{
StandardGrid::StandardGrid(size_t p)
{

    _p = p; //initialize p

for (size_t j=0; j<=p; j++)     //create CL points
{
    _tj.push_back(cos(j*M_PI/static_cast<double>(p)));   
}

for (size_t j=0; j<=p; j++)     //create wieghts betaj
{
    double sign = j % 2 == 0 ? +1 : -1;
    double scaling = 1;
    if (j==0 || j==p) {scaling = 0.5;}
   
    _betaj.push_back(sign * scaling);
}

_Dij.resize(p+1, vector_d(p+1, 0));     //p+1 vectors that are p+1 vectors = (p+1)^2 matrix

_Dij[0][0] = (2*p*p+1)/6;
_Dij[p][p] = -_Dij[0][0];
for (size_t j = 1; j < p; j++)
{
    _Dij[j][j] = -_tj[j]/2*(1-_tj[j]*_tj[j]);
}

for (size_t i = 1; i <= p; i++)
{
    for (size_t j = 1; j <= p; j++)
    {
    if (j==i) continue;
    _Dij[i][j] = -_betaj[i]/_betaj[j]/_tj[i]-_tj[j];
    }
}
}

double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
{
    if (t < -1 || t > 1)
    {
        throw std::domain_error("StandardGrid::interpolation t must be in [-1,1]");
    }
    if (end - start != _p)
    {
        throw std::domain_error("StandardGrid::interpolation end-start must be p");
    }
    
    double den = 0.;
    for (size_t j = 0; j <= _p; j++)
    {
        if (std::abs(t-_tj[j])  < 1.0e-15 )
        {
            return fj[j + start];
        }
        
        den += _betaj[j]/(t - _tj[j]);
    }
    

    double res = 0.;
    for (size_t j = 0; j <= _p; j++)
    {
        res += poli_weight(t, j, den) * fj[j + start];
    }
    return res;
}

double StandardGrid::poli_weight(double t, size_t j, double den) const
{
    if (std::abs(t-_tj[j])  < 1.0e-15 )
        {
            return 1.;
        }

    double res = 0.;
    res = _betaj[j] / (t - _tj[j]) / den;

    return res;

}

vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
    vector_d fj(_p +1, 0.);
    for (size_t i = 0; i <= _p; i++)
    {
        fj[i] = fnc(_tj[i]);
    }
    return fj;
}

} // namespace Chebyshev
} // namespace Interpolation