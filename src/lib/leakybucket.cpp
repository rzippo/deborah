//
// Created by raffa on 22/02/2021.
//

#include <sstream>
#include <iomanip>
#include "leakybucket.h"

deborah::LeakyBucket::LeakyBucket(double sigma, double rho)
{
    this->sigma = sigma;
    this->rho = rho;
}

deborah::LeakyBucket::LeakyBucket(const deborah::ParametrizedLeakyBucket& plb, const simplex::Solution& sol)
{
    if(plb.sigma.getNumVariables() != sol.getNumVariables())
        throw;

    double s = plb.sigma.getConstant();
    for (int var_idx = 0; var_idx < plb.sigma.getNumVariables(); ++var_idx)
    {
        s += plb.sigma.getCoefficient(var_idx) * sol.getCoefficient(var_idx);
    }

    this->sigma = s;
    this->rho = plb.rho;
}

double deborah::LeakyBucket::getSigma() const
{
    return sigma;
}

double deborah::LeakyBucket::getRho() const
{
    return rho;
}

std::string deborah::LeakyBucket::PrintJson() const
{
    std::stringstream ss;
    ss
        << "{"
        << "\"sigma\":" << std::fixed << std::setprecision(6) << sigma << ","
        << "\"rho\":" << std::fixed << std::setprecision(6) << rho
        << "}";

    return ss.str();
}
