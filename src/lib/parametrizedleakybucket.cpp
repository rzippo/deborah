#include "parametrizedleakybucket.h"
#include <iostream>

namespace deborah
{

ParametrizedLeakyBucket::ParametrizedLeakyBucket(uint64_t num_variables) : sigma(num_variables)
{
    sigma.Zero();
    rho = 0.0;
    nvars = num_variables;
}

ParametrizedLeakyBucket::ParametrizedLeakyBucket(const ParametrizedLeakyBucket &lb) : sigma(lb.nvars)
{
    sigma = lb.sigma;
    rho = lb.rho;
    nvars = sigma.getNumVariables();
}

ParametrizedLeakyBucket::~ParametrizedLeakyBucket()
{
}

ParametrizedLeakyBucket &ParametrizedLeakyBucket::operator=(const ParametrizedLeakyBucket &lb)
{
    if (&lb != this)
    {
        sigma = lb.sigma;
        rho = lb.rho;
    }
    return *this;
}

}
