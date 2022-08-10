//
// Created by raffa on 22/02/2021.
//

#ifndef DEBORAH_SRC_LIB_LEAKYBUCKET_H_
#define DEBORAH_SRC_LIB_LEAKYBUCKET_H_

#include <string>

#include "parametrizedleakybucket.h"
#include "simplex/solution.h"

namespace deborah {

/// Non-parametric LeakyBucket
/// Can be computed from a parametric l.b. and a set of parameter values
class LeakyBucket
{
public:
    LeakyBucket(double sigma, double rho);
    LeakyBucket(const ParametrizedLeakyBucket& plb, const simplex::Solution& sol);

    [[nodiscard]] double getSigma() const;
    [[nodiscard]] double getRho() const;
    [[nodiscard]] std::string PrintJson() const;

private:
    double sigma;
    double rho;
};

}

#endif //DEBORAH_SRC_LIB_LEAKYBUCKET_H_
