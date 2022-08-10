//
// Created by raffa on 22/02/2021.
//

#include "pseudoaffine.h"

#include <utility>
#include <sstream>
#include <iomanip>

PseudoAffine::PseudoAffine(std::vector<LeakyBucket> stages, double delay)
{
    this->stages = std::move(stages);
    this->delay = delay;
}

PseudoAffine::PseudoAffine(const ParametrizedPseudoAffine &ppa, const simplex::Solution &sol)
{
    double d = ppa.delay.getConstant();
    for (int var_idx = 0; var_idx < ppa.delay.getNumVariables(); ++var_idx)
    {
        d += ppa.delay.getCoefficient(var_idx) * sol.getCoefficient(var_idx);
    }

    this->delay = d;

    for (const auto& plb : ppa.stages)
    {
        this->stages.emplace_back(plb, sol);
    }
}

std::string PseudoAffine::PrintJson() const
{
    std::stringstream ss;
    ss
        << "{"
        << "\"delay\":" << std::fixed << std::setprecision(6) << delay << ","
        << "\"stages\":[";

    auto lbs = stages.size();
    for (int lb_idx = 0; lb_idx < lbs; ++lb_idx)
    {
        auto lb = stages[lb_idx];
        ss << lb.PrintJson();

        if (lbs > 1 && lb_idx < lbs - 1)
            ss << " , ";
    }
    ss << "]}";

    return ss.str();
}
