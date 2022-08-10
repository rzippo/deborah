/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   *
 *   luca.bisti@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "parametrizedpseudoaffine.h"

#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

namespace deborah
{

ParametrizedPseudoAffine::ParametrizedPseudoAffine(uint64_t num_variables) : delay(num_variables)
{
    poly_variables = num_variables;
    stages.clear();
}

ParametrizedPseudoAffine::~ParametrizedPseudoAffine()
{
    stages.clear();
}

ParametrizedPseudoAffine& ParametrizedPseudoAffine::operator=(const ParametrizedPseudoAffine &p)
{
    if (&p != this)
    {
        Zero();
        delay = p.delay;
        for (int32_t i = 0; i < stages.size(); i++)
        {
            stages.push_back(p.stages[i]);
        }
    }
    return *this;
}

bool ParametrizedPseudoAffine::MakeLatencyRate(LinearSegment &ls)
{
    stages.clear();
    delay.Zero();
    delay.setConstant(ls.x);
    ParametrizedLeakyBucket lb(poly_variables);
    lb.rho = ls.rate;
    stages.push_back(lb);
    return true;
}

/* Convolve the current curve with the passed one:
 * sum delays and add leaky bucket stages
 */
bool ParametrizedPseudoAffine::Convolve(ParametrizedPseudoAffine &p)
{
    delay = delay + p.delay;
    for (const auto& stage : p.stages)
    {
        stages.push_back(stage);
    }

    return true;
}

/// Clears this ParametrizedPseudoAffine
void ParametrizedPseudoAffine::Zero()
{
    stages.clear();
    delay.Zero();
}

void ParametrizedPseudoAffine::AddStage(ParametrizedLeakyBucket &lb)
{
    stages.push_back(lb);
}

void ParametrizedPseudoAffine::Print()
{
    std::cout << stages.size() << " leaky-bucket stage(s):" << std::endl;

    std::cout << "Delay = ";
    delay.Print();
    std::cout << std::endl;

    for (int32_t i = 0; i < stages.size(); i++)
    {
        std::cout << "Stage " << std::setw(2) << std::setfill('0') << i << ": sigma = ";
        stages[i].sigma.Print();
        std::cout
            << std::endl
            << "            rho = " << std::fixed << std::setprecision(2) << stages[i].rho
            << std::endl;
    }
}

std::string ParametrizedPseudoAffine::PrintJson()
{
    std::stringstream ss;

    ss
        << "\t{"
        << " \"delay\": " << delay.PrintJson() << ","
        << " \"stages\" : [ " ;

    auto lbs = stages.size();
    for (int lb_idx = 0; lb_idx < lbs; ++lb_idx)
    {
        auto lb = stages[lb_idx];
        ss
            << "{ \"rho\": " << std::fixed << std::setprecision(6) << lb.rho << ","
            << " \"sigma\": " << lb.sigma.PrintJson() << " }";

        if (lbs > 1 && lb_idx < lbs - 1)
        {
            ss << " , ";
        }
    }

    ss << " ] }";

    return ss.str();
}

}
