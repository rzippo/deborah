#include "solution.h"

#include <iostream>
#include <stdio.h>
#include <iomanip>

namespace simplex
{

void Solution::Zero()
{
    std::fill(variableCoefficients.begin(), variableCoefficients.end(), 0.0);
    constant = 0.0;
}

Solution::Solution(uint64_t num_variables)
{
    variableCoefficients = std::vector<double>(num_variables);
    Zero();
}

Solution::Solution(std::vector<double> init_coefficients, double init_constant)
{
    variableCoefficients = init_coefficients;
    constant = init_constant;
}

Solution::Solution(const deborah::Polynomial &p)
{
    variableCoefficients = p.GetCoefficients();
    constant = p.getConstant();
}

uint64_t Solution::getNumVariables() const
{
    return variableCoefficients.size();
}

/*
Solution::Polynomial(uint64_t num_variables, LinearSegment &ls)
{
	num_coeff = num_variables + 1;
	coeff = new double [num_coeff];
	Zero();
	if(ls.y == 0.0)
		SetB(ls.x);
	else SetB(ls.y);
	
}
*/

void Solution::setCoefficient(uint64_t idx, double a)
{
    variableCoefficients[idx] = a;
}

void Solution::SumVariableCoefficient(uint64_t idx, double a)
{
    variableCoefficients[idx] += a;
}

double Solution::getCoefficient(uint64_t idx) const
{
    return variableCoefficients[idx];
}

void Solution::setConstant(double b)
{
    constant = b;
}

double Solution::getConstant() const
{
    return constant;
}

bool Solution::isConstant()
{
    for (double variableCoefficient : variableCoefficients)
    {
        if (variableCoefficient != 0.0)
        {
            return false;
        }
    }

    return true;
}

Solution &Solution::operator=(const Solution &p)
{
    if (&p != this)
    {
        variableCoefficients = p.variableCoefficients;
        constant = p.constant;
    }
    return *this;
}

Solution Solution::operator+(const Solution &p)
{
    Solution result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] + p.variableCoefficients[i];
    }

    result.constant = this->constant + p.constant;

    return result;
}

Solution Solution::operator+(double b)
{
    Solution result = *this;
    result.constant += b;

    return result;
}

Solution Solution::operator-(const Solution &p)
{
    Solution result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] - p.variableCoefficients[i];
    }

    result.constant = this->constant - p.constant;

    return result;
}

Solution Solution::operator-(double b)
{
    Solution result = *this;
    result.constant -= b;

    return result;
}

Solution Solution::operator*(const double f)
{
    Solution result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] * f;
    }

    result.constant = this->constant * f;

    return result;
}

Solution Solution::operator/(const double f)
{
    Solution result(variableCoefficients.size());

    if (f == 0.0)
    {
        for (int32_t i = 0; i < variableCoefficients.size(); i++)
        {
            result.variableCoefficients[i] = 999999999999.9999999999999;
        }

        result.constant = 999999999999.9999999999999;
    }
    else
    {
        for (int32_t i = 0; i < variableCoefficients.size(); i++)
        {
            result.variableCoefficients[i] = this->variableCoefficients[i] / f;
        }

        result.constant = this->constant / f;
    }

    return result;
}

void Solution::Negate()
{
    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        variableCoefficients[i] = -variableCoefficients[i];
    }

    constant = -constant;
}

void Solution::Print()
{
    //if(num_coeff == 0) return;
    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        std::cout
            << std::fixed << std::setprecision(2) << variableCoefficients[i]
            << "*X" << i << " + ";
    }

    std::cout << std::fixed << std::setprecision(2) << constant << std::flush;
}

}
