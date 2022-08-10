#include "objective.h"

#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

namespace simplex
{

void Objective::Zero()
{
    std::fill(variableCoefficients.begin(), variableCoefficients.end(), 0.0);
    constant = 0.0;
}

Objective::Objective(uint64_t num_variables)
{
    variableCoefficients = std::vector<double>(num_variables);
    Zero();
}

Objective::Objective(std::vector<double> init_coefficients, double init_constant)
{
    variableCoefficients = init_coefficients;
    constant = init_constant;
}

Objective::Objective(const deborah::Polynomial &p)
{
    variableCoefficients = p.GetCoefficients();
    constant = p.getConstant();
}

uint64_t Objective::getNumVariables() const
{
    return variableCoefficients.size();
}

/*
Objective::Polynomial(uint64_t num_variables, LinearSegment &ls)
{
	num_coeff = num_variables + 1;
	coeff = new double [num_coeff];
	Zero();
	if(ls.y == 0.0)
		SetB(ls.x);
	else SetB(ls.y);
	
}
*/

void Objective::setCoefficient(uint64_t idx, double a)
{
    variableCoefficients[idx] = a;
}

void Objective::SumVariableCoefficient(uint64_t idx, double a)
{
    variableCoefficients[idx] += a;
}

double Objective::getCoefficient(uint64_t idx) const
{
    return variableCoefficients[idx];
}

void Objective::setConstant(double b)
{
    constant = b;
}

double Objective::getConstant() const
{
    return constant;
}

bool Objective::isConstant()
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

Objective &Objective::operator=(const Objective &p)
{
    if (&p != this)
    {
        variableCoefficients = p.variableCoefficients;
        constant = p.constant;
    }
    return *this;
}

Objective Objective::operator+(const Objective &p)
{
    Objective result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] + p.variableCoefficients[i];
    }

    result.constant = this->constant + p.constant;

    return result;
}

Objective Objective::operator+(const deborah::Polynomial &p)
{
    Objective result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] + p.getCoefficient(i);
    }

    result.constant = this->constant + p.getConstant();

    return result;
}

Objective Objective::operator+(double b)
{
    Objective result = *this;
    result.constant += b;

    return result;
}

Objective Objective::operator-(const Objective &p)
{
    Objective result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] - p.variableCoefficients[i];
    }

    result.constant = this->constant - p.constant;

    return result;
}

Objective Objective::operator-(double b)
{
    Objective result = *this;
    result.constant -= b;

    return result;
}

Objective Objective::operator*(const double f)
{
    Objective result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] * f;
    }

    result.constant = this->constant * f;

    return result;
}

Objective Objective::operator/(const double f)
{
    Objective result(variableCoefficients.size());

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

void Objective::Negate()
{
    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        variableCoefficients[i] = -variableCoefficients[i];
    }

    constant = -constant;
}

void Objective::Print()
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

std::string Objective::PrintJson()
{
    std::stringstream ss;

    ss << "{ "
       << " \"constant\": " << std::fixed << std::setprecision(6) << getConstant() << ","
        << " \"coefficients\": [ ";

    auto vars = getNumVariables();
    for (int var_idx = 0; var_idx < vars; ++var_idx)
    {
        ss << std::fixed << std::setprecision(6) << getCoefficient(var_idx);

        if (vars > 1 && var_idx < vars - 1)
        {
            ss << " , ";
        }
    }

    ss << " ]}";

    return ss.str();
}

}
