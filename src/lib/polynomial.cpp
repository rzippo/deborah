#include "polynomial.h"
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

namespace deborah
{

void Polynomial::Zero()
{
    std::fill(variableCoefficients.begin(), variableCoefficients.end(), 0.0);
    constant = 0.0;
}

Polynomial::Polynomial(uint64_t num_variables)
{
    variableCoefficients = std::vector<double>(num_variables);
    Zero();
}

Polynomial::Polynomial(std::vector<double> init_coefficients, double init_constant)
{
    variableCoefficients = init_coefficients;
    constant = init_constant;
}

Polynomial::Polynomial(const Polynomial &p)
{
    variableCoefficients = p.variableCoefficients;
    constant = p.constant;
}

uint64_t Polynomial::getNumVariables() const
{
    return variableCoefficients.size();
}

/*
Polynomial::Polynomial(uint64_t num_variables, LinearSegment &ls)
{
	num_coeff = num_variables + 1;
	coeff = new double [num_coeff];
	Zero();
	if(ls.y == 0.0)
		SetB(ls.x);
	else SetB(ls.y);
	
}
*/

void Polynomial::SetCoefficient(uint64_t idx, double a)
{
    variableCoefficients[idx] = a;
}

void Polynomial::SumVariableCoefficient(uint64_t idx, double a)
{
    variableCoefficients[idx] += a;
}

std::vector<double> Polynomial::GetCoefficients() const
{
    return std::vector<double>(variableCoefficients);
}

double Polynomial::getCoefficient(uint64_t idx) const
{
    return variableCoefficients[idx];
}

void Polynomial::SetConstant(double b)
{
    constant = b;
}

double Polynomial::getConstant() const
{
    return constant;
}

bool Polynomial::isConstant()
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

Polynomial &Polynomial::operator=(const Polynomial &p)
{
    if (&p != this)
    {
        variableCoefficients = p.variableCoefficients;
        constant = p.constant;
    }
    return *this;
}

Polynomial Polynomial::operator+(const Polynomial &p)
{
    Polynomial result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] + p.variableCoefficients[i];
    }

    result.constant = this->constant + p.constant;

    return result;
}

Polynomial Polynomial::operator+(double b)
{
    Polynomial result = *this;
    result.constant += b;

    return result;
}

Polynomial Polynomial::operator-(const Polynomial &p)
{
    Polynomial result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] - p.variableCoefficients[i];
    }

    result.constant = this->constant - p.constant;

    return result;
}

Polynomial Polynomial::operator-(double b)
{
    Polynomial result = *this;
    result.constant -= b;

    return result;
}

Polynomial Polynomial::operator*(const double f)
{
    Polynomial result(variableCoefficients.size());

    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        result.variableCoefficients[i] = this->variableCoefficients[i] * f;
    }

    result.constant = this->constant * f;

    return result;
}

Polynomial Polynomial::operator/(const double f)
{
    Polynomial result(variableCoefficients.size());

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

void Polynomial::Negate()
{
    for (int32_t i = 0; i < variableCoefficients.size(); i++)
    {
        variableCoefficients[i] = -variableCoefficients[i];
    }

    constant = -constant;
}

void Polynomial::Print()
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

std::string Polynomial::PrintJson()
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
