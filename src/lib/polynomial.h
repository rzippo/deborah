#ifndef DEBORAHPOLY_H
#define DEBORAHPOLY_H

#include <vector>
#include <string>
#include "linearsegment.h"

namespace deborah
{

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class Polynomial
{
 public:
  explicit Polynomial(uint64_t num_variables);

  Polynomial(std::vector<double> init_coefficients, double init_constant);

  Polynomial(const Polynomial &p);
  //Polynomial(uint64_t num_variables, LinearSegment &s);

  void Zero();

  std::vector<double> GetCoefficients() const;

  double getCoefficient(uint64_t idx) const;

  void SetCoefficient(uint64_t idx, double a);

  double getConstant() const;

  void SetConstant(double b);

  void SumVariableCoefficient(uint64_t idx, double a);

  uint64_t getNumVariables() const;

  bool isConstant();

  void Negate();

  Polynomial operator+(const Polynomial &p);

  Polynomial operator+(double b);

  Polynomial operator-(const Polynomial &p);

  Polynomial operator-(double b);

  Polynomial operator*(double f);

  Polynomial operator/(double f);

  Polynomial &operator=(const Polynomial &p);

  void Print();

  std::string PrintJson();

 private:
  std::vector<double> variableCoefficients;
  double constant;
};

}

#endif
