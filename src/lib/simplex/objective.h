#ifndef SIMPLEX_OBJECTIVE_H
#define SIMPLEX_OBJECTIVE_H

#include "../polynomial.h"

#include <vector>

namespace simplex
{

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class Objective
{
 public:
  explicit Objective(uint64_t num_variables);

  Objective(std::vector<double> init_coefficients, double init_constant);

  explicit Objective(const deborah::Polynomial &p);
  //Polynomial(uint64_t num_variables, LinearSegment &s);

  void Zero();

  double getCoefficient(uint64_t idx) const;

  void setCoefficient(uint64_t idx, double a);

  [[nodiscard]]
  double getConstant() const;

  void setConstant(double b);

  void SumVariableCoefficient(uint64_t idx, double a);

  uint64_t getNumVariables() const;

  bool isConstant();

  void Negate();

  Objective operator+(const Objective &p);
  Objective operator+(const deborah::Polynomial &p);
  Objective operator+(double b);

  Objective operator-(const Objective &p);

  Objective operator-(double b);

  Objective operator*(double f);

  Objective operator/(double f);

  Objective &operator=(const Objective &p);

  void Print();

  std::string PrintJson();

 private:
  std::vector<double> variableCoefficients;
  double constant;
};

}

#endif //SIMPLEX_OBJECTIVE_H
