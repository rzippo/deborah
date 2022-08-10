#ifndef SIMPLEX_SOLUTION_H
#define SIMPLEX_SOLUTION_H

#include <vector>
//#include "../linearsegment.h"
#include "../polynomial.h"

namespace simplex
{

class Solution
{
 public:
  explicit Solution(uint64_t num_variables);

  Solution(std::vector<double> init_coefficients, double init_constant);

  explicit Solution(const deborah::Polynomial &p);
  //Polynomial(uint64_t num_variables, LinearSegment &s);

  void Zero();

  double getCoefficient(uint64_t idx) const;

  void setCoefficient(uint64_t idx, double a);

  double getConstant() const;

  void setConstant(double b);

  void SumVariableCoefficient(uint64_t idx, double a);

  uint64_t getNumVariables() const;

  bool isConstant();

  void Negate();

  Solution operator+(const Solution &p);

  Solution operator+(double b);

  Solution operator-(const Solution &p);

  Solution operator-(double b);

  Solution operator*(double f);

  Solution operator/(double f);

  Solution &operator=(const Solution &p);

  void Print();

 private:
  std::vector<double> variableCoefficients;
  double constant;
};

}

#endif //SIMPLEX_SOLUTION_H
