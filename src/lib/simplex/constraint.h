#ifndef SIMPLEX_CONSTRAINT_H
#define SIMPLEX_CONSTRAINT_H

#include "../polynomial.h"

namespace simplex
{

#define CPOLY_TNODEID_NULL    0xFFFFFFFF

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class Constraint
{
 public:
  Constraint(uint64_t num_variables);

  Constraint(uint64_t num_variables, uint64_t tnode_id);

  Constraint(const Constraint &cp);

  ~Constraint();

  Constraint &operator=(const Constraint &cp);

  Constraint &operator=(const deborah::Polynomial &cp);

  bool isOxymoron();

  bool isTautology();

  void InvertDirection();

  bool isSatisfied(deborah::Polynomial &p);

  void Print();

  deborah::Polynomial Inequality;
  bool greater_than_zero;
  uint64_t nvars;
  uint64_t tnode_id;
};

}

#endif //SIMPLEX_CONSTRAINT_H
