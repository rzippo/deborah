//
// Created by raffa on 22/02/2021.
//

#ifndef DEBORAH_SRC_LIB_PSEUDOAFFINE_H_
#define DEBORAH_SRC_LIB_PSEUDOAFFINE_H_

#include <vector>
#include <string>

#include "leakybucket.h"
#include "parametrizedpseudoaffine.h"

namespace deborah
{

class PseudoAffine
{
 public:
  PseudoAffine(std::vector<LeakyBucket> stages, double delay);
  PseudoAffine(const ParametrizedPseudoAffine& ppa, const simplex::Solution& sol);

  [[nodiscard]] std::string PrintJson() const;

 private:
  std::vector<LeakyBucket> stages;
  double delay;
};

}

#endif //DEBORAH_SRC_LIB_PSEUDOAFFINE_H_
