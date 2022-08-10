#ifndef DEBORAHPSEUDOAFFINE_H
#define DEBORAHPSEUDOAFFINE_H

#include "polynomial.h"
#include "linearsegment.h"
#include "parametrizedleakybucket.h"

#include "simplex/objective.h"

#include <vector>

using namespace deborah;

namespace deborah
{

/*
typedef struct ParametrizedLeakyBucket {
    Polynomial *sigma;
    double rho;
} ParametrizedLeakyBucket;
*/

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class ParametrizedPseudoAffine
{
 public:
  explicit ParametrizedPseudoAffine(uint64_t variables);

  bool MakeLatencyRate(LinearSegment &ls);

  bool Convolve(ParametrizedPseudoAffine &p);

  void Zero();

  void AddStage(ParametrizedLeakyBucket &lb);

  ParametrizedPseudoAffine& operator=(const ParametrizedPseudoAffine &p);

  ~ParametrizedPseudoAffine();

  void Print();

  std::string PrintJson();

  /// Delay expression, value depends on s parameters
  simplex::Objective delay;

  /// Leaky bucket stages
  std::vector<ParametrizedLeakyBucket> stages;

 private:
  uint64_t poly_variables;

};

}

#endif
