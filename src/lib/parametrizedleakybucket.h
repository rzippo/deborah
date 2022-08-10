#ifndef DEBORAHLEAKYBUCKET_H
#define DEBORAHLEAKYBUCKET_H

#include "polynomial.h"

namespace deborah
{

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class ParametrizedLeakyBucket
{
 public:
  explicit ParametrizedLeakyBucket(uint64_t num_variables);

  ParametrizedLeakyBucket(const ParametrizedLeakyBucket &lb);

  ~ParametrizedLeakyBucket();

  ParametrizedLeakyBucket &operator=(const ParametrizedLeakyBucket &lb);

  /// Burst expression, value depends on s parameters
  Polynomial sigma;

  /// Long term rate, does not depend on s parameters
  double rho;

  uint64_t nvars;
};

}

#endif
