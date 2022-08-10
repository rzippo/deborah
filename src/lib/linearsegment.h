#ifndef DEBORAHLINEARSEGMENT_H
#define DEBORAHLINEARSEGMENT_H

#include <float.h>
#include <cstdint>

namespace deborah
{

/**
This class represent a linear segment of a curve
For instance, it can represent a latency-rate service curve

	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/


#define LINEARSEGMENT_POSITIVE_INFINITE DBL_MAX
#define LINEARSEGMENT_NEGATIVE_INFINITE -DBL_MAX
#define LINEARSEGMENT_EPSILON 0.000000001

class LinearSegment
{
 public:
  LinearSegment();

  LinearSegment(double x0, double y0, double r);

  LinearSegment(double x0, double r);

  LinearSegment(const LinearSegment &ls);

  ~LinearSegment();

  LinearSegment &operator=(const LinearSegment &s2);

  void setLatencyRate(double l, double r);

  void setAffine(double b, double r);

  void add(LinearSegment &s2, double xpos, bool leftopen);

  void sub(LinearSegment &s2, double xpos, bool leftopen);

  void Zero();

  //Return the value of the function at the given x-coordinate
  double value(double x1);

  //Return the x-coordinate of the intersection between this segment and the specified one
  double get_X_intersection(LinearSegment &ls);

  //Modifies the current segment by convolving it with the specified one
  bool convolve(LinearSegment &ls);

  bool convolve(double latency, double rate);

  double x;
  double y;
  double rate;
  bool leftopen;
  int32_t id;
  double gap;

  void Print();

};

}

#endif
