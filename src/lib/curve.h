#ifndef DEBORAHCURVE_H
#define DEBORAHCURVE_H

#include <cmath>
#include <vector>
#include "linearsegment.h"

namespace deborah
{

/**
This class implements a piecewise linear curve object and some common operations
that can be performed on it, such as convolution, summation, etc.
The curve is assumed to be monotone non-decreasing.

	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/

class Curve
{
 public:

  Curve();

  Curve(int32_t nsegs);

  ~Curve();

  int32_t num_segments() const;

  void Zero();

  bool add(Curve &c2); //TODO: rewrite
  bool add_segment(LinearSegment &s);

  bool add_segment(double x, double y, double rate);

  bool create_latency_rate(double latency, double rate);

  bool create_line_burst(double x1, double x2, double rate, double burst, uint64_t flow_id);

  bool create_burst_line(double x1, double x2, double rate, double burst, uint64_t flow_id);

  bool create_token_bucket(double burst, double rate);

  double f(double x1);

  double f(double x1, bool right_limit);

  double f_inv(double y1, bool rightmost);

  double f_inv2(double y1, bool rightmost);

  int32_t getSegmentDefining(double x);

  double getFirstBit(uint64_t flow);

  double getLastBit(uint64_t flow);

  bool ConvolveWithLatencyRate(Curve &beta);

  double getLatency();

  bool removeLatency();

  bool shiftRight(double L);

  LinearSegment &getSeg(uint64_t s);

  double getSlope(double x_start, bool right);

  double getTotalTraffic();

  int32_t has_BPX_at(double x);

  int32_t find_gaps_at(double x, int32_t &num_gaps);

  int32_t has_BPY_at(double y);

  int32_t add_BPX(double x);

  int32_t add_BPY(double y);

  bool recompute_slopes();

  bool remove_colinear_segments();

  bool remove_infinitesimal_segments();

  bool move_burst_id(uint64_t flow_id, bool forward);

  Curve &operator=(const Curve &c2);

  static const int32_t Error = -2147483647 - 1;

  static const double Epsilon;
  static const double PlusInfinity;
  static const double MinusInfinity;

  /// True if x is within Curve::Epsilon from 0
  static bool IsZero(double x)
  {
      return std::fabs(x) < Epsilon;
  }

  /// True if the values are within Curve::Epsilon from each other
  static bool AreEqual(double x, double y)
  {
      return IsZero(x - y);
  }

  /// True if x is above y by more than Curve::Epsilon
  static bool IsGreater(double x, double y)
  {
      return x > y + Epsilon;
  }

  /// True if x is below y by more than Curve::Epsilon
  static bool IsLower(double x, double y)
  {
      return x < y - Epsilon;
  }

  /// True if x is Double::PlusInfinity
  static bool IsPlusInfinity(double x)
  {
      return x == PlusInfinity;
  }

  /// True if x is Double::MinusInfinity
  static bool IsMinusInfinity(double x)
  {
      return x == MinusInfinity;
  }

  void Print();

  bool SanityCheck();

 private:
  std::vector<LinearSegment> segments;

  int32_t getSegmentFirstAtValue(double y1);

  bool convolution(Curve &c2);
};

}

#endif
