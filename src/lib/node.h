#ifndef DEBORAHNODE_H
#define DEBORAHNODE_H

#include "curve.h"

namespace deborah
{

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/

class Node
{
 public:
  Node();

  ~Node();

  Node(uint64_t src, uint64_t exit, double theta, double R);

  double latency;
  double rate;

  void clear_curves(uint64_t n_flows);

  //Cumulative Departure Function of all traffic from the node.
  Curve CDF;

  //Cumulative Arrival Function of all traffic to the node.
  Curve CAF;

  //Cumulative Departure Functions of each flow from the node.
  //Flows not traversing the node have Zero().
  Curve *CDF_i;

  //Cumulative Arrival Functions of each flow to the node.
  //Flows not traversing the node have Zero().
  Curve *CAF_i;

  //Service function of the node
  Curve beta;

  bool set_output_rates(uint64_t num_flows, std::vector<double> out_rates, bool rate_red);

  double get_output_rate(uint64_t flow);

  bool set_output_rate(uint64_t flow, double rate);

  bool is_rate_reducing();

 private:
  double *max_output_rates;
  bool rate_reduction;

};

}

#endif
