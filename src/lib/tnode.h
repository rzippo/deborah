#ifndef DEBORAHTNODE_H
#define DEBORAHTNODE_H

// enable the following to dump the PI and PI_EQ curves at each tnode during computation of LUDB
//#define TNODE_VERBOSE_LUDB

// enable this to print debug messages in the experimental LUDB computation function
//#define TNODE_DEBUG_LUDB_EXPERIMENTAL

// enable this to dump the constraints set for each simplex performed
//#define TNODE_DEBUG_LUDB_CONSTRAINTS


#define TNODE_LUDB_EPSILON    0.0000001
#define TNODE_LUDB_INFINITY    9999999999.0

#define LUDB_MODE_EXACT        -100
#define LUDB_MODE_HEURISTIC 0x7FFFFFFF

#include "flow.h"
#include "node.h"
#include "curve.h"
#include "linearsegment.h"
#include "parametrizedpseudoaffine.h"

#include "simplex/constraint.h"
#include "simplex/objective.h"
#include "simplex/solution.h"

#include <vector>
#include <memory>

namespace deborah
{

typedef std::vector<uint64_t> IndexSet;
//typedef std::vector<simplex::Constraint> ConstraintSet;
typedef std::vector<IndexSet> IndexSetVector;

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class TNode
{
 public:
  TNode();

  TNode(uint64_t num_flows);

  TNode(uint64_t flow_id, uint64_t num_tnodes, uint64_t num_leaves);

  bool has_leaf();

  uint64_t num_children();

  ~TNode();

  // computes a specific instance of the PI
  ParametrizedPseudoAffine &getPI(const IndexSet &idx_set, std::vector<simplex::Constraint> &c_set);

  // computes the equivalent service curve for the given PI
  ParametrizedPseudoAffine &computeEquivalent(ParametrizedPseudoAffine &p, const IndexSet &set, std::vector<simplex::Constraint> &c_set);

  // compute a specific instance of the equivalent service curve at this node
  ParametrizedPseudoAffine computeEquivalentServiceCurve(const IndexSet &idx_set);

  // compute the Delay Bound at this node
  simplex::Objective& getDelayBound(const IndexSet &idx_set, std::vector<simplex::Constraint> &constraints);

  // compute LUDB (old-fashioned algorithm, serial enumeration)
  double computeLUDB(simplex::Objective &delay, simplex::Solution &solution, int32_t nvars, IndexSetVector *good_combos);

  // compute LUDB (new implementation, recursive enumeration)
  uint64_t computeLUDB(
      simplex::Objective &delay, simplex::Solution &solution, int32_t nvars, IndexSetVector &good_combos,
      int32_t max_good_combos, bool eta
  );

  void Print(bool subtree);

  void PrintIndexSet(const IndexSet &set, bool compact);

  // index and pointer to a Flow object, representing flow (h,k)
  uint64_t flow;
  Flow *pflow;

  // array of TNode children, representing the set S(h,k)
  std::vector<std::shared_ptr<TNode>> children;
  uint64_t num_tnodes;

  // array of node indexes, representing set C(h,k)
  uint64_t *leaves;
  uint64_t num_leaves;

  // The PI(h,k) pseudo-affine service curve
  ParametrizedPseudoAffine *pi;
  ParametrizedPseudoAffine *pi_eq;
  // The PI_C(h,k) rate-latency service curve, obtained as (X)Bj
  LinearSegment pi_c;

  // todo: what is an 'index'?
  // Structured set of sets of indexes for the current node
  IndexSet index_set;

  // Each element is the number of indexes belonging to the corresponding node
  // the first node is the current one, then children follow
  std::vector<uint64_t> index_schema;

  // Each element is the TNode number associated to the corresponding index
  IndexSet index_nodes;

  // Methods devoted to handle the recursive set data structure
  bool getCombo(IndexSet &set, uint64_t n);

  int32_t incCombo(IndexSet &set, int32_t index);

  uint64_t getNumCombos();

  // Return the Ith subset of the node's index set
  bool getSubset(IndexSet &set, int32_t n_child);

  bool getSubset(const IndexSet &set, int32_t n_child, IndexSet &subset);

};

}

#endif
