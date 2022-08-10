/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   								   *
 *   luca.bisti@iet.unipi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef DEBORAHTANDEM_H
#define DEBORAHTANDEM_H

#include "node.h"
#include "flow.h"
#include "tnode.h"
#include "rng.h"
#include "deborah-globals.h"
#include "pseudoaffine.h"

#include <set>
#include <cstddef>
#include <memory>

namespace deborah
{

// This is the return value used to signal an error in LUDB-related functions
// (which deal with non-negative results only)
#define TANDEM_LUDB_ERROR_NESTINGTREE    -10000.0
#define TANDEM_LB_ERROR_INVALID_PARAM    -10000.0
#define TANDEM_LB_ERROR_INTERNAL    -20000.0

// When the following is defined, C++ exceptions handling is enabled in certain parts of the code (small performance penalty)
#define TANDEM_CATCH_EXCEPTIONS

// When the following is defined, the "--lb-random-combos" switch will be obeyed
#define TANDEM_LB_RANDOM_COMBOS

#define TANDEM_LUDB_MAXCUTS    0x7FFFFFFF

struct balanced_tree_descriptor
{
  int32_t order;
  int32_t depth;
  double node_latency_min;
  double node_latency_max;
  double node_rate_overprov_min;
  double node_rate_overprov_max;
  double flow_rate_min;
  double flow_rate_max;
  double flow_burst_min;
  double flow_burst_max;
};

// Also referred as cutset
typedef std::vector<uint64_t> CutNodeList;

struct OutputArrivalCurveParameters
{
  double burst;
  double rate;
};

// ----- Class declaration -----

class Tandem
{
 public:
  Tandem();

  Tandem(const Tandem &t2);

  Tandem(uint64_t nodes, uint64_t flows);

  Tandem(const std::vector<Node>& nodes, const std::vector<Flow>& flows);

  ~Tandem();

  bool CreateTree(const balanced_tree_descriptor &bt_desc, RNG &rng);

  bool CreateFullTandem(const balanced_tree_descriptor &bt_desc, RNG &rng);

  bool Load(const std::string &config_filename);

  bool Save(const std::string &config_filename);

  void Clear();

  void Print();

  bool BuildNestingTree();

  void PrintNestingTree();

  bool isSinkTree();

  bool isNested();

  bool checkProvisioning();

  std::tuple<int32_t, double, double> getUnderProvisionedNode();

  bool analyzeProvisioning(double &min, int32_t &nodemin, double &max, int32_t &nodemax);

  bool scaleProvisioning(double flow_mult, double node_mult);

  void PrintAllCombinations();

  void PrintSolution(simplex::Solution &solution, uint64_t dead_var);

  uint64_t LUDB_TotalCombos();

  Node *getNode(uint64_t index);

  Flow *getFlow(uint64_t index);

  Flow *getFlowId(uint64_t uid);

  uint64_t getFlowIdx(uint64_t uid);

  uint64_t NumNodes() const;

  uint64_t NumFlows() const;

  bool SetTaggedFlow(uint64_t tf);

  uint64_t GetTaggedFlow() const;

  // Non-nested tandem splitting methods
  std::vector<CutNodeList> compute_cutsets();

  void filter_cutsets(std::vector<CutNodeList> &cutsets, int32_t deltalen);

  std::vector<Tandem> cut(CutNodeList &cutset);

  int32_t join(std::vector<Tandem> &cuts);

  uint64_t aggregate_flows(bool merge_tf);

  bool compute_output_arrival_curve(uint64_t fid, Flow *fdest, int32_t mode);

  void output_arrival_curves(Tandem &t_next, int32_t mode);

  double LUDB(int32_t max_good_combos, bool deterministicOutput = false);

  double LUDB_quiet(int32_t max_good_combos, uint64_t &n_simplexes);

  double LUDB_quiet(const std::shared_ptr<TNode> &root, int32_t max_good_combos, uint64_t &n_simplexes);

  double LUDB_trivial();

  [[nodiscard]] double LUDB_getLast() const
  {
      return cached_ludb;
  };

  void LUDB_setLast(double ludb)
  {
      cached_ludb = ludb;
  };

  // CDF of foi computation, only for nested tandems
  OutputArrivalCurveParameters GetFoiOuputArrivalCurve(int32_t mode);

  PseudoAffine GetEquivalentServiceCurve(int32_t max_good_combos);

  // LowerBound computation
  double LowerBound(
      bool quiet = false, uint64_t combo_range_start = 0, uint64_t combo_range_len = 0,
      float percentage = 100.0, bool deterministic = false
  );

  // Optimized LowerBound computation with reduced number of combinations to try (UNSTABLE)
  double LowerBound_experimental(float percentage = 100.0);

  double LowerBound_experimental_quiet(uint64_t &num_combos, float percentage = 100.0);

  bool compute_output_rates(uint64_t tagged); // temporarily public method
  uint64_t
  compute_lowerbound_rules(uint64_t tagged, std::vector<std::byte> &burst_time); // temporarily public method

 private:
  uint64_t num_nodes;
  uint64_t num_flows;
  uint64_t tagged_flow;

  std::vector<Node> nodes;
  std::vector<Flow> flows;

  void ComputeNestingLevels();

  int32_t nesting_level;

  std::shared_ptr<TNode> nesting_tree;

  std::shared_ptr<TNode> build_nesting_subtree(uint64_t flow_id);

  double compute_delay(std::vector<bool> &combo);

  bool bit(uint64_t x, uint64_t n);

  bool compute_CDF_i(uint64_t node_id, uint64_t flow_id);

  bool isNodeCrossedByFlow(uint64_t node, uint64_t flow);

  bool isNodeEnteredByFlow(uint64_t node, uint64_t flow);

  bool isNodeLeftByFlow(uint64_t node, uint64_t flow);

  bool isInterferingFlow(uint64_t node, uint64_t flow, uint64_t tagged);

#define MAXDEPS 32
  struct depsse
  {
    uint64_t dep[MAXDEPS];
  };

  struct compstatus
  {
    uint64_t cres;
    depsse cdep;
    uint64_t cnode;
  };

  bool
  cutcomp_slow(bool *res, uint64_t ndeps, std::set<uint64_t> nodedeps[], std::set<std::set<uint64_t> > &allres);

  void cutcomp(compstatus *stat, std::vector<depsse> deps, std::set<std::set<uint64_t> > &allres);

  // Cached data
  double cached_ludb;
};

}
#endif
