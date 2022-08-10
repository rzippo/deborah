#pragma once

#include "lib/tandem.h"

extern uint32_t RNG_SEED;

bool MakeTreeConfig(char *filename, const int32_t order, const int32_t depth);

bool MakeFullTandemConfig(char *filename, const int32_t num_nodes, const int32_t flows_percent);

bool LUDB_evaluation(deborah::Tandem &tandem, int32_t order, int32_t depth, int32_t max_good_combos);

bool LowerBound_evaluation(deborah::Tandem &tandem, int32_t order, int32_t depth, bool experimental, float percentage);

double LUDB_non_nested_MSA(
    deborah::Tandem &tandem, int32_t mode, std::vector<CutNodeList> *precomputed_cutsets, uint64_t &totsimpl,
    uint64_t &num_cuts
);

double LUDB_non_nested_STA(
    deborah::Tandem &tandem, int32_t mode, std::vector<CutNodeList> *precomputed_cutsets, uint64_t &totsimpl,
    uint64_t &num_cuts
);

double
OutputCurve_non_nested_STA(
    deborah::Tandem &tandem, int32_t mode, std::vector<CutNodeList> *precomputed_cutsets, uint64_t &totsimpl,
    uint64_t &num_cuts
);

double
EquivalentServiceCurve_non_nested_STA(
    deborah::Tandem &tandem, int32_t mode, std::vector<CutNodeList> *precomputed_cutsets, uint64_t &totsimpl,
    uint64_t &num_cuts
);

bool
LUDB_non_nested_evaluation(int32_t nnodes, int32_t flows_percent, int32_t mode, uint64_t maxcuts, bool algo_classic);

double DelayBound_per_node(deborah::Tandem &t1, int32_t ludb_mode);
