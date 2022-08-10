/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   				   *
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "deborah-core.h"
#include "lib/simplex/simplex.h"
#include "lib/timing.h"

using namespace std;
using namespace deborah;

// Enable the following to split benchmarking for cuts computation and LUDB processing
//#define LUDB_NNESTED_EVAL_BENCHMARK

// Initial seed for Random Number Generator
uint32_t RNG_SEED = 1000;
// Random Number Generator (global)
static RNG rng(RNG_SEED);

// T-Student coefficients
#define T_STUDENT_SIZE 30
double t_student95[T_STUDENT_SIZE] = {12.7062, 4.3027, 3.1825, 2.7765, 2.5706, 2.4469, 2.3646, 2.3060, 2.2622, 2.2281,
                                      2.2010, 2.1788, 2.1604, 2.1448,
                                      2.1315, 2.1199, 2.1098, 2.1009, 2.0930, 2.0860, 2.0796, 2.0739, 2.0687, 2.0639,
                                      2.0595, 2.0555, 2.0518, 2.0484, 2.0452, 2.0423};

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*	This method computes the confidence interval for the 95th percentile
	Stores the mean value and the semi-width of the confidence interval (delta)
*/
bool confidence95(std::vector<double> &vett, const int32_t n, double &mean, double &delta)
{
    double var;
    double tstud_coeff;

    if (n < 2)
    {
        return false;
    }
    if (n - 1 >= T_STUDENT_SIZE)
    {
        tstud_coeff = t_student95[T_STUDENT_SIZE - 1];
    }
    else
    {
        tstud_coeff = t_student95[n - 2];
    }

    // compute mean
    mean = 0.0;
    for (int32_t i = 0; i < n; i++)
    {
        mean += vett[i];
    }
    mean /= n;

    // compute variance
    var = 0.0;
    for (int32_t i = 0; i < n; i++)
    {
        var += pow(vett[i] - mean, 2.0);
    }
    var /= n - 1;

    // compute semi-interval
    delta = tstud_coeff * sqrt(var / n);

    return true;
}

/*	This method loads an array of LUDB values and simplex counts
	to avoid recomputing some stats
*/
int32_t LUDB_LoadPrecomputedValues(
    const std::string &filename,
    std::vector<double> &ludb,
    std::vector<uint64_t> &simplexes,
    int32_t num_trials
)
{
    int32_t trial;
    int32_t n;
    for (int32_t i = 0; i < num_trials; i++)
    {
        ludb[i] = 0.0;
        simplexes[i] = 0;
    }

    std::ifstream f;
    f.open(filename, std::ifstream::in | std::ifstream::binary);
    if (f.fail())
    {
        return -1;
    }

    for (trial = 0; trial < num_trials && !f.eof(); trial++)
    {
        f >> ludb[trial] >> simplexes[trial];
    }

    f.close();
    return trial;
}

/*	This method dumps an array of LUDB values and simplex counts to disk, for reuse
*/
bool LUDB_SavePrecomputedValues(
    const std::string &filename, std::vector<double> &ludb, std::vector<uint64_t> &simplexes,
    int32_t num_trials
)
{
    //FILE *f;
    std::ofstream f;

    f.open(filename, std::ofstream::out | std::ofstream::binary);
    if (f.fail())
    {
        return false;
    }

    for (int32_t i = 0; i < num_trials; i++)
    {
        f << ludb[i] << simplexes[i];
    }

    f.close();
    return true;
}

/*	This method evaluates the accuracy of the LUDB-heuristic method compared to
	the LUDB-exact algorithm.
	A number of runs are performed on a balanced-tree topology.
*/
bool LUDB_evaluation(Tandem &tandem, int32_t order, int32_t depth, int32_t max_good_combos)
{
    int32_t num_trials = 50;
    bool compute_exact = true;
    int32_t mgc_limit = 100;
    int32_t ludb_exact_preload;

    std::vector<uint64_t> simplexes_exact(num_trials);
    std::vector<double> ludb_exact(num_trials);
    std::vector<double> ludb_diff(num_trials);
    std::vector<uint64_t> ludb_simplexes(num_trials);
    std::vector<double> ludb_simplexes_percent(num_trials);
    double mean, delta;

    balanced_tree_descriptor bt_desc = {
        .order = order,
        .depth = depth,
        .node_latency_min = 0.0,
        .node_latency_max = 0.0,
        .node_rate_overprov_min = 1.01,
        .node_rate_overprov_max = 1.50,
        .flow_rate_min = 10.0,
        .flow_rate_max = 100.0,
        .flow_burst_min = 100.0,
        .flow_burst_max = 1000.0,
    };

    tandem.CreateTree(bt_desc, rng);
    if (tandem.NumFlows() > NMAX)
    {
        std::cout
            << "LUDB_Evaluation() error: too many t-nodes (" << tandem.NumFlows()
            << "), max " << NMAX
            << " allowed by simplex library."
            << std::endl;

        return false;
    }
    tandem.BuildNestingTree();

    std::cout << "LUDB Evaluation():  performing " << num_trials << " runs per simulation." << std::endl;
    std::cout << "max_good_combos=" << max_good_combos << std::endl;
    std::cout << "order=" << order << ", depth=" << depth << "  --> " << tandem.NumFlows() << " flows, max "
              << tandem.LUDB_TotalCombos() << " simplexes" << std::endl;
    std::cout << "flow_burst: " << bt_desc.flow_burst_min << " - " << bt_desc.flow_burst_max << std::endl;
    std::cout << "flow_rate: " << bt_desc.flow_rate_min << " - " << bt_desc.flow_rate_max << std::endl;
    std::cout << "node_latency: " << bt_desc.node_latency_min << " - " << bt_desc.node_latency_max << std::endl;
    std::cout
        << "node_rate_overprov: " << std::setprecision(2) << bt_desc.node_rate_overprov_min
        << " - " << std::setprecision(2) << bt_desc.node_rate_overprov_max
        << std::endl;
    std::cout << std::endl;

    std::stringstream fname_ss;
    fname_ss << "ludb-exact-o" << order << "-d" << depth << "-seed" << RNG_SEED << ".dat";
    auto fname = fname_ss.str();

    ludb_exact_preload = LUDB_LoadPrecomputedValues(fname, ludb_exact, simplexes_exact, num_trials);
    if (ludb_exact_preload > 0)
    {
        std::cout << "Loaded precomputed LUDB values for " << ludb_exact_preload << " runs (" << fname << ")"
                  << std::endl;
    }

    if (tandem.LUDB_TotalCombos() > 1000000)
    {
        compute_exact = false;
    }

    rng = RNG(RNG_SEED);
    if (compute_exact)
    {
        std::cout << "LUDB_exact\tSimplexes_Exact" << std::endl;
    }
    else
    {
        std::cout << "LUDB_MGC" << mgc_limit << "\tSimplexes_MGC" << mgc_limit << std::endl;
    }
    for (int32_t trial = 0; trial < num_trials; trial++)
    {
        tandem.CreateTree(bt_desc, rng);
        tandem.BuildNestingTree();
        if (simplexes_exact[trial] == 0)
        {
            if (compute_exact)
            {
                ludb_exact[trial] = tandem.LUDB_quiet(LUDB_MODE_EXACT, simplexes_exact[trial]);
            }
            else
            {
                ludb_exact[trial] = tandem.LUDB_quiet(mgc_limit, simplexes_exact[trial]);
            }
        }
        std::cout << std::setprecision(12) << ludb_exact[trial] << "\t" << simplexes_exact[trial] << std::endl;
    }
    std::cout << std::endl;
    if (ludb_exact_preload < num_trials)
    {
        if (LUDB_SavePrecomputedValues(fname, ludb_exact, simplexes_exact, num_trials))
        {
            std::cout << "Saved precomputed LUDB values for " << num_trials << " runs (" << fname << ")" << std::endl;
        }
    }

    std::cout
        << "Max_GC" << "\t" << "Mean_Err%" << "\t" << "Int_semiwidth" << "\t"
        << "Max_Err%" << "\t" << "Match%" << "\t\t" << "Mean_Simpl" << "\t" << "Int_semiwidth"
        << std::endl;
    for (int32_t k = 0; k < max_good_combos + 1; k++)
    {
        int32_t nc;
        double ludb_heuristic;
        double ludb_diff_max = 0.0;
        rng = RNG(RNG_SEED);
        uint64_t matches = 0;
        for (int32_t trial = 0; trial < num_trials; trial++)
        {
            tandem.CreateTree(bt_desc, rng);
            tandem.BuildNestingTree();
            if (compute_exact)
            {
                nc = (k < max_good_combos) ? k + 1 : LUDB_MODE_HEURISTIC;
            }
            else
            {
                nc = (k < max_good_combos) ? k + 1 : mgc_limit;
            }
            ludb_heuristic = tandem.LUDB_quiet(nc, ludb_simplexes[trial]);
            ludb_diff[trial] = (ludb_heuristic - ludb_exact[trial]) / ludb_exact[trial];
            ludb_diff_max = std::max(ludb_diff_max, ludb_diff[trial]);
            ludb_simplexes_percent[trial] = (double) ludb_simplexes[trial] / (double) simplexes_exact[trial];
            if (std::fabs(ludb_heuristic - ludb_exact[trial]) < TNODE_LUDB_EPSILON)
            {
                matches++;
            }
        }
        confidence95(ludb_diff, num_trials, mean, delta);

        std::cout
            << std::setw(3) << std::setfill('0') << nc << "\t"
            << std::setprecision(8) << mean << "\t"
            << std::setprecision(8) << delta << "\t"
            << std::setprecision(8) << ludb_diff_max << "\t"
            << std::setprecision(4) << (double) matches / (double) num_trials << "\t\t";

        confidence95(ludb_simplexes_percent, num_trials, mean, delta);
        std::cout
            << std::fixed << std::setprecision(8) << mean << "\t"
            << std::fixed << std::setprecision(8) << delta
            << std::endl;
    }
    return true;
}

/*	This method evaluates the lower-bound delay using either the exhaustive or the reduced method.
	A number of runs are performed on a balanced-tree topology.
*/
bool LowerBound_evaluation(Tandem &tandem, int32_t order, int32_t depth, bool experimental, float percentage)
{
    int32_t num_trials = 50;

    std::vector<double> lb_exact(num_trials);
    std::vector<uint64_t> lb_combos(num_trials);
    uint64_t lb_maxcombos;
    double mean, lb_max;
    uint64_t max_flows = 8 * sizeof(uint64_t);
    bool dump_config = false;

    balanced_tree_descriptor bt_desc;
    bt_desc.order = order;
    bt_desc.depth = depth;
    bt_desc.flow_burst_min = 100.0;
    bt_desc.flow_burst_max = 1000.0;
    bt_desc.flow_rate_min = 10.0;
    bt_desc.flow_rate_max = 100.0;
    bt_desc.node_rate_overprov_min = 1.01;
    bt_desc.node_rate_overprov_max = 1.50;
    bt_desc.node_latency_min = 0.1;
    bt_desc.node_latency_max = 0.1;

    tandem.CreateTree(bt_desc, rng);

    if (tandem.NumFlows() >= max_flows)
    {
        std::cout << "LowerBound_Evaluation() error: too many flows (" << tandem.NumFlows() << "), max " << max_flows
                  << " allowed." << std::endl;
        return false;
    }

    std::cout
        << "LowerBound Evaluation():  performing " << num_trials << " run(s) per simulation using "
        << (experimental ? "experimental" : "exhaustive") << " algorithm." << std::endl;
    std::cout << "order=" << order << ", depth=" << depth << "  --> " << tandem.NumFlows() << " flows" << std::endl;
    std::cout << "flow_burst: " << bt_desc.flow_burst_min << " - " << bt_desc.flow_burst_max << "" << std::endl;
    std::cout << "flow_rate: " << bt_desc.flow_rate_min << " - " << bt_desc.flow_rate_max << "" << std::endl;
    std::cout << "node_latency: " << bt_desc.node_latency_min << " - " << bt_desc.node_latency_max << "" << std::endl;
    std::cout << "node_rate_overprov: "
              << std::fixed << std::setprecision(2) << bt_desc.node_rate_overprov_min << " - "
              << std::fixed << std::setprecision(2) << bt_desc.node_rate_overprov_max
              << std::endl;
    std::cout << "random combinations percentage: " << std::fixed << std::setprecision(3) << percentage << " %"
              << std::endl;
    std::cout << std::endl;

    // determine the total number of combinations that should be evaluated
    lb_maxcombos = (uint64_t) pow(2.0, (double) tandem.NumFlows());
    mean = 0.0;
    lb_max = 0.0;

    // initialize random number generator
    rng = RNG(RNG_SEED);

    std::cout << "Trial\tLowerBound\tCombos\tMaxCombos" << std::endl;
    for (int32_t trial = 0; trial < num_trials; trial++)
    {
        tandem.CreateTree(bt_desc, rng);
        if (dump_config)
        {
            // save the generated tandem configuration for later debugging purposes

            std::stringstream fname_ss;
            fname_ss
                << "lb-o" << bt_desc.order
                << "-d" << bt_desc.depth
                << "_" << std::setw(2) << std::setfill('0') << trial
                << ".conf";
            tandem.Save(fname_ss.str());
        }
        // build nesting tree to initialize required data structures
        tandem.BuildNestingTree();
        if (experimental)
        {
            lb_exact[trial] = tandem.LowerBound_experimental_quiet(lb_combos[trial], percentage);
        }
        else
        {
            lb_exact[trial] = tandem.LowerBound(true, 0, 0, percentage);
            lb_combos[trial] = (uint64_t) pow(2.0, (double) tandem.NumFlows());
        }
        std::cout << trial << "\t" << lb_exact[trial] << "\t" << lb_combos[trial] << "\t" << lb_maxcombos << std::endl;

        // update some simple stats
        mean += lb_exact[trial];
        if (lb_exact[trial] - lb_max >= Curve::Epsilon)
        {
            lb_max = lb_exact[trial];
        }
    }
    mean = mean / num_trials;
    std::cout << std::endl;
    std::cout << "Mean=\t" << mean << "	Max=\t" << lb_max << std::endl;

    return true;
}

/*	This method saves a configuration file for a tandem representing
	a balanced tree with the given parameters
*/
bool MakeTreeConfig(char *filename, const int32_t order, const int32_t depth)
{
    Tandem t0;
    bool res;
    balanced_tree_descriptor bt_desc;
    bt_desc.order = order;
    bt_desc.depth = depth;
    bt_desc.flow_burst_min = 100.0;
    bt_desc.flow_burst_max = 1000.0; //150.0;
    bt_desc.flow_rate_min = 10.0;
    bt_desc.flow_rate_max = 100.0; //70.0;
    bt_desc.node_rate_overprov_min = 1.01;
    bt_desc.node_rate_overprov_max = 1.50; //1.1;
    bt_desc.node_latency_min = 0.0;
    bt_desc.node_latency_max = 0.0;

    rng = RNG(RNG_SEED);

    res = t0.CreateTree(bt_desc, rng);
    if (!res)
    {
        return false;
    }
    res = t0.Save(filename);
    return res;
}

bool MakeFullTandemConfig(char *filename, const int32_t num_nodes, const int32_t flows_percent)
{
    Tandem tandem;
    bool res;
    balanced_tree_descriptor bt_desc;
    bt_desc.order = 0;
    bt_desc.depth = flows_percent; // allocate nn% of all possible flows [N*(N+1)/2]
    bt_desc.flow_burst_min = 100.0;
    bt_desc.flow_burst_max = 1000.0;
    bt_desc.flow_rate_min = 10.0;
    bt_desc.flow_rate_max = 100.0;
    bt_desc.node_rate_overprov_min = 1.01;
    bt_desc.node_rate_overprov_max = 1.50;
    bt_desc.node_latency_min = 0.0;
    bt_desc.node_latency_max = 0.0;

    rng = RNG(RNG_SEED);

    bt_desc.order = num_nodes;
    int32_t retries = 10;
    do
    {
        tandem.CreateFullTandem(bt_desc, rng);
        retries--;
    }
    while (tandem.isNested() && retries > 0);
    if (retries == 0)
    {
        return false;
    }

    res = tandem.Save(filename);
    return res;
}

/* Compares two cutsets and returns the number of cuts before the two sets diverge.
 * Used to measure overlap between cutsets and reuse previous computations.
 */
uint64_t compare_cutsets(CutNodeList &c1, CutNodeList &c2)
{
    uint64_t l1 = c1.size();
    uint64_t l2 = c2.size();
    uint64_t len = std::min(l1, l2);
    for (uint64_t i = 0; i < len; i++)
    {
        if (c1[i] != c2[i])
        {
            return i;
        }
    }
    return len;
}

//#define LUDB_NON_NESTED_DEBUG 2

/*	Computes the delay bound for the non-nested tandem 'tandem'
 * 	'precomputed_cutsets' points to a pre-computed set of cuts; if NULL, an exhaustive computation will be performed
 *  'num_cuts' is initialized so that cut sets longer than (min_len + num_cuts) will be skipped
 */
/// Computes LUDB for a non-nested tandem using Multiple Sub-tandem Analysis (MSA)
/// For each set of cuts:
/// - compute the LUDB at each sub-tandem
/// - compute the delay for the whole tandem as the sum of these delays
/// \param tandem
/// \param mode todo: what are the options?
/// \param precomputed_cutsets If not nullptr, sets of cuts to be used. If nullptr, the sets are computed with tandem.compute_cutsets()
/// \param totsimpl todo: document
/// \param num_cuts In input, used to limit the sets of cuts considered to those not longer than longer than (min_len + num_cuts).
/// It is then updated with the resulting # of sets of cuts.
/// \return The computed LUDB
double LUDB_non_nested_MSA(
    Tandem &tandem,
    int32_t mode,
    std::vector<CutNodeList> *precomputed_cutsets,
    uint64_t &totsimpl,
    uint64_t &num_cuts
)
{
    double min_ludb = TNODE_LUDB_INFINITY;
    double ludb;
    totsimpl = 0;
    int32_t subtandems_total = 0, subtandems_reused = 0;
    std::vector<uint64_t> best_results;

    // compute all possible sets of cuts
    std::vector<CutNodeList> cutsets;
    cutsets = precomputed_cutsets == nullptr ? tandem.compute_cutsets() : *precomputed_cutsets;

    if (num_cuts != TANDEM_LUDB_MAXCUTS)
    {
        tandem.filter_cutsets(cutsets, num_cuts);
    }
    num_cuts = cutsets.size();

    // remember the last set of cuts to reuse LUDB computations
    std::vector<Tandem> last_cuts;
    // number of sub-tandem which could be reused from the previous step
    uint64_t cuts_reuse = 0;

    // sanity check of set of cuts found
    if (cutsets.empty())
    {
        std::cout << "Fatal error: the tandem is non-nested but no cuts have been found (bug?)" << std::endl;
        return DEBORAH_ERROR_SAFECHECK;
    }

    if (precomputed_cutsets == nullptr)
    {
        std::cout << "Number of sets of cuts to try: " << cutsets.size() << std::endl;
    }

    // for each set of cuts, compute end-to-end delay as the sum of intermediate LUDBs
    for (uint64_t cutset_idx = 0; cutset_idx < cutsets.size(); cutset_idx++)
    {
        ludb = 0.0;
        if (cutset_idx > 0)
        {
            // compare with the last cutset
            cuts_reuse = compare_cutsets(cutsets[cutset_idx], cutsets[cutset_idx - 1]);

            // safety check
            if (cuts_reuse == cutsets[cutset_idx].size())
            {
                std::cout << "*** CUTS_REUSE RANGE EXCEPTION ***" << std::endl;
                cuts_reuse = 0;
            }
        }
        else
        {
            cuts_reuse = 0;
        }

        // Cache current set of cut_tandem, to reuse computations in following iterations
        std::vector<Tandem> current_cuts;

        // Cut of the original tandem according to the i-th set of cut_tandem
        std::vector<Tandem> cut_tandem = tandem.cut(cutsets[cutset_idx]);

        // update stats about computation savings
        subtandems_total += cut_tandem.size();
        subtandems_reused += cuts_reuse;

#ifdef LUDB_NON_NESTED_DEBUG
        std::cout << std::endl;
        std::cout
            << "Cut set #" << std::setw(2) << setfill('0') << cutset_idx
            << " contains " << cut_tandem.size()
            << " sub-tandem(s)  (reusable=" << cuts_reuse << "):"
            << std::endl;

        for (uint64_t i = 0; i < cut_tandem.size(); i++)
            cut_tandem[i].Print();
#endif

        // for each sub-tandem, compute 1) the LUDB, and 2) the relevant output arrival curves
        for (uint64_t subtandem_idx = 0; subtandem_idx < cut_tandem.size(); subtandem_idx++)
        {
            double subtandem_ludb;
#if LUDB_NON_NESTED_DEBUG > 1
            std::cout << "********************************** sub-tandem " << subtandem_idx <<": **********************************" << std::endl;
#endif
            // check if we can reuse numeric results from the previous cut_tandem
            if (subtandem_idx < cuts_reuse)
            {
                // retrieve cached LUDB values
                Tandem subtandem = last_cuts[subtandem_idx];
                subtandem_ludb = subtandem.LUDB_getLast();
#ifdef LUDB_NON_NESTED_DEBUG
                std::cout
                    << "T" << subtandem_idx
                    << ".LUDB_cached = " << std::setprecision(3) << subtandem_ludb
                    << std::endl;
#endif
                ludb += subtandem_ludb;
                current_cuts.push_back(subtandem);

                // nothing else to do, continue with next sub-tandem
                continue;
            }
            else if ((subtandem_idx > 0) && (subtandem_idx == cuts_reuse))
            {
                // for this tandem we have to recompute LUDB using the proper O.A.Curves
                Tandem &subtandem = cut_tandem[subtandem_idx];
                Tandem &cached_subtandem = last_cuts[subtandem_idx];

                // copy output arrival curves from the previous cut_tandem
                for (uint64_t f = 0; f < subtandem.NumFlows(); f++)
                {
                    Flow *fp1 = subtandem.getFlow(f);
                    uint64_t fid = fp1->uid;
                    if (fid == FLOW_INVALID_ID)
                    {
                        continue;
                    }
                    Flow *fp2 = cached_subtandem.getFlowId(fid);
                    if (fp2 == nullptr)
                    {
                        continue;
                    }
                    fp1->burst = fp2->burst;
                    fp1->rate = fp2->rate;
                }
            }

            // compute LUDB for the tagged flow
            Tandem subtandem = cut_tandem[subtandem_idx];
            // save tandem
            current_cuts.push_back(subtandem);
            // aggregate flows to the maximum extent, including the tagged flow
            subtandem.aggregate_flows(true);
            subtandem.BuildNestingTree();
            //subtandem.Print(); subtandem.PrintNestingTree();
            uint64_t nsimpl = 0;
            // first try with the trivial cases
            subtandem_ludb = subtandem.LUDB_trivial();  //todo: look here
            // if not a supported case, perform normal LUDB computation
            if (subtandem_ludb < 0.0)
            {
                subtandem_ludb = subtandem.LUDB_quiet(mode, nsimpl); //todo: look here
            }
            current_cuts[subtandem_idx].LUDB_setLast(subtandem_ludb); //todo: look here

#ifdef LUDB_NON_NESTED_DEBUG
            std::cout << "T" << subtandem_idx << ".LUDB = " << std::setprecision(3) << subtandem_ludb << std::endl;
#endif
            ludb += subtandem_ludb; // add the sub-tandem's delay to the global LUDB
            totsimpl += nsimpl;

            // compute output arrival curves
            // except for the last sub-tandem
            if (subtandem_idx == cut_tandem.size() - 1)
                continue;

#if LUDB_NON_NESTED_DEBUG > 1
            std::cout << ">>> compute output arrival curves at exit of Tandem #" << subtandem_idx << std::endl;
#endif
            // recursive helper function that does the job
            cut_tandem[subtandem_idx].output_arrival_curves(cut_tandem[subtandem_idx + 1], mode);
#if LUDB_NON_NESTED_DEBUG > 1
            std::cout << "<<< done\n\n" << std::endl;
#endif
        }

        if (precomputed_cutsets == nullptr)
        {
            std::cout
                << "Cuts set " << std::setw(3) << std::setfill('0') << cutset_idx
                << "  ==>  LUDB= " << std::fixed << std::setprecision(4) << std::setfill('0') << ludb
                << " \t cuts_len= " << cutsets[cutset_idx].size()
                << " \tcuts={ ";

            for (auto cut : cutsets[cutset_idx])
            {
                std::cout << cut << " ";
            }

            std::cout << "}" << std::endl;
        }

        // keep track of the minimum LUDB found amongst all possible cut_tandem sets
        if (ludb <= min_ludb)
        {
            if (ludb < min_ludb)
            {
                best_results.clear();
                min_ludb = ludb;
            }
            best_results.push_back(cutset_idx);
        }

        // commit saved data for reuse
        last_cuts.clear();
        last_cuts = current_cuts;
    }

    if ((cutsets.size() > 1) && (precomputed_cutsets == nullptr))
    {
        double reusage = (double) subtandems_reused / (double) subtandems_total;
        std::cout
            << "Sub-tandem computations reusage ratio = " << std::setprecision(4) << reusage
            << "  (" << subtandems_reused << " / " << subtandems_total << ")"
            << std::endl;
    }

    // print the IDs of the cuts that correspond to the minimum LUDB found
    std::cout << "Best set(s) of cuts:  { ";
    for (auto best_result : best_results)
    {
        std::cout << std::setw(3) << std::setfill('0') << best_result << " ";
    }
    std::cout << "}" << std::endl;

    return min_ludb;
}

/// Computes LUDB for a non-nested tandem using Single Tandem Analysis (STA)
/// For each set of cuts:
/// - compute arrival curves at cuts' boundaries
/// - compute LUDB for the whole tandem (rather than the sum of per-cut LUDBs)
/// \param tandem
/// \param mode todo: what are the options?
/// \param precomputed_cutsets If not nullptr, sets of cuts to be used. If nullptr, the sets are computed with tandem.compute_cutsets()
/// \param totsimpl todo: document
/// \param num_cuts In input, used to limit the sets of cuts considered to those not longer than longer than (min_len + num_cuts).
/// It is then updated with the resulting # of sets of cuts.
/// \return The computed LUDB
double LUDB_non_nested_STA(
    Tandem &tandem,
    int32_t mode,
    std::vector<CutNodeList> *precomputed_cutsets,
    uint64_t &totsimpl,
    uint64_t &num_cuts
)
{
    double min_ludb = TNODE_LUDB_INFINITY;
    double ludb;

    totsimpl = 0;
    int32_t subtandems_total = 0, subtandems_reused = 0;
    std::vector<uint64_t> best_results;

    // compute all possible sets of cuts
    std::vector<CutNodeList> cutsets;
    cutsets = precomputed_cutsets == nullptr ? tandem.compute_cutsets() : *precomputed_cutsets;

    if (num_cuts != TANDEM_LUDB_MAXCUTS)
    {
        tandem.filter_cutsets(cutsets, num_cuts);
    }
    num_cuts = cutsets.size();

    // remember the last set of cuts to reuse LUDB computations
    std::vector<Tandem> last_cuts;
    // number of sub-tandem which could be reused from the previous step
    uint64_t cuts_reuse = 0;

    // sanity check of set of cuts found
    if (cutsets.empty())
    {
        std::cout << "Fatal error: the tandem is non-nested but no cuts have been found (bug?)" << std::endl;
        return DEBORAH_ERROR_SAFECHECK;
    }

    if (precomputed_cutsets == nullptr)
    {
        std::cout << "Number of sets of cuts to try: " << cutsets.size() << std::endl;
    }

    // for each set of cuts, compute the global end-to-end delay
    for (uint64_t cutset_idx = 0; cutset_idx < cutsets.size(); cutset_idx++)
    {
        ludb = 0.0;
        if (cutset_idx > 0)
        {
            // compare with the last cutset
            cuts_reuse = compare_cutsets(cutsets[cutset_idx], cutsets[cutset_idx - 1]);

            // safety check
            if (cuts_reuse == cutsets[cutset_idx].size())
            {
                std::cout << "*** CUTS_REUSE RANGE EXCEPTION ***" << std::endl;
                cuts_reuse = 0;
            }
        }
        else
        {
            cuts_reuse = 0;
        }

        // Cache current set of cut_tandem, to reuse computations in following iterations
        std::vector<Tandem> current_cuts;

        // cut the original tandem according to the i-th set of cut_tandem
        std::vector<Tandem> cut_tandem = tandem.cut(cutsets[cutset_idx]);

        // update stats about computation savings
        subtandems_total += cut_tandem.size();
        subtandems_reused += cuts_reuse;

#ifdef LUDB_NON_NESTED_DEBUG
        std::cout << "\nUsing alternative algorithm (single LUDB over joined tandem):" << std::endl;

        std::cout
            << "Cut set #" << std::setw(2) << setfill('0') << cutset_idx
            << " contains " << cut_tandem.size()
            << " sub-tandem(s)  (reusable=" << cuts_reuse << "):"
            << std::endl;

        for (uint64_t i = 0; i < cut_tandem.size(); i++)
            cut_tandem[i].Print();
#endif

        // for each sub-tandem, compute the relevant output arrival curves
        for (uint64_t subtandem_idx = 0; subtandem_idx < cut_tandem.size(); subtandem_idx++)
        {
#if LUDB_NON_NESTED_DEBUG > 1
            std::cout << "********************************** sub-tandem " << subtandem_idx << ": **********************************" << std::endl;
#endif
            // check if we can reuse numeric results from the previous cut
            if (subtandem_idx < cuts_reuse)
            {
                // retrieve cached LUDB values
                Tandem subtandem = last_cuts[subtandem_idx];
                current_cuts.push_back(subtandem);

                // nothing else to do, continue with next sub-tandem
                continue;
            }
            else if ((subtandem_idx > 0) && (subtandem_idx == cuts_reuse))
            {
                // for this tandem we have to recompute LUDB using the proper O.A.Curves
                Tandem &subtandem = cut_tandem[subtandem_idx];
                Tandem &cached_subtandem = last_cuts[subtandem_idx];

                // copy output arrival curves from the previous cut_tandem
                for (uint64_t f = 0; f < subtandem.NumFlows(); f++)
                {
                    Flow *fp1 = subtandem.getFlow(f);
                    uint64_t fid = fp1->uid;
                    if (fid == FLOW_INVALID_ID)
                        continue;

                    Flow *fp2 = cached_subtandem.getFlowId(fid);
                    if (fp2 == nullptr)
                        continue;

                    fp1->burst = fp2->burst;
                    fp1->rate = fp2->rate;
                }
            }

            Tandem subtandem = cut_tandem[subtandem_idx];
            current_cuts.push_back(subtandem);

            // compute new output arrival curves
            // except for the last sub-tandem
            if (subtandem_idx == cut_tandem.size() - 1)
                continue;

#if LUDB_NON_NESTED_DEBUG > 1
            std::cout << ">>> compute output arrival curves at exit of Tandem #" << subtandem_idx << std::endl;
#endif
            // recursive helper function that does the job
            cut_tandem[subtandem_idx].output_arrival_curves(cut_tandem[subtandem_idx + 1], mode);

#if LUDB_NON_NESTED_DEBUG > 1
            std::cout << "<<< done\n\n" << std::endl;
#endif
        }

        // now re-assemble sub-tandems into a global nested tandem
        Tandem reassembled_tandem;
        if (!reassembled_tandem.join(current_cuts))
        {
            std::cout << "Internal error: joined tandem is not nested, skipping!" << std::endl;
            last_cuts.clear();
            last_cuts = current_cuts;
            continue;
        }

        // compute LUDB for the global
        reassembled_tandem.aggregate_flows(false);
        reassembled_tandem.BuildNestingTree();

#if LUDB_NON_NESTED_DEBUG > 0
        std::cout << "Reassembled tandem:" << std::endl;
        reassembled_tandem.Print();
#endif

        uint64_t nsimpl = 0;
        ludb = reassembled_tandem.LUDB_trivial();
        if (ludb < 0.0)
            ludb = reassembled_tandem.LUDB_quiet(mode, nsimpl); //todo: look here

        // LUDB done
        if (precomputed_cutsets == nullptr)
        {
            std::cout
                << "Cuts set " << std::setw(3) << std::setfill('0') << cutset_idx
                << "  ==>  LUDB= " << std::fixed << std::setprecision(4) << std::setfill('0') << ludb
                << " \t cuts_len= " << cutsets[cutset_idx].size()
                << " \tcuts={ ";

            for (auto cut : cutsets[cutset_idx])
            {
                std::cout << cut << " ";
            }

            std::cout << "}" << std::endl;
        }

        // keep track of the minimum LUDB found amongst all possible cut sets
        if (ludb <= min_ludb)
        {
            if (ludb < min_ludb)
            {
                best_results.clear();
                min_ludb = ludb;
            }
            best_results.push_back(cutset_idx);
        }

        // commit saved data for reuse
        last_cuts.clear();
        last_cuts = current_cuts;
    }

    if ((cutsets.size() > 1) && (precomputed_cutsets == nullptr))
    {
        double reusage = (double) subtandems_reused / (double) subtandems_total;
        std::cout
            << "Sub-tandem computations reusage ratio = " << std::setprecision(4) << reusage
            << "  (" << subtandems_reused << " / " << subtandems_total << ")"
            << std::endl;
    }

    // print the IDs of the cuts that correspond to the minimum LUDB found
    std::cout << "Best set(s) of cuts:  { ";
    for (auto best_result : best_results)
    {
        std::cout << std::setw(3) << std::setfill('0') << best_result << " ";
    }
    std::cout << "}" << std::endl;

    return min_ludb;
}

/// Computes LUDB and CDFs for a non-nested tandem using Single Tandem Analysis (STA)
/// For each set of cuts:
/// - compute arrival curves at cuts' boundaries
/// - compute LUDB and output arrival curve for the whole tandem
/// \param tandem
/// \param mode todo: what are the options?
/// \param precomputed_cutsets If not nullptr, sets of cuts to be used. If nullptr, the sets are computed with tandem.compute_cutsets()
/// \param totsimpl todo: document
/// \param num_cuts In input, used to limit the sets of cuts considered to those not longer than longer than (min_len + num_cuts).
/// It is then updated with the resulting # of sets of cuts.
/// \return The computed LUDB, while output arrival curves are printed to standard output
double OutputCurve_non_nested_STA(
    Tandem &tandem,
    int32_t mode,
    std::vector<CutNodeList> *precomputed_cutsets,
    uint64_t &totsimpl,
    uint64_t &num_cuts
)
{
    double min_ludb = TNODE_LUDB_INFINITY;
    double ludb;

    totsimpl = 0;
    int32_t subtandems_total = 0, subtandems_reused = 0;

    std::vector<uint64_t> best_results;
    std::vector<OutputArrivalCurveParameters> best_CDFs;

    // compute all possible sets of cuts
    std::vector<CutNodeList> cutsets;
    cutsets = precomputed_cutsets == nullptr ? tandem.compute_cutsets() : *precomputed_cutsets;

    if (num_cuts != TANDEM_LUDB_MAXCUTS)
    {
        tandem.filter_cutsets(cutsets, num_cuts);
    }
    num_cuts = cutsets.size();

    // remember the last set of cuts to reuse LUDB computations
    std::vector<Tandem> last_cuts;
    // number of sub-tandem which could be reused from the previous step
    uint64_t cuts_reuse = 0;

    // sanity check of set of cuts found
    if (cutsets.empty())
    {
        std::cout << "Fatal error: the tandem is non-nested but no cuts have been found (bug?)" << std::endl;
        return DEBORAH_ERROR_SAFECHECK;
    }

    if (precomputed_cutsets == nullptr)
    {
        std::cout << "Number of sets of cuts to try: " << cutsets.size() << std::endl;
    }

    // for each set of cuts, compute the global end-to-end delay
    for (uint64_t cutset_idx = 0; cutset_idx < cutsets.size(); cutset_idx++)
    {
        ludb = 0.0;
        if (cutset_idx > 0)
        {
            // compare with the last cutset
            cuts_reuse = compare_cutsets(cutsets[cutset_idx], cutsets[cutset_idx - 1]);

            // safety check
            if (cuts_reuse == cutsets[cutset_idx].size())
            {
                std::cout << "*** CUTS_REUSE RANGE EXCEPTION ***" << std::endl;
                cuts_reuse = 0;
            }
        }
        else
        {
            cuts_reuse = 0;
        }

        // Cache current set of cut_tandem, to reuse computations in following iterations
        std::vector<Tandem> current_cuts;

        // cut the original tandem according to the i-th set of cut_tandem
        std::vector<Tandem> cut_tandem = tandem.cut(cutsets[cutset_idx]);

        // update stats about computation savings
        subtandems_total += cut_tandem.size();
        subtandems_reused += cuts_reuse;

        // for each sub-tandem, compute the relevant output arrival curves
        for (uint64_t subtandem_idx = 0; subtandem_idx < cut_tandem.size(); subtandem_idx++)
        {
            // check if we can reuse numeric results from the previous cut
            if (subtandem_idx < cuts_reuse)
            {
                // retrieve cached LUDB values
                Tandem subtandem = last_cuts[subtandem_idx];
                current_cuts.push_back(subtandem);

                // nothing else to do, continue with next sub-tandem
                continue;
            }
            else if ((subtandem_idx > 0) && (subtandem_idx == cuts_reuse))
            {
                // for this tandem we have to recompute LUDB using the proper O.A.Curves
                Tandem &subtandem = cut_tandem[subtandem_idx];
                Tandem &cached_subtandem = last_cuts[subtandem_idx];

                // copy output arrival curves from the previous cut
                for (uint64_t f = 0; f < subtandem.NumFlows(); f++)
                {
                    Flow *fp1 = subtandem.getFlow(f);
                    uint64_t fid = fp1->uid;
                    if (fid == FLOW_INVALID_ID)
                        continue;

                    Flow *fp2 = cached_subtandem.getFlowId(fid);
                    if (fp2 == nullptr)
                        continue;

                    fp1->burst = fp2->burst;
                    fp1->rate = fp2->rate;
                }
            }

            Tandem subtandem = cut_tandem[subtandem_idx];
            current_cuts.push_back(subtandem);

            // compute new output arrival curves
            // except for the last sub-tandem
            if (subtandem_idx == cut_tandem.size() - 1)
                continue;

            // recursive helper function that does the job
            cut_tandem[subtandem_idx].output_arrival_curves(cut_tandem[subtandem_idx + 1], mode);
        }

        // now re-assemble sub-tandems into a global nested tandem
        Tandem reassembled_tandem;
        if (!reassembled_tandem.join(current_cuts))
        {
            std::cout << "Internal error: joined tandem is not nested, skipping!" << std::endl;
            last_cuts.clear();
            last_cuts = current_cuts;
            continue;
        }

        // compute LUDB for the global
        reassembled_tandem.aggregate_flows(false);
        reassembled_tandem.BuildNestingTree();

        uint64_t nsimpl = 0;
        ludb = reassembled_tandem.LUDB_trivial();
        if (ludb < 0.0)
        {
            ludb = reassembled_tandem.LUDB_quiet(mode, nsimpl);
        }

        // LUDB done
        if (precomputed_cutsets == nullptr)
        {
            std::cout
                << "Cuts set " << std::setw(3) << std::setfill('0') << cutset_idx
                << "  ==>  LUDB= " << std::fixed << std::setprecision(4) << std::setfill('0') << ludb
                << " \t cuts_len= " << cutsets[cutset_idx].size()
                << " \tcuts={ ";

            for (auto cut : cutsets[cutset_idx])
            {
                std::cout << cut << " ";
            }

            std::cout << "}" << std::endl;
        }

        // keep track of the minimum LUDB found amongst all possible cut sets
        if (ludb <= min_ludb)
        {
            OutputArrivalCurveParameters foiCdfParameters =
                reassembled_tandem.GetFoiOuputArrivalCurve(mode);

            if (ludb < min_ludb)
            {
                best_results.clear();
                best_CDFs.clear();
                min_ludb = ludb;
            }
            best_results.push_back(cutset_idx);
            best_CDFs.push_back(foiCdfParameters);
        }

        // commit saved data for reuse
        last_cuts.clear();
        last_cuts = current_cuts;
    }

    if ((cutsets.size() > 1) && (precomputed_cutsets == nullptr))
    {
        double reusage = (double) subtandems_reused / (double) subtandems_total;
        std::cout
            << "Sub-tandem computations reusage ratio = " << std::setprecision(4) << reusage
            << "  (" << subtandems_reused << " / " << subtandems_total << ")"
            << std::endl;
    }

    // print the IDs of the cuts that correspond to the minimum LUDB found
    std::cout << "Best set(s) of cuts:  { ";
    for (auto best_result : best_results)
    {
        std::cout << std::setw(3) << std::setfill('0') << best_result << " ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Output arrival curve parameters: " << endl <<
              "[" << endl;

    auto cdfs = best_CDFs.size();
    for (uint64_t c = 0; c < cdfs; c++)
    {
        OutputArrivalCurveParameters param = best_CDFs[c];
        std::cout
            << "\t{ \"burst\" : " << std::fixed << std::setprecision(6) << param.burst
            << " , \"rate\" : " << std::fixed << std::setprecision(6) << param.rate
            << " }";

        if (cdfs > 1 && c < cdfs - 1)
        {
            std::cout << ',' << endl;
        }
        else
        {
            std::cout << endl;
        }
    }

    std::cout << "]" << endl;

    return min_ludb;
}

/// Computes an equivalent service curve for a non-nested tandem.
/// The tandem cut is selected using Single Tandem Analysis (STA),
/// but for the given cut the eq. curve is selected randomly
/// \param tandem
/// \param mode todo: what are the options?
/// \param precomputed_cutsets If not nullptr, sets of cuts to be used. If nullptr, the sets are computed with tandem.compute_cutsets()
/// \param totsimpl todo: document
/// \param num_cuts In input, used to limit the sets of cuts considered to those not longer than longer than (min_len + num_cuts).
/// It is then updated with the resulting # of sets of cuts.
/// \return The computed LUDB, while output arrival curves are printed to standard output
double EquivalentServiceCurve_non_nested_STA(
    Tandem &tandem,
    int32_t mode,
    std::vector<CutNodeList> *precomputed_cutsets,
    uint64_t &totsimpl,
    uint64_t &num_cuts
)
{
    double min_ludb = TNODE_LUDB_INFINITY;
    double ludb;

    totsimpl = 0;
    int32_t subtandems_total = 0, subtandems_reused = 0;

    std::vector<uint64_t> best_results;
    std::vector<PseudoAffine> best_EQs;

    // compute all possible sets of cuts
    std::vector<CutNodeList> cutsets;
    cutsets = precomputed_cutsets == nullptr ? tandem.compute_cutsets() : *precomputed_cutsets;

    if (num_cuts != TANDEM_LUDB_MAXCUTS)
    {
        tandem.filter_cutsets(cutsets, num_cuts);
    }
    num_cuts = cutsets.size();

    // remember the last set of cuts to reuse LUDB computations
    std::vector<Tandem> last_cuts;
    // number of sub-tandem which could be reused from the previous step
    uint64_t cuts_reuse = 0;

    // sanity check of set of cuts found
    if (cutsets.empty())
    {
        std::cout << "Fatal error: the tandem is non-nested but no cuts have been found (bug?)" << std::endl;
        return DEBORAH_ERROR_SAFECHECK;
    }

    if (precomputed_cutsets == nullptr)
    {
        std::cout << "Number of sets of cuts to try: " << cutsets.size() << std::endl;
    }

    // for each set of cuts, compute the global end-to-end delay
    for (uint64_t cutset_idx = 0; cutset_idx < cutsets.size(); cutset_idx++)
    {
        ludb = 0.0;
        if (cutset_idx > 0)
        {
            // compare with the last cutset
            cuts_reuse = compare_cutsets(cutsets[cutset_idx], cutsets[cutset_idx - 1]);

            // safety check
            if (cuts_reuse == cutsets[cutset_idx].size())
            {
                std::cout << "*** CUTS_REUSE RANGE EXCEPTION ***" << std::endl;
                cuts_reuse = 0;
            }
        }
        else
        {
            cuts_reuse = 0;
        }

        // Cache current set of cut_tandem, to reuse computations in following iterations
        std::vector<Tandem> current_cuts;

        // cut the original tandem according to the i-th set of cut_tandem
        std::vector<Tandem> cut_tandem = tandem.cut(cutsets[cutset_idx]);

        // update stats about computation savings
        subtandems_total += cut_tandem.size();
        subtandems_reused += cuts_reuse;

        // for each sub-tandem, compute the relevant output arrival curves
        for (uint64_t subtandem_idx = 0; subtandem_idx < cut_tandem.size(); subtandem_idx++)
        {
            // check if we can reuse numeric results from the previous cut
            if (subtandem_idx < cuts_reuse)
            {
                // retrieve cached LUDB values
                Tandem subtandem = last_cuts[subtandem_idx];
                current_cuts.push_back(subtandem);

                // nothing else to do, continue with next sub-tandem
                continue;
            }
            else if ((subtandem_idx > 0) && (subtandem_idx == cuts_reuse))
            {
                // for this tandem we have to recompute LUDB using the proper O.A.Curves
                Tandem &subtandem = cut_tandem[subtandem_idx];
                Tandem &cached_subtandem = last_cuts[subtandem_idx];

                // copy output arrival curves from the previous cut
                for (uint64_t f = 0; f < subtandem.NumFlows(); f++)
                {
                    Flow *fp1 = subtandem.getFlow(f);
                    uint64_t fid = fp1->uid;
                    if (fid == FLOW_INVALID_ID)
                        continue;

                    Flow *fp2 = cached_subtandem.getFlowId(fid);
                    if (fp2 == nullptr)
                        continue;

                    fp1->burst = fp2->burst;
                    fp1->rate = fp2->rate;
                }
            }

            Tandem subtandem = cut_tandem[subtandem_idx];
            current_cuts.push_back(subtandem);

            // compute new output arrival curves
            // except for the last sub-tandem
            if (subtandem_idx == cut_tandem.size() - 1)
                continue;

            // recursive helper function that does the job
            cut_tandem[subtandem_idx].output_arrival_curves(cut_tandem[subtandem_idx + 1], mode);
        }

        // now re-assemble sub-tandems into a global nested tandem
        Tandem reassembled_tandem;
        if (!reassembled_tandem.join(current_cuts))
        {
            std::cout << "Internal error: joined tandem is not nested, skipping!" << std::endl;
            last_cuts.clear();
            last_cuts = current_cuts;
            continue;
        }

        // compute LUDB for the global
        reassembled_tandem.aggregate_flows(false);
        reassembled_tandem.BuildNestingTree();

        uint64_t nsimpl = 0;
        ludb = reassembled_tandem.LUDB_trivial();
        if (ludb < 0.0)
        {
            ludb = reassembled_tandem.LUDB_quiet(mode, nsimpl);
        }

        // LUDB done
        if (precomputed_cutsets == nullptr)
        {
            std::cout
                << "Cuts set " << std::setw(3) << std::setfill('0') << cutset_idx
                << "  ==>  LUDB= " << std::fixed << std::setprecision(4) << std::setfill('0') << ludb
                << " \t cuts_len= " << cutsets[cutset_idx].size()
                << " \tcuts={ ";

            for (auto cut : cutsets[cutset_idx])
            {
                std::cout << cut << " ";
            }

            std::cout << "}" << std::endl;
        }

        // keep track of the minimum LUDB found amongst all possible cut sets
        if (ludb <= min_ludb)
        {
            auto eq_curve =
                reassembled_tandem.GetEquivalentServiceCurve(num_cuts);

            if (ludb < min_ludb)
            {
                best_results.clear();
                best_EQs.clear();
                min_ludb = ludb;
            }
            best_results.push_back(cutset_idx);
            best_EQs.push_back(eq_curve);
        }

        // commit saved data for reuse
        last_cuts.clear();
        last_cuts = current_cuts;
    }

    if ((cutsets.size() > 1) && (precomputed_cutsets == nullptr))
    {
        double reusage = (double) subtandems_reused / (double) subtandems_total;
        std::cout
            << "Sub-tandem computations reusage ratio = " << std::setprecision(4) << reusage
            << "  (" << subtandems_reused << " / " << subtandems_total << ")"
            << std::endl;
    }

    // print the IDs of the cuts that correspond to the minimum LUDB found
    std::cout << "Best set(s) of cuts:  { ";
    for (auto best_result : best_results)
    {
        std::cout << std::setw(3) << std::setfill('0') << best_result << " ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Equivalent service curves: " << endl <<
              "[" << endl;

    auto eqs = best_EQs.size();
    for (uint64_t eq_idx = 0; eq_idx < eqs; eq_idx++)
    {
        PseudoAffine eq = best_EQs[eq_idx];
        std::cout << eq.PrintJson();

        if (eqs > 1 && eq_idx < eqs - 1)
        {
            std::cout << ", " << endl;
        }
        else
        {
            std::cout << endl;
        }
    }

    std::cout << "]" << endl;

    return min_ludb;
}

/*	Run a series of simulations to assess the performance of the LUDB analysis
 * 	in non-nested tandems.
 * 	'nnodes' is the length of the tandem (number of nodes)
 * 	'flows_percent' is the percentage (1 - 100) of how many flows will be allocated
 * 	'mode' specifies the LUDB computation mode (exact or heuristic)
 * 	'maxcuts' limits the length of the analyzed cuts to 'min_length+maxcuts'
 */
bool
LUDB_non_nested_evaluation(int32_t nnodes, int32_t flows_percent, int32_t mode, uint64_t maxcuts, bool algo_classic)
{
    balanced_tree_descriptor bt_desc;
    bt_desc.order = 0;
    bt_desc.depth = flows_percent; // allocate nn% of all possible flows [N*(N+1)/2]
    bt_desc.flow_burst_min = 100.0;
    bt_desc.flow_burst_max = 1000.0;
    bt_desc.flow_rate_min = 10.0;
    bt_desc.flow_rate_max = 100.0;
    bt_desc.node_rate_overprov_min = 1.01;
    bt_desc.node_rate_overprov_max = 1.50;
    bt_desc.node_latency_min = 0.0;
    bt_desc.node_latency_max = 0.0;

    //BasicTimer t0;
    NanoTimer t0;

    std::cout
        << "LUDB_non_nested_evaluation(): tandem length = " << nnodes
        << " nodes, " << bt_desc.depth
        << "% flow coverage, algorithm = " << (algo_classic ? "MTA" : "STA")
        << std::endl;

    std::cout << "                              LUDB computation algorithm = ";

    if (mode == LUDB_MODE_EXACT)
    {
        std::cout << "exact" << std::endl;
    }
    else if (mode == LUDB_MODE_HEURISTIC)
    {
        std::cout << "heuristic, all results" << std::endl;
    }
    else
    {
        std::cout << "heuristic, best " << mode << " results" << std::endl;
    }

    //fprintf(stderr, "num_nodes,perc_flows,num_flows,num_cutsets,time,ludb,num_simplexes\n\n");
    rng = RNG(RNG_SEED);
    Tandem tandem;
    uint64_t nsimpl = 0;
    uint64_t num_cuts = maxcuts;
    double ludb = 0.0;
    double t_compcuts, t_ludb;

    for (int32_t n = nnodes; n <= nnodes; n++)
    {
        bt_desc.order = n;
        int32_t retries = 10;
        do
        {
            tandem.CreateFullTandem(bt_desc, rng);
            retries--;
        }
        while (tandem.isNested() && retries > 0);
        uint64_t nflows = tandem.NumFlows();
        //tandem.getFlow(1)->burst *= 200;
        tandem.Print();

#ifdef LUDB_NNESTED_EVAL_BENCHMARK
        t0.mark(true);
        std::vector<CutNodeList> cuts = tandem.compute_cutsets();
        t_compcuts = t0.mark(false);
        num_cuts = _MIN(cuts.size(), num_cuts);

        t0.mark(true);
        // Comment next if() to skip LUDB and measure cuts computation time only
        //if(algo_classic) ludb = LUDB_non_nested_MSA(tandem, mode, &cuts, nsimpl, num_cuts);
        //else ludb = LUDB_non_nested_STA(tandem, mode, &cuts, nsimpl, num_cuts);
        ludb = 1.0;
        t0.mark(false);

        for (int32_t i = 0; i < num_cuts; i++)
        {
            std::cout << "set #%02i: { ", i);
            for (int32_t j = 0; j < cuts[i].size(); j++)
                std::cout << "%i ", cuts[i][j]);
            std::cout << "}\t len=%i\n", cuts[i].size());
        }
        //fprintf(stderr, "%i,%i,%i,%i,%.3lf,%.3lf,%lf,%Lu\n", n, flows_percent, nflows, num_cuts, t_compcuts,0.0,0.0,0);
        //continue;
#else
        t_compcuts = -1.0;
        t0.mark(true);
        if (algo_classic)
        {
            ludb = LUDB_non_nested_MSA(tandem, mode, nullptr, nsimpl, num_cuts);
        }
        else
        {
            ludb = LUDB_non_nested_STA(tandem, mode, nullptr, nsimpl, num_cuts);
        }
        t_ludb = t0.mark(false);
        std::vector<CutNodeList> cuts;
#endif

        // compare with per-node analysis
        cuts.clear();
        CutNodeList allnodes;
        uint64_t nsimpl2 = 0;
        uint64_t ncuts2 = TANDEM_LUDB_MAXCUTS;
        for (uint64_t i = 1; i <= tandem.NumNodes(); i++)
        {
            allnodes.push_back(i);
        }
        cuts.push_back(allnodes);
        double ludb_pernode = LUDB_non_nested_MSA(tandem, mode, &cuts, nsimpl2, ncuts2);
        double ratio = ludb_pernode / ludb;
        std::cout
            << "LUDB= " << ludb
            << " \tpnLUDB= " << ludb_pernode
            << " \tpnLUDB/LUDB= " << ratio
            << std::endl;

        std::cerr
            << n << ","
            << flows_percent << ","
            << nflows << ","
            << num_cuts << ","
            << std::setprecision(3) << t_compcuts << ","
            << std::setprecision(3) << t_ludb << ","
            << ratio << ","
            << nsimpl
            << std::endl;
    }
    return true;
}

/*	Compute a delay bound for tandem 't1' using per-node analysis
 */
double DelayBound_per_node(Tandem &t1, int32_t ludb_mode)
{
    std::vector<CutNodeList> cuts;
    CutNodeList allnodes;
    uint64_t nsimpl = 0;
    uint64_t ncuts = TANDEM_LUDB_MAXCUTS;
    for (uint64_t i = 1; i <= t1.NumNodes(); i++)
    {
        allnodes.push_back(i);
    }
    cuts.push_back(allnodes);
    double ludb_pernode = LUDB_non_nested_MSA(t1, ludb_mode, &cuts, nsimpl, ncuts);
    return ludb_pernode;
}
