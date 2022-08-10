//
// Created by raffa on 17/03/2019.
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <cinttypes>
#include <sstream>

#include "deborah-core.h"
#include "lib/tandem.h"
#include "lib/curve.h"
#include "lib/polynomial.h"
#include "lib/rng.h"
#include "lib/simplex/simplex.h"
#include "lib/timing.h"

using namespace std;
using namespace deborah;

/* ******************************************************************************************************* */
/* **********************************   DEBORAH application main code   ********************************** */
/* ******************************************************************************************************* */

// todo: update interface, clearer switch between MSA and STA

/*	Print command line syntax
*/
int32_t print_help()
{
    std::cout
        << "Usage: deborah <config file> [options]\n" << std::endl
        << "A tiny Network Calculus tool to compute LUDB and Lower-bounds in tandems.\n" << std::endl
        << "       Options:   --ludb, --noludb     :  enable/disable computation of LUDB" << std::endl
        << "                  --foi-output-curve   :  compute foi's output arrival curve" << std::endl
        << "                  --eq-service-curve   :  compute equivalent service curve for the foi" << std::endl
        << "                  --ludb-heuristic <N> :  compute LUDB using the heuristic algorithm (max combos returned = N)"
        << std::endl
        << "                  --ludb-evaluate <tree_order> <tree_depth> <max_good_combos>" << std::endl
        << "                                       :  run LUDB computations on a balanced tree tandem with random provisioning"
        << std::endl
        << "                  --ludb-nnested-evaluate <N> <F>" << std::endl
        << "                                       :  run LUDB on a non-nested tandem with N nodes and F%% flows out of the max possible"
        << std::endl
        << "                  --ludb-nnested-sta   :  use STA algorithm to compute LUDB for non-nested tandems"
        << std::endl
        << "                  --ludb-nnested-single:  same as --ludb-nnested-sta" << std::endl
        << "                  --ludb-cuts-len  <N> :  limit max cuts length to <min_len+N> nodes" << std::endl
        << "                  --per-node           :  compute a delay bound using per-node analysis" << std::endl
        << "                  --lb, --nolb         :  enable/disable computation of Lower-bound" << std::endl
        << "                  --lb-experimental    :  turn on experimental features in Lower-bound computation (UNSTABLE)"
        << std::endl
        << "                  --lb-evaluate <tree_order> <tree_depth>" << std::endl
        << "                                       :  run Lower-bound computations on a balanced tree tandem with random provisioning"
        << std::endl
        << "                  --lb-random-combo <P>:  percentage of combos to be computed for Lower-bound (default 100.0)"
        << std::endl
        << "                  --tagged <N>         :  assume flow #N as tagged flow (by default it's the longest one)"
        << std::endl
        << "                  --scale-rates <F> <N>:  scale provisioned flow and node rates by F and N times, respectively"
        << std::endl
        << "                  --gen-tree <tree_order> <tree_depth> <filename> :  generate balanced-tree scenario config file"
        << std::endl
        << "                                          if <tree_order> is 1, then a sink-tree topology is generated."
        << std::endl
        << "                  --gen-nnested <nodes> <flows percentage> <filename> :  generate non-nested tandem config file"
        << std::endl
        << "                  --rng-seed <N>       :  set seed for internal random number generator." << std::endl
        << "                  --det-output         :  output is deterministic, for testing. Hides header and runtime."
        << std::endl
        << std::endl;
    return 1;
}

/*	Dumb compare function between LUDB and LowerBound
*/
double compare_results(double ludb, double lowerbound, double ludb_pn, bool c_ludb, bool c_lowerbound)
{
    double delta;
    double rob;
    bool c_pernode = (ludb_pn > 0.0);

    if (!c_ludb && !c_lowerbound && !c_pernode)
    {
        return 0.0;
    }

    std::cout << "DEBORAH results:\n-------------------------" << std::endl;
    delta = ludb - lowerbound;
    if (c_ludb)
    {
        std::cout << "      LUDB = " << std::setprecision(6) << ludb << std::endl;
    }
    if (c_lowerbound)
    {
        std::cout << "LowerBound = " << std::setprecision(6) << lowerbound << std::endl;
    }
    if (c_pernode)
    {
        std::cout << "   PN-LUDB = " << ludb_pn << std::endl;
    }
    if (c_ludb && c_pernode)
    {
        std::cout << "  PN-Ratio = " << ludb_pn / ludb << std::endl;
    }
    if (!c_ludb || !c_lowerbound)
    {
        return 0.0;
    }
    std::cout << "     delta = " << delta << std::endl;
    rob = 1.0 - (lowerbound / ludb);
    std::cout << "       ROB = " << rob << std::endl;
    if ((delta >= 0.0) && (delta < Curve::Epsilon))
    {
        std::cout << "The bounds computed are tight ;-)" << std::endl;
    }
    else if (delta < 0.0)
    {
        std::cout << "Ooops, LUDB is lower than LowerBound :-[" << std::endl;
    }
    return delta;
}

struct Arguments
{
    std::string conf_file = "";
    uint64_t TaggedFlow = -1;
    bool compute_ludb = false;
    bool compute_output_curve = false;
    bool compute_eq_service_curve = false;
    bool compute_ludb_heuristic = false;
    bool compute_ludb_nnested_alt = false;
    int32_t ludb_max_good_combos = LUDB_MODE_HEURISTIC;
    int32_t ludb_tree_order = 2;
    int32_t ludb_tree_depth = 3;
    bool ludb_evaluation = false;
    bool ludb_nnested_evaluation = false;
    bool compute_ludb_pernode = false;
    bool lowerbound_evaluation = false;
    bool compute_lowerbound = false;
    bool lowerbound_experimental = false;
    float lowerbound_percentage = 100.0;
    uint64_t ludb_maxcuts = TANDEM_LUDB_MAXCUTS;
    double flowrate_mult = 1.0, noderate_mult = 1.0;
    uint64_t combo_range_start = 0, combo_range_len = 0;
    bool deterministicOutput = false;
};

Arguments parseConfiguration(int32_t argc, char *argv[])
{
    Arguments args;

    // parse command line arguments
    for (int32_t i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--noludb"))
        {
            args.compute_ludb = false;
        }
        else if (!strcmp(argv[i], "--nolb"))
        {
            args.compute_lowerbound = false;
        }
        else if (!strcmp(argv[i], "--lb"))
        {
            args.compute_lowerbound = true;
        }
        else if (!strcmp(argv[i], "--foi-output-curve"))
        {
            args.compute_output_curve = true;
        }
        else if (!strcmp(argv[i], "--eq-service-curve"))
        {
            args.compute_eq_service_curve = true;
        }
        else if (!strcmp(argv[i], "--lb-experimental"))
        {
            args.lowerbound_experimental = true;
        }
        else if (!strcmp(argv[i], "--ludb-evaluate"))
        {
            args.ludb_evaluation = true;
            if (i < argc - 3)
            {
                std::istringstream(argv[i + 1]) >> args.ludb_tree_order;
                std::istringstream(argv[i + 2]) >> args.ludb_tree_depth;
                std::istringstream(argv[i + 3]) >> args.ludb_max_good_combos;
                i += 3;
            }
        }
        else if (!strcmp(argv[i], "--ludb-nnested-evaluate"))
        {
            args.ludb_nnested_evaluation = true;
            args.ludb_max_good_combos = LUDB_MODE_EXACT;
            args.ludb_tree_depth = 100; // default: allocate all flows (100%)
            if (i < argc - 2)
            {
                std::istringstream(argv[i + 1]) >> args.ludb_tree_order;
                std::istringstream(argv[i + 2]) >> args.ludb_tree_depth;
                i += 2;
            }
            if (args.ludb_tree_order < 3)
            {
                args.ludb_tree_order = 3;
            }
        }
        else if (!strcmp(argv[i], "--lb-evaluate"))
        {
            args.lowerbound_evaluation = true;
            if (i < argc - 2)
            {
                std::istringstream(argv[i + 1]) >> args.ludb_tree_order;
                std::istringstream(argv[i + 2]) >> args.ludb_tree_depth;
                i += 2;
            }
        }
        else if ((!strcmp(argv[i], "--lb-random-combo")) && (i < argc - 1))
        {
#ifdef TANDEM_LB_RANDOM_COMBOS
            std::istringstream(argv[i + 1]) >> args.lowerbound_percentage;
#endif
            i++;
        }
        else if (!strcmp(argv[i], "--gen-tree"))
        {
            if (i < argc - 3)
            {
                std::istringstream(argv[i + 1]) >> args.ludb_tree_order;
                std::istringstream(argv[i + 2]) >> args.ludb_tree_depth;
                MakeTreeConfig(
                    argv[i + 3],
                    args.ludb_tree_order,
                    args.ludb_tree_depth
                );
            }
            std::cout
                << "Balanced tree (order=" << args.ludb_tree_order
                << ", depth=" << args.ludb_tree_depth
                << ") written to \"" << argv[i + 3] << "\""
                << std::endl;
            exit(0);
        }
        else if (!strcmp(argv[i], "--gen-nnested"))
        {
            if (i < argc - 3)
            {
                std::istringstream(argv[i + 1]) >> args.ludb_tree_order;
                std::istringstream(argv[i + 2]) >> args.ludb_tree_depth;
                MakeFullTandemConfig(
                    argv[i + 3],
                    args.ludb_tree_order,
                    args.ludb_tree_depth
                );
            }

            std::cout
                << "Non-nested tandem (nodes=" << args.ludb_tree_order <<
                ", flows percentage=" << args.ludb_tree_depth <<
                ") written to \"" << argv[i + 3] << "\""
                << std::endl;

            exit(0);
        }
        else if (!strcmp(argv[i], "--combo-range"))
        {
            if (i < argc - 2)
            {
                std::istringstream(argv[i + 1]) >> args.combo_range_start;
                std::istringstream(argv[i + 2]) >> args.combo_range_len;
                i += 2;
            }
        }
        else if (!strcmp(argv[i], "--ludb"))
        {
            args.compute_ludb = true;
            args.ludb_max_good_combos = LUDB_MODE_EXACT;
        }
        else if (!strcmp(argv[i], "--ludb-heuristic"))
        {
            args.compute_ludb_heuristic = args.compute_ludb = true;
            args.ludb_max_good_combos = -1;
            if (i < argc - 1)
            {
                auto ss = std::istringstream(argv[i + 1]);
                ss >> args.ludb_max_good_combos;
                if (ss.good())
                {
                    i++;
                }
            }
            if (args.ludb_max_good_combos <= 0)
            {
                args.ludb_max_good_combos = LUDB_MODE_HEURISTIC;
            }
        }
        else if (!strcmp(argv[i], "--ludb-cuts-len"))
        {
            if (i < argc - 1)
            {
                std::istringstream(argv[i + 1]) >> args.ludb_maxcuts;
                i++;
            }
        }
        else if (!strcmp(argv[i], "--ludb-nnested-single") || !strcmp(argv[i], "--ludb-nnested-sta"))
        {
            args.compute_ludb_nnested_alt = true;
        }
        else if (!strcmp(argv[i], "--per-node"))
        {
            args.compute_ludb_pernode = true;
            args.ludb_max_good_combos = LUDB_MODE_EXACT;
        }
        else if ((!strcmp(argv[i], "--tagged")) && (i < argc - 1))
        {
            std::istringstream(argv[i + 1]) >> args.TaggedFlow;
            i++;
        }
        else if ((!strcmp(argv[i], "--threads")) && (i < argc - 1)) // todo: wtf, pls do this
        {
            //sscanf_s(argv[i+1],"%i",&NumThreads);
            i++;
        }
        else if ((!strcmp(argv[i], "--rng-seed")) && (i < argc - 1))
        {
            std::istringstream(argv[i + 1]) >> RNG_SEED;
            i++;
        }
        else if (!strcmp(argv[i], "--scale-rates"))
        {
            if (i < argc - 2)
            {
                std::istringstream(argv[i + 1]) >> args.flowrate_mult;
                std::istringstream(argv[i + 2]) >> args.noderate_mult;
                i += 2;
                if (args.flowrate_mult <= 0.0)
                {
                    args.flowrate_mult = 1.0;
                }
                if (args.noderate_mult <= 0.0)
                {
                    args.noderate_mult = 1.0;
                }
            }
        }
        else if (!strcmp(argv[i], "--det-output"))
        {
            args.deterministicOutput = true;
        }
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
        {
            exit(print_help());
        }
            // last cases: either an unrecognized option, or the scenario file has not been specified
        else if (argv[i][0] == '-')
        {
            std::cout << "Error: unrecognized option \"" << argv[i] << "\"\n";
            std::cout << "       Please run \"" << argv[0] << " --help\" for a list of valid command line options.\n";
            exit(DEBORAH_ERROR_CLI);
        }
        else
        {
            args.conf_file = argv[1];
        }
    }

    return args;
}

int32_t main(int32_t argc, char *argv[])
{
    Arguments args = parseConfiguration(argc, argv);

    if (!args.deterministicOutput)
    {
        cout << "--=[ DEBORAH v" << DEBORAH_VERSION
             << " ]=--   (C) 2008-2011 Dip. Ing. dell'Informazione, University of Pisa (Italy)   ";
        cout << "(build: " << __DATE__ << ", " << __TIME__ << ")" << endl
             << endl;
    }

    Tandem t1;
    double lowerbound = 0.0, ludb = 0.0, ludb_pn = 0.0;

    if (args.ludb_evaluation)
    {
        LUDB_evaluation(
            t1,
            args.ludb_tree_order,
            args.ludb_tree_depth,
            args.ludb_max_good_combos
        );
        return 0;
    }
    if (args.lowerbound_evaluation)
    {
        LowerBound_evaluation(
            t1,
            args.ludb_tree_order,
            args.ludb_tree_depth,
            args.lowerbound_experimental,
            args.lowerbound_percentage
        );
        return 0;
    }
    if (args.ludb_nnested_evaluation)
    {
        LUDB_non_nested_evaluation(
            args.ludb_tree_order,
            args.ludb_tree_depth,
            args.ludb_max_good_combos,
            args.ludb_maxcuts,
            !args.compute_ludb_nnested_alt
        );
        return 0;
    }

    // check if a configuration file has been specified
    if (args.conf_file.empty())
    {
        std::cout << "Error: a scenario configuration file must be specified (or any -xxxx-evaluate switch)."
                  << std::endl;
        print_help();
        return -1;
    }

    // load the configuration file
    if (!args.deterministicOutput)
    {
        std::cout << "Parsing scenario configuration file \"" << args.conf_file << "\"... ";
    }
    else
    {
        std::cout << "Parsing scenario configuration file... ";
    }

    if (!t1.Load(args.conf_file))
    {
        std::cout << "error!" << std::endl;
        return DEBORAH_ERROR_CONFIG;
    }
    else
    {
        std::cout << "ok." << std::endl;
    }

    if ((args.flowrate_mult != 1.0) || (args.noderate_mult != 1.0))
    {
        double rmin, rmax;
        int32_t nmin, nmax;
        t1.analyzeProvisioning(rmin, nmin, rmax, nmax);
        std::cout << "Current over-provisioning status: min=" << rmin << " at node " << nmin + 1 << ",  max=" << rmax
                  << " at node " << nmax + 1 << "\n";
        std::cout << "Rescaling rates for over-provisioning: flows ratio=" << args.flowrate_mult << ",  nodes ratio="
                  << args.noderate_mult << "\n";

        bool res = t1.scaleProvisioning(args.flowrate_mult, args.noderate_mult);
        if (!res)
        {
            std::cout << "Error: Provisioning re-scaling parameters determine violation of constraints" << std::endl;
            return DEBORAH_ERROR_PROV;
        }
        t1.analyzeProvisioning(rmin, nmin, rmax, nmax);
        std::cout
            << "New over-provisioning status: min=" << rmin
            << " at node " << nmin + 1
            << "  max=" << rmax
            << " at node " << nmax + 1
            << std::endl;
    }

    // check if provisioning at the nodes is correct
    if (!t1.checkProvisioning())
    {
        auto details = t1.getUnderProvisionedNode();
        if(args.deterministicOutput)
        {
            // Use old error message to avoid breaking tests
            std::cout << "Error: Provisioning constraints not respected at node " << std::get<0>(details) << "\n";
        }
        else
        {
            std::cout
                << "Error: Provisioning constraints not respected at node " << std::get<0>(details)
                << ": flows expect " << std::get<1>(details)
                << " , node provides " << std::get<2>(details) << std::endl;
        }
        return DEBORAH_ERROR_PROV;
    }

    // check if user has specified a tagged flow explicitly
    if (args.TaggedFlow >= 0)
    {
        t1.SetTaggedFlow(args.TaggedFlow);
    }

    t1.Print();

    // check whether the tandem is nested
    bool nested = t1.isNested();

    // build the nesting tree
    t1.BuildNestingTree();

    if (nested)
    {
        std::cout << "\nAssociated nesting tree:" << endl;
        std::cout << "----------------------------------------------------" << endl;
        t1.PrintNestingTree();
        std::cout << "----------------------------------------------------" << endl;
    }

    // check if this tandem is equivalent to a sink-tree
    if (t1.isSinkTree())
    {
        std::cout << "This is a sink-tree." << endl;
    }

    std::cout << std::endl;

    // check if dimensional limits of the simplex library are respected
    if (args.compute_ludb && (t1.NumFlows() > NMAX))
    {
        std::cout << "LUDB computation disabled: too many t-nodes (" << t1.NumFlows() << "), max " << NMAX
                  << " allowed by simplex library." << std::endl;
        std::cout << "(see \"simplex.h\" in source code)" << std::endl;
        args.compute_ludb = false;
    }

    // before proceeding, make sure we actually have some job to do
    if (
        !args.compute_ludb && !args.compute_lowerbound &&
        !args.compute_ludb_pernode && !args.compute_output_curve &&
        !args.compute_eq_service_curve
    )
    {
        std::cout << "Nothing to do, exiting...\n" << std::endl;
        return 0;
    }

    // run the output curve computation, using LUDB
    if (args.compute_output_curve)
    {
        if (nested)
        {
            OutputArrivalCurveParameters param = t1.GetFoiOuputArrivalCurve(args.ludb_max_good_combos);
            std::cout
                << "Output arrival curve parameters: " << std::endl
                << "[" << std::endl
                << "\t{ \"burst\" : " << std::fixed << std::setprecision(6) << param.burst
                << " , \"rate\" : " << std::fixed << std::setprecision(6) << param.rate
                << " }" << std::endl
                << "]" << std::endl;
        }
        else
        {
            std::cout << "Running STA algorithm for output arrival curve computation" << std::endl;

            uint64_t nsimpl = 0;
            uint64_t ncuts = args.ludb_maxcuts;

            OutputCurve_non_nested_STA(t1, args.ludb_max_good_combos, nullptr, nsimpl, ncuts);
        }
    }

    // run the equivalent service curve computation, using LUDB
    if (args.compute_eq_service_curve)
    {
        if (nested)
        {
            PseudoAffine eq = t1.GetEquivalentServiceCurve(args.ludb_max_good_combos);
            std::cout
                << "Equivalent service curve: " << std::endl
                << "[" << std::endl
                << "\t" << eq.PrintJson() << std::endl
                << "]" << std::endl;
        }
        else
        {
            std::cout << "Running STA algorithm for output arrival curve computation" << std::endl;

            uint64_t nsimpl = 0;
            uint64_t ncuts = args.ludb_maxcuts;

            EquivalentServiceCurve_non_nested_STA(t1, args.ludb_max_good_combos, nullptr, nsimpl, ncuts);
        }
    }

    // run the LUDB algorithm
    if (args.compute_ludb)
    {
        if (nested)
        {
            ludb = t1.LUDB(args.ludb_max_good_combos, args.deterministicOutput);
        }
        else
        {
            std::cout
                << std::endl
                << "This tandem is not nested, using \""
                << (args.compute_ludb_nnested_alt ? "STA" : "MTA")
                << "\" algorithm."
                << std::endl;

            uint64_t nsimpl = 0;
            uint64_t ncuts = args.ludb_maxcuts;

            if (args.compute_ludb_nnested_alt)
            {
                ludb = LUDB_non_nested_STA(
                    t1,
                    args.ludb_max_good_combos,
                    nullptr,
                    nsimpl,
                    ncuts
                );
            }
            else
            {
                ludb = LUDB_non_nested_MSA(
                    t1,
                    args.ludb_max_good_combos,
                    nullptr,
                    nsimpl,
                    ncuts
                );
            }
        }
    }

    // run the LowerBound algorithm
    if (args.compute_lowerbound)
    {
        lowerbound = t1.LowerBound(
            false,
            args.combo_range_start,
            args.combo_range_len,
            args.lowerbound_percentage,
            args.deterministicOutput
        );
    }
    else if (args.lowerbound_experimental)
    {
        lowerbound = t1.LowerBound_experimental(args.lowerbound_percentage);
    }

    // run the per-node delay bound analysis
    if (args.compute_ludb_pernode)
    {
        ludb_pn = DelayBound_per_node(
            t1,
            args.ludb_max_good_combos
        );
    }

    // compare the delay values obtained from the two algorithms
    std::cout << std::endl;
    compare_results(
        ludb,
        lowerbound,
        ludb_pn,
        args.compute_ludb,
        args.compute_lowerbound || args.lowerbound_experimental
    );

    std::cout << std::endl;
    return 0;
}
