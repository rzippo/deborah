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
#include "tandem.h"
#include "simplex/simplex.h"
#include "timing.h"
#include <sstream>
#include <cmath>
#include <ctime>
#include <set>
#include <algorithm>
#include <memory>
#include <fstream>

// The following directive enables sanity checks in the LowerBound functions (recommended for testing)
#define TANDEM_DEBUG_CHECKS

// The following directive enables very verbose debug messages during the LowerBound computation
//#define TANDEM_DEBUG_LOWERBOUND

// The following directive enables debugging information in the LowerBound compute_CDF_I() method
//#define TANDEM_DEBUG_CDF_I

// The following directive enables debugging in LowerBound combinations cutting code
//#define TANDEM_DEBUG_LBOPT

// The directive below enables the printing of some info about each iteration during the LUDB computation
//#define TANDEM_VERBOSE_LUDB

// The directive below enables very verbose debug messages during the LUDB computation
//#define TANDEM_DEBUG_LUDB

// The directive below enables debug prints about the tandem cutting functions
//#define TANDEM_DEBUG_CUTS

// The directive below enables verbose debug prints about the output arrival curves computing functions
//#define TANDEM_DEBUG_OUTPUT_ARR_CURVES

//#define TANDEM_DEBUG_TIMINGS

#ifdef TANDEM_CATCH_EXCEPTIONS

#include <exception>
#include <iostream>
#include <iomanip>

#endif

namespace deborah
{

/* ************************************************************************************************************* */
/* ********************                     General purpose methods                         ******************** */
/* ************************************************************************************************************* */

Tandem::Tandem()
{
    num_nodes = 0;
    num_flows = 0;
    nesting_level = 0;
    nodes = {};
    flows = {};
    tagged_flow = 0;
    nesting_tree = nullptr;
    cached_ludb = -1.0;
}

Tandem::Tandem(uint64_t n_nodes, uint64_t n_flows)
{
    num_nodes = n_nodes;
    num_flows = n_flows;
    nesting_level = 0;
    nodes = std::vector<Node>(num_nodes);
    flows = std::vector<Flow>(num_flows);
    tagged_flow = 0;
    nesting_tree = nullptr;
    cached_ludb = -1.0;
}

Tandem::Tandem(const std::vector<Node>& v_nodes, const std::vector<Flow>& v_flows)
{
    num_nodes = v_nodes.size();
    num_flows = v_flows.size();
    nesting_level = 0;
    nodes = v_nodes;
    flows = v_flows;
    tagged_flow = 0;
    nesting_tree = nullptr;
    cached_ludb = -1.0;
}

Tandem::Tandem(const Tandem &t2)
{
    nesting_level = 0;
    nesting_tree = nullptr;
    num_nodes = t2.num_nodes;
    num_flows = t2.num_flows;
    tagged_flow = t2.tagged_flow;
    cached_ludb = t2.cached_ludb;

    if (num_nodes > 0)
    {
        nodes = std::vector<Node>(num_nodes);
        for (uint64_t i = 0; i < num_nodes; i++)
        {
            nodes[i] = t2.nodes[i];
        }
    }
    else
    {
        nodes = {};
    }

    if (num_flows > 0)
    {
        flows = std::vector<Flow>(num_flows);
        for (uint64_t i = 0; i < num_flows; i++)
        {
            flows[i] = t2.flows[i];
        }
    }
    else
    {
        flows = {};
    }
}

Tandem::~Tandem()
{
    Clear();
}

uint64_t Tandem::NumNodes() const
{
    return num_nodes;
}

uint64_t Tandem::NumFlows() const
{
    return num_flows;
}

void Tandem::Clear()
{
    nesting_tree = nullptr;
    nodes = {};
    flows = {};
    num_nodes = 0;
    num_flows = 0;
    nesting_level = 0;
    cached_ludb = -1.0;
}

Node *Tandem::getNode(uint64_t index)
{
    if (index < num_nodes)
    {
        return &nodes[index];
    }
    else
    {
        return nullptr;
    }
}

// Return flow at #index
Flow *Tandem::getFlow(uint64_t index)
{
    if (index < num_flows)
    {
        return &flows[index];
    }
    else
    {
        return nullptr;
    }
}

// Find a flow by UID
Flow *Tandem::getFlowId(uint64_t id)
{
    for (uint64_t i = 0; i < num_flows; i++)
    {
        if (flows[i].uid == id)
        {
            return &flows[i];
        }
    }
    // not found
    return nullptr;
}

// Find a flow by UID
uint64_t Tandem::getFlowIdx(uint64_t id)
{
    for (uint64_t i = 0; i < num_flows; i++)
    {
        if (flows[i].uid == id)
        {
            return i;
        }
    }
    // not found
    return FLOW_INVALID_ID;
}

/*	This method creates a balanced-tree topology.
	'bt.order' is the number of children for each t-node
	'bt.depth' is the number of nesting levels
	If bt.order is 1, then a sink-tree is generated.
*/
bool Tandem::CreateTree(const balanced_tree_descriptor &bt, RNG &rng)
{
    uint64_t curr_flows;
    uint64_t curr_node, nodes_span;
    bool make_sink_tree = false;

    if ((bt.depth <= 0) || (bt.order <= 0))
    {
        return false;
    }

    Clear();

    if (bt.order < 2)
    {
        make_sink_tree = true;
    }

    if (!make_sink_tree)
    {
        // balanced tree
        num_nodes = (uint64_t) pow(bt.order, bt.depth);
        num_flows = 1;
        for (int32_t i = 1; i <= bt.depth; i++)
        {
            num_flows += (uint64_t) pow(bt.order, i);
        }
    }
    else
    {
        // sink tree
        num_nodes = bt.depth;
        num_flows = num_nodes;
    }

    nodes = std::vector<Node>(num_nodes);
    if (nodes.empty())
    {
        return false;
    }
    flows = std::vector<Flow>(num_flows);
    if (flows.empty())
    {
        Clear();
        return false;
    }

    // create flows
    if (make_sink_tree)
    {
        // create sink tree
        for (uint32_t i = 0; i < num_flows; i++)
        {
            flows[i].src = i;
            flows[i].exit = num_flows - 1;
            flows[i].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
            flows[i].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
            flows[i].uid = i;
        }
    }
    else
    {
        // create balanced tree
        num_flows = 0;
        for (int32_t i = 0; i <= bt.depth; i++)
        {
            curr_flows = (uint64_t) pow(bt.order, i); // number of tnodes at this level
            nodes_span = num_nodes / curr_flows;
            curr_node = 0;
            for (uint64_t j = 0; j < curr_flows; j++)
            {
                flows[num_flows].src = curr_node;
                curr_node += nodes_span;
                flows[num_flows].exit = curr_node - 1;
                flows[num_flows].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
                flows[num_flows].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
                flows[num_flows].uid = num_flows;
                num_flows++;
            }
        }
    }

    // create nodes
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        double node_rate_min = 0.0;
        for (uint32_t j = 0; j < num_flows; j++)
        {
            if (isNodeCrossedByFlow(i, j))
            {
                node_rate_min += flows[j].rate;
            }
        }
        nodes[i].latency = rng.uniform(bt.node_latency_min, bt.node_latency_max);
        nodes[i].rate = node_rate_min * rng.uniform(bt.node_rate_overprov_min, bt.node_rate_overprov_max);
    }

    tagged_flow = 0;

    ComputeNestingLevels();

    return true;
}

/*	This method creates a tandem comprising 'n' nodes and 'k' flows
 *  for all possible combinations of (src, dst) nodes.
 *  Note: the resulting tandem may be non-nested.
 */
bool Tandem::CreateFullTandem(const balanced_tree_descriptor &bt, RNG &rng)
{
    bool random = false;
    Clear();
    if (bt.order < 3)
    {
        return false;
    }

    // allocate nodes
    num_nodes = bt.order;
    nodes = std::vector<Node>(num_nodes);

    // allocate flows
    num_flows = num_nodes * (num_nodes - 1) / 2;
    if (bt.depth > 0 && bt.depth < 100)
    {
        num_flows *= (float) bt.depth / 100.0;
        random = true;
    }
    else if (bt.depth < 0 && bt.depth > -num_flows)
    {
        num_flows = -bt.depth;
        random = false;
    }
    flows = std::vector<Flow>(num_flows);

    // configure flows
    uint64_t fid = 0;

    if (random)
    {
        while (fid < num_flows)
        { // only a subset of flows
            if (fid == 0)
            {
                // set tagged flow manually
                flows[fid].src = 0;
                flows[fid].exit = num_nodes - 1;
            }
            else
            {
                flows[fid].src = rng.uniform(0, num_nodes - 1);
                flows[fid].exit = rng.uniform(flows[fid].src + 1, num_nodes);
                if (flows[fid].src == flows[fid].exit)
                {
                    continue;
                }
                bool found = false;
                for (uint64_t i = 0; i < fid; i++)
                {
                    if (flows[i].src == flows[fid].src && flows[i].exit == flows[fid].exit)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                {
                    continue;
                }
            }
            flows[fid].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
            flows[fid].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
            flows[fid].uid = fid;
            fid++;
        }
    }
    else
    {
        for (int32_t flen = num_nodes; flen >= 2; flen--)
        { // all possible flows
            int32_t nf = num_nodes - flen + 1;
            for (int32_t i = 0; i < nf; i++)
            {
                if (fid < num_flows)
                {
                    flows[fid].src = i;
                    flows[fid].exit = i + flen - 1;
                    flows[fid].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
                    flows[fid].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
                    flows[fid].uid = fid;
                }
                fid++;
            }
        }
    }

    tagged_flow = 0;
    // sanity check
    if (fid != num_flows)
    {
        std::cout
            << "Tandem::CreateFullTandem(): BUG: " << fid
            << " flows created, expected count = " << num_flows
            << std::endl;
    }

    // configure nodes
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        double node_rate_min = 0.0;
        for (uint32_t j = 0; j < num_flows; j++)
        {
            if (isNodeCrossedByFlow(i, j))
            {
                node_rate_min += flows[j].rate;
            }
        }
        nodes[i].latency = rng.uniform(bt.node_latency_min, bt.node_latency_max);
        nodes[i].rate = node_rate_min * rng.uniform(bt.node_rate_overprov_min, bt.node_rate_overprov_max);
    }

    return true;
}

/*	This method sets the current object according to the configuration file
 *  In case of error, it returns an empy object
 */
bool Tandem::Load(const std::string &config_filename)
{
    int32_t tagged = -1;

    int32_t n_nodes = 0, n_flows = 0;

    std::ifstream f;
    f.open(config_filename, std::ifstream::in);

    if (f.fail())
    {
        return false;
    }

    //Delete previous data, if present
    Clear();

    while (!f.eof())
    {
        uint64_t src_node, exit_node;
        double theta, R, sigma, rho;

        std::string line;
        std::getline(f, line);

        if (line.empty())
        {
            continue;
        }

        if (line[0] == '#')
        {
            continue;
        }

        std::stringstream tokenizer(line);
        std::string directive;
        tokenizer >> directive;

        if (directive == "TANDEM")
        {
            tokenizer >> num_nodes >> num_flows;
            if ((num_nodes <= 0) || (num_flows <= 0))
            {
                break;
            }
            //Allocate structures
            nodes = std::vector<Node>(num_nodes);
            flows = std::vector<Flow>(num_flows);
            continue;
        }

        if (directive == "NODE")
        {
            uint64_t node_id;
            if (n_nodes >= num_nodes)
            {
                continue;
            }
            theta = 0.0;
            R = 0.0;

            tokenizer >> std::fixed >> node_id >> theta >> R;

            if (node_id - 1 > num_nodes)
            {
                std::cout << "Warning: invalid node ID " << node_id << " in configuration file, skipped";
                continue;
            }
            //Add node to tandem
            nodes[node_id - 1].latency = theta;
            nodes[node_id - 1].rate = R;
            n_nodes++;
            continue;
        }

        if (directive == "FLOW" || directive == "TFLOW")
        {
            if (n_flows >= num_flows)
            {
                continue;
            }
            // check for tagged-flow directive
            if ((directive[0] == 'T') || (directive[0] == 't'))
            {
                tagged = n_flows;
            }
            src_node = -1;
            exit_node = -1;
            sigma = 0.0;
            rho = 0.0;

            tokenizer >> std::fixed >> src_node >> exit_node >> sigma >> rho;

            if ((src_node == 1) && (exit_node == num_flows))
            {
                tagged_flow = n_flows;
            }
            //Add flow to tandem
            flows[n_flows].src = src_node - 1;
            flows[n_flows].exit = exit_node - 1;
            flows[n_flows].burst = sigma;
            flows[n_flows].rate = rho;
            flows[n_flows].uid = n_flows;
            n_flows++;
            continue;
        }
    }

    //Close configuration file
    f.close();

    //Check results
    bool result = true;
    if ((n_nodes <= 0) || (n_flows <= 0))
    {
        result = false;
    }
    if ((n_nodes < num_nodes) || (n_flows < num_flows))
    {
        std::cout << "Error: config file is not consistent: missing declarations for nodes and/or flows"
                  << std::endl;
        result = false;
    }

    if (tagged >= 0)
    {
        tagged_flow = tagged;
    }

    //Now update the auxiliary data structures (nesting tree, etc)
    ComputeNestingLevels();

    //In case of error, empty the current Tandem
    if (!result)
    {
        Clear();
    }

    return result;
}

/*	Writes the current tandem topology to a configuration file
*/
bool Tandem::Save(const std::string &config_filename)
{
    std::ofstream config;

    // some preliminary checks
    if (config_filename.empty())
    {
        return false;
    }
    if (!num_nodes || !num_flows)
    {
        return false;
    }

    config.open(config_filename, std::ofstream::out);
    if (config.fail())
    {
        return false;
    }

    // write header
    config << "# DEBORAH CONFIGURATION FILE v1.2\n\n";
    config << "# TANDEM <nodes #> <flows #>\n";
    config << "TANDEM " << num_nodes << " " << num_flows << "\n\n";

    // write nodes
    config << "# NODE <node #> <latency> <rate>\n";
    for (int32_t i = 0; i < num_nodes; i++)
    {
        config << "NODE	" << i + 1 << "\t" << nodes[i].latency << "\t" << nodes[i].rate << "\n";
    }
    config << "\n";

    // write flows
    config << "# FLOW <source node> <sink node> <burst> <rate>\n";
    for (int32_t i = 0; i < num_flows; i++)
    {
        config << "FLOW	" << flows[i].src + 1 << "\t" << flows[i].exit + 1 << "\t" << flows[i].burst << "\t"
               << flows[i].rate << "\n";
    }
    config << "\n";

    config.close();
    return true;
}

void Tandem::Print()
{
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "Tandem: " << num_nodes << " nodes, " << num_flows << " flows, nesting_level=" << nesting_level
              << "\n";
    if (num_nodes > 0 && num_flows > 0)
    {
        std::cout << "        Tagged flow = (" << flows[tagged_flow].src + 1 << "," << flows[tagged_flow].exit + 1
                  << ")  flow_id=#" << tagged_flow << "\n";
        std::cout << "        Tandem is" << (isNested() ? " " : " not ") << "nested." << std::endl;
        for (int32_t i = 0; i < num_nodes; i++)
        {
            std::cout
                << "NODE " << std::setw(2) << std::setfill('0') << i + 1
                << ": latency=" << std::fixed << std::setprecision(2) << nodes[i].latency
                << ", rate=" << std::fixed << std::setprecision(2) << nodes[i].rate
                << std::endl;
        }
        for (int32_t i = 0; i < num_flows; i++)
        {
            std::cout
                << "FLOW " << std::setw(2) << std::setfill('0') << i
                << " (" << flows[i].src + 1
                << "," << flows[i].exit + 1
                << "): burst=" << flows[i].burst
                << ", rate=" << flows[i].rate
                << "; nesting_level=" << flows[i].nesting_level
                << "; uid=" << std::hex << flows[i].uid << std::dec
                << "\n";
        }
    }
    std::cout << "------------------------------------------------------------" << std::endl;
}

bool Tandem::isSinkTree()
{
    if (num_flows != num_nodes)
    {
        return false;
    }
    for (int32_t i = 0; i < num_flows; i++)
    {
        if ((flows[i].src != i) || (flows[i].exit != num_nodes - 1))
        {
            return false;
        }
    }
    return true;
}

bool Tandem::isNested()
{
    // check that no two flows intersect each other
    for (uint32_t i = 0; i < num_flows; i++)
    {
        for (uint32_t j = 0; j < num_flows; j++)
        {
            if (i == j)
            {
                continue;
            }
            if (flows[i].intersects(flows[j]))
            {
                return false;
            }
        }
    }
    return true;
}

bool Tandem::checkProvisioning()
{
    double rate;
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        rate = 0.0;
        for (uint32_t j = 0; j < num_flows; j++)
        {
            if (isNodeCrossedByFlow(i, j))
            {
                rate += flows[j].rate;
            }
        }
        if (rate > nodes[i].rate)
        {
            return false;
        }
    }
    return true;
}

/// Dual of checkProvisioning(), used to detect the point of error.
/// \return a tuple containing 1) position of the under-provisioned node 2) the rate expected by flows 3) the rate provided by the node
std::tuple<int32_t, double, double> Tandem::getUnderProvisionedNode()
{
    double rate;
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        rate = 0.0;
        for (uint32_t j = 0; j < num_flows; j++)
        {
            if (isNodeCrossedByFlow(i, j))
            {
                rate += flows[j].rate;
            }
        }
        if (rate > nodes[i].rate)
        {
            auto ret = std::tuple(i + 1, rate, nodes[i].rate);
            return ret;
        }
    }

    // This is reached if there is no under-provisioned node
    return std::tuple(0, 0.0, 0.0);
}

bool Tandem::scaleProvisioning(double flow_mult, double node_mult)
{
    if ((flow_mult <= 0.0) || (node_mult <= 0.0))
    {
        return false;
    }
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        nodes[i].rate *= node_mult;
        for (uint32_t j = 0; j < num_flows; j++)
        {
            if (isNodeCrossedByFlow(i, j))
            {
                flows[j].rate *= flow_mult;
            }
        }
    }
    return checkProvisioning();
}

bool Tandem::analyzeProvisioning(double &min_prov, int32_t &nodemin, double &max_prov, int32_t &nodemax)
{
    if (!num_nodes || !num_flows)
    {
        return false;
    }
    double rate, overprov;
    min_prov = 10000000.0;
    max_prov = 0.0;
    nodemin = -1;
    nodemax = -1;
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        rate = 0.0;
        for (uint32_t j = 0; j < num_flows; j++)
        {
            if (isNodeCrossedByFlow(i, j))
            {
                rate += flows[j].rate;
            }
        }
        overprov = nodes[i].rate / rate;
        if (overprov < min_prov)
        {
            min_prov = overprov;
            nodemin = i;
        }
        if (overprov > max_prov)
        {
            max_prov = overprov;
            nodemax = i;
        }
    }
    return true;
}

/* ************************************************************************************************************* */
/* ********************                       LUDB related methods                          ******************** */
/* ************************************************************************************************************* */

/* This method computes the nesting level of each flow in the tandem
 * For simplicity, values are stored in each flow object
 */
void Tandem::ComputeNestingLevels()
{
    if ((num_flows <= 0) || (flows.empty()))
    {
        return;
    }

    nesting_level = 0;
    for (int32_t i = 0; i < num_flows; i++)
    {
        flows[i].nesting_level = 0;
        for (int32_t j = 0; j < num_flows; j++)
        {
            if (flows[i].is_nested(flows[j]))
            {
                flows[i].nesting_level++;
                if (flows[i].nesting_level > nesting_level)
                {
                    nesting_level = flows[i].nesting_level;
                }
            }
        }
    }
}

/*	This method builds the complete nesting tree for this tandem
 */
bool Tandem::BuildNestingTree()
{
    if (flows.empty())
    {
        return false;
    }
    if (nodes.empty())
    {
        return false;
    }

    nesting_tree = build_nesting_subtree(tagged_flow);
    if (nesting_tree == nullptr)
    {
        return false;
    }

    return true;
}

/*	Auxiliary method which builds a sub-tree
*/
std::shared_ptr<TNode> Tandem::build_nesting_subtree(uint64_t flow_id)
{
    if (flow_id >= num_flows)
    {
        return nullptr;
    }

    //pointer to considered flow (h,k)
    Flow *f = &flows[flow_id];
    int32_t nest_level = f->nesting_level;

    std::shared_ptr<TNode> tnode = std::make_shared<TNode>();
    tnode->flow = flow_id;
    tnode->pflow = f;

    uint64_t n_leaves;
    uint64_t n_tnodes;

    bool contained;

    //pointer to flow(s) directly nested into (h,k)
    Flow *nf;
    //temp pointer to a node
    Node *pn;

    //Find number of directly nested flows
    n_tnodes = 0;
    for (int32_t i = 0; i < num_flows; i++)
    {
        if (flows[i].nesting_level == (nest_level + 1))
        {
            if (flows[i].is_nested(*f))
            {
                n_tnodes++;
            }
        }
    }
    tnode->num_tnodes = n_tnodes;

    //Create child nodes for directly nested flows
    if (tnode->num_tnodes > 0)
    {
        tnode->children = std::vector<std::shared_ptr<TNode>>(tnode->num_tnodes);
        n_tnodes = 0;
        for (int32_t i = 0; i < num_flows; i++)
        {
            if (flows[i].nesting_level == (nest_level + 1))
            {
                if (flows[i].is_nested(*f))
                {
                    tnode->children[n_tnodes] = build_nesting_subtree(i);
                    //check errors in the costruction of the subtree
                    if (tnode->children[n_tnodes] == nullptr)
                    {
                        return nullptr;
                    }
                    n_tnodes++;
                }
            }
        }
    }
    else
    {
        tnode->children = {};
    }

    //Find number of leaves
    n_leaves = 0;
    for (int32_t i = f->src; i <= f->exit; i++)
    {
        contained = false;
        for (int32_t j = 0; j < tnode->num_tnodes; j++)
        {
            nf = &flows[tnode->children[j]->flow];
            if (nf->contains(i))
            {
                contained = true;
            }
        }
        if (!contained)
        {
            n_leaves++;
        }
    }
    tnode->num_leaves = n_leaves;

    //Create the actual leaves
    if (tnode->num_leaves > 0)
    {
        tnode->leaves = new uint64_t[tnode->num_leaves];
        n_leaves = 0;
        for (int32_t i = f->src; i <= f->exit; i++)
        {
            contained = false;
            for (int32_t j = 0; j < tnode->num_tnodes; j++)
            {
                nf = &flows[tnode->children[j]->flow];
                if (nf->contains(i))
                {
                    contained = true;
                }
            }
            if (!contained)
            {
                tnode->leaves[n_leaves] = i;
                n_leaves++;
            }
        }
    }

    // Compute service curve PI_C(h,k) for leaves
    if (tnode->has_leaf())
    {
        pn = getNode(tnode->leaves[0]);
        tnode->pi_c.setLatencyRate(pn->latency, pn->rate);
        for (int32_t i = 1; i < tnode->num_leaves; i++)
        {
            pn = getNode(tnode->leaves[i]);
            tnode->pi_c.convolve(pn->latency, pn->rate);
        }
    }

    // Initialize service curves PI(h,k) and PI_EQ(h,k) for the current t-node
    if (tnode->has_leaf())
    {
        // Initialize the service curve at this "child flow" with the PI_C
        tnode->pi = new ParametrizedPseudoAffine(num_flows);
        tnode->pi_eq = new ParametrizedPseudoAffine(num_flows);
        tnode->pi->MakeLatencyRate(tnode->pi_c);
    }
    else
    {
        // Initialize an empty PI and PI_EQ
        tnode->pi = new ParametrizedPseudoAffine(num_flows);
        tnode->pi_eq = new ParametrizedPseudoAffine(num_flows);
    }
    // When you compute PI_EQ you introduce new indexes for the current node!

    // Build the set of sets of indexes of the current t-node
    tnode->index_set.clear();
    tnode->index_schema.clear();
    tnode->index_nodes.clear();
    tnode->index_nodes.push_back(tnode->flow);
    // first add the (single) index for our own t-node
    tnode->index_schema.push_back(1);
    // initialize the index limit as the number of stages of the equivalent service curve of this t-node
    uint64_t index_limit = tnode->has_leaf() ? 1 : 0;
    // then add the sets of indexes from every child t-node
    for (int32_t i = 0; i < tnode->num_tnodes; i++)
    {
        std::vector<uint64_t> child_schema = tnode->children[i]->index_schema;
        IndexSet child_set = tnode->children[i]->index_set;
        for (int32_t j = 0; j < child_set.size(); j++)
        {
            tnode->index_set.push_back(child_set[j]);
            tnode->index_nodes.push_back(tnode->children[i]->index_nodes[j]);
        }
        tnode->index_schema.push_back(child_set.size());
        // add the number of stages of the child's equivalent service curve
        // note: if the child has no leaf, then the actual number of stages of his PI is
        // child_set[0]-1 because its index limit takes into account the extra case "max()=0"
        if (tnode->children[i]->has_leaf())
        {
            index_limit += child_set[0];
        }
        else
        {
            index_limit += child_set[0] - 1;
        }
    }
    // if we have no leaf, count ONE MORE choice representing "case max()=0"
    if (!tnode->has_leaf())
    {
        index_limit++;
    }
    // head-insert the determined index limit for our t-node
    tnode->index_set.insert(tnode->index_set.begin(), index_limit);

    return tnode;
}

void Tandem::PrintNestingTree()
{
    if (nesting_tree != nullptr)
    {
        nesting_tree->Print(true);
    }
    else
    {
        std::cout << "Nesting Tree not computed for this tandem!" << std::endl;
    }
}

uint64_t Tandem::LUDB_TotalCombos()
{
    if (nesting_tree == nullptr)
    {
        return 0;
    }
    return nesting_tree->getNumCombos();
}

void Tandem::PrintAllCombinations()
{
    if (nesting_tree == nullptr)
    {
        return;
    }

    IndexSet set;
    uint64_t num = nesting_tree->getNumCombos();
    for (uint64_t i = 0; i < num; i++)
    {
        nesting_tree->getCombo(set, i);
        std::cout << std::setw(3) << std::setfill('0') << i << ": [ ";
        for (auto j : set)
        {
            std::cout << j;
        }
        std::cout << "]\n";
    }

    // alternate enumeration method
    std::cout << "--- alternate method\n";
    nesting_tree->getCombo(set, 0);
    for (uint64_t i = 0; i < num; i++)
    {
        int32_t h;
        std::cout << std::setfill('0') << std::setw(3) << i << ": [ ";
        for (auto j : set)
        {
            std::cout << j << " ";
        }
        std::cout << "]" << std::endl;
        h = nesting_tree->incCombo(set, -1);
    }
}

/// Performs the LUDB computation.
/// Optimized to reduce the number of simplexes solved.
/// \param max_good_combos
/// \param deterministicOutput
/// \return The computed LUDB
double Tandem::LUDB(int32_t max_good_combos, bool deterministicOutput)
{
    std::vector<simplex::Constraint> constraints;
    uint64_t n_simplexes;
    simplex::Solution solution(num_flows);
    simplex::Objective delay(num_flows);
    IndexSetVector good_combos;
    clock_t t_start, t_end;

    if (nesting_tree == nullptr)
    {
        return -1.0;
    }
    if (num_flows > NMAX)
    {
        return -2.0;
    }

    std::cout
        << ((max_good_combos != LUDB_MODE_EXACT) ? "Heuristic" : "Exact")
        << " LUDB algorithm started:"
        << std::endl;

    if (max_good_combos != LUDB_MODE_EXACT)
    {
        std::cout << "Maximum number of combinations returned by each t-node: " << max_good_combos << std::endl;
    }
    std::cout << std::endl << "----------------------------------------------------" << std::endl;

    double t_elapsed;
    t_start = clock();
    n_simplexes = nesting_tree->computeLUDB(
        delay, solution, num_flows, good_combos, max_good_combos,
        true
    );
    t_end = clock();
    t_elapsed = ((double) t_end - (double) t_start) / (double) CLOCKS_PER_SEC;

    std::cout << "LUDB = " << std::setprecision(6) << solution.getConstant() << std::endl;
    if (!deterministicOutput)
    {
        std::cout << "Computation time: " << std::setprecision(3) << t_elapsed << " seconds" << std::endl;
    }
    std::cout << "Simplexes computed: " << n_simplexes << std::endl;
    std::cout << "Delay bound expression: ";
    delay.Print();
    std::cout << std::endl;
    std::cout << "Solution:" << std::endl;
    PrintSolution(solution, tagged_flow);
    std::cout << std::endl;
    cached_ludb = solution.getConstant();
    return solution.getConstant();
}

/// Non-printing version of the LUDB computation, to be used as part of other algorithms.
/// \param max_good_combos
/// \param n_simplexes
/// \return The computed LUDB
double Tandem::LUDB_quiet(int32_t max_good_combos, uint64_t &n_simplexes)
{
    return LUDB_quiet(nesting_tree, max_good_combos, n_simplexes);
}

/// Non-printing version of the LUDB computation, to be used as part of other algorithms.
/// \param max_good_combos
/// \param n_simplexes
/// \return The computed LUDB
double Tandem::LUDB_quiet(const std::shared_ptr<TNode> &root, int32_t max_good_combos, uint64_t &n_simplexes)
{
    std::vector<simplex::Constraint> constraints;
    simplex::Solution solution(num_flows);
    simplex::Objective delay(num_flows);
    IndexSetVector good_combos;

    if (root == nullptr)
    {
        return TANDEM_LUDB_ERROR_NESTINGTREE;
    }

    n_simplexes = root->computeLUDB(delay, solution, num_flows, good_combos, max_good_combos, false);

    cached_ludb = solution.getConstant();
    return solution.getConstant();
}

OutputArrivalCurveParameters Tandem::GetFoiOuputArrivalCurve(int32_t mode)
{
    Tandem tt = *this; // a copy of the original sub-tandem
    Flow *foi = tt.getFlow(tt.GetTaggedFlow());

    Flow outFlow = *foi;
    tt.aggregate_flows(false);
    tt.BuildNestingTree();

    if (!tt.compute_output_arrival_curve(foi->uid, &outFlow, mode))
    {
        exit(DEBORAH_ERROR_SAFECHECK);
    }

    OutputArrivalCurveParameters parameters{
        .burst = outFlow.burst,
        .rate = outFlow.rate
    };

    return parameters;
}

PseudoAffine Tandem::GetEquivalentServiceCurve(int32_t max_good_combos)
{
    std::vector<simplex::Constraint> constraints;
    simplex::Solution solution(num_flows);
    simplex::Objective delay(num_flows);
    IndexSetVector good_combos;

    if (nesting_tree == nullptr)
    {
        throw;
    }

    nesting_tree->computeLUDB(delay, solution, num_flows, good_combos, max_good_combos, false);

    auto indexSet = good_combos[0]; //pick a random one
    auto pa = nesting_tree->computeEquivalentServiceCurve(indexSet);

    return PseudoAffine(pa, solution);
}

/*	This method prints a solution vector for the LUDB algorithm
*/
void Tandem::PrintSolution(simplex::Solution &solution, uint64_t dead_var)
{
    int32_t nvars = solution.getNumVariables();

    for (int32_t i = 0; i < nvars; i++)
    {
        if (i == dead_var)
        {
            continue;
        }

        auto flowString = flows[i].toString();
        std::cout
            << "S" << flowString << " = "
            << std::fixed << std::setprecision(6) << solution.getCoefficient(i)
            << std::endl;
    }
}

bool Tandem::SetTaggedFlow(uint64_t tf)
{
    if (tf >= num_flows)
    {
        return false;
    }
    tagged_flow = tf;
    return true;
}

uint64_t Tandem::GetTaggedFlow() const
{
    return tagged_flow;
}

/*	Joins a vector of sub-tandems into the original tandem
 * 	nearly reverting the operation of cut()
 * 	Return: one if the resulting tandem is nested, zero if not.
 * 	Negative values in case of error.
 */
int32_t Tandem::join(std::vector<Tandem> &cuts)
{
    Clear();

    uint64_t size = cuts.size();
    if (!size)
    {
        return DEBORAH_ERROR_SIZE;
    }
    for (uint64_t i = 0; i < size; i++)
    {
        num_nodes += cuts[i].num_nodes;
        num_flows += cuts[i].num_flows;
    }
    if (!num_nodes || (num_flows <= size - 1))
    {
        return DEBORAH_ERROR_SAFECHECK;
    }
    num_flows -= size - 1; // count tagged flow only once
    flows = std::vector<Flow>(num_flows);
    nodes = std::vector<Node>(num_nodes);

    //	if (!flows || !nodes)
    //		return DEBORAH_ERROR_ALLOC;

    // reserve flow #0 for tagged flow
    uint64_t curr_flow = 1;
    uint64_t curr_node = 0;
    for (uint64_t i = 0; i < size; i++)
    {
        Tandem &ti = cuts[i];
        uint64_t tf = ti.tagged_flow;
        // compile flows
        for (uint64_t j = 0; j < ti.num_flows; j++)
        {
            if (j == tf)
            {
                continue;
            } // skip tagged flow
            flows[curr_flow].burst = ti.flows[j].burst;
            flows[curr_flow].rate = ti.flows[j].rate;
            //flows_arr[curr_flow].uid = ti.flows_arr[j].uid;
            flows[curr_flow].uid = curr_flow;
            flows[curr_flow].src = curr_node + ti.flows[j].src;
            flows[curr_flow].exit = curr_node + ti.flows[j].exit;
            curr_flow++;
        }
        // compile nodes
        for (uint64_t j = 0; j < cuts[i].num_nodes; j++)
        {
            nodes[curr_node].latency = cuts[i].nodes[j].latency;
            nodes[curr_node].rate = cuts[i].nodes[j].rate;
            curr_node++;
        }
    }
    // finally set the tagged flow
    tagged_flow = 0;
    Flow *tf = cuts[0].getFlow(cuts[0].tagged_flow);
    flows[tagged_flow].burst = tf->burst;
    flows[tagged_flow].rate = tf->rate;
    flows[tagged_flow].src = tf->src;
    flows[tagged_flow].exit = num_nodes - 1;
    //flows_arr[tagged_flow].uid = tf->uid;
    flows[tagged_flow].uid = 0;

    if (isNested())
    {
        ComputeNestingLevels();
        return 1;
    }

    return 0;
}

/* Compute the exhaustive set of cuts for this (non-nested) tandem
 * Cuts longer than 'max_length' nodes will be pruned
 */
std::vector<CutNodeList> Tandem::compute_cutsets()
{
    std::vector<CutNodeList> result;
    // for non-nested tandems only
    if (isNested())
    {
        return result;
    }
    std::set<std::vector<uint64_t>> cutpos;

    BasicTimer tcut, tcomp;

    tcut.mark(true);
    // fill-in the conflicts matrix
    for (uint64_t i = 0; i < num_flows; i++)
    {
        for (uint64_t j = 0; j < num_flows; j++)
        {
            uint64_t s1 = flows[i].src;
            uint64_t e1 = flows[i].exit;
            uint64_t s2 = flows[j].src;
            uint64_t e2 = flows[j].exit;
            if (s1 < s2 && s2 <= e1 && e1 < e2)
            {
                //std::cout << "intersection %02i with " << i,j << ": ";
                // flows intersect between nodes [s2,e1] included
                std::vector<uint64_t> alt;
                for (uint64_t k = s2; k <= e1 + 1; k++)
                {
                    alt.push_back(k);
                    //std::cout << "" << k << " ";
                    //node_tag[k]++;
                }
                //std::cout << std::endl;
                cutpos.insert(alt);
            }
        }
    }

    /*
          for(int32_t i=0; i<num_deps; i++)
          {
              std::cout << "d" << i << ":  ";
              for(int32_t j=0; j<cutpos[i].size(); j++) std::cout << "" << cutpos[i][j] << " ";
              std::cout << std::endl;
          }
          */

    // setup the node vs dependency matrix
    uint64_t num_deps = cutpos.size(); // number of dependencies
    std::vector<std::set<uint64_t> > nodedeps(num_nodes);

    // fill-in
    for (uint64_t n = 0; n < num_nodes; n++)
    {
        //std::cout << "node " << n << ":  ";
        uint64_t d = 0;
        for (const auto &cutpo : cutpos)
        {
            if (std::find(cutpo.begin(), cutpo.end(), n) != cutpo.end())
            {
                nodedeps[n].insert(d);
                //std::cout << "d" << d << " ";
            }
            d++;
        }
        //std::cout << std::endl;
    }

#ifdef TANDEM_DEBUG_TIMINGS
    tcut.mark(false);
std::cout
<< "compute_cutsets(): t_dep=" << std::setprecision(3) << tcut.time()
<< "  n_deps=" << num_deps
<< std::endl;
#endif

    // compute the sets of cuts
#ifdef TANDEM_CUTCOMP_OLDVERSION
    bool res[num_nodes];
for (uint64_t i = 0; i < num_nodes; i++)
res[i] = false;
std::set<uint64_t> cdeps;
#endif

    std::set<std::set<uint64_t>> allres;

    compstatus stat{
        .cres = 0x0000,
        .cnode = 0
    };

    if (num_deps > MAXDEPS * 32)
    {
        std::cout << "Tandem::compute_cutsets(): ERROR, max " << MAXDEPS * 32 << " dependencies supported ("
                  << num_deps << " found)." << std::endl;
        return result;
    }
    std::vector<depsse> alldeps(num_nodes);

    int32_t a = num_deps / MAXDEPS;
    int32_t b = num_deps % MAXDEPS;
    //std::cout << "num_deps=%i   a=%i b=" << num_deps,a,b << "\n";
    for (uint64_t i = 0; i < num_nodes; i++)
    {
        for (int32_t j = 0; j <= a; j++)
        {
            alldeps[i].dep[j] = 0;
        }
        for (int32_t k = 0; k < 32 - b; k++)
        {
            alldeps[i].dep[a] >>= 1;
            alldeps[i].dep[a] |= 0x80000000;
        }
        for (int32_t j = a + 1; j < MAXDEPS; j++)
        {
            alldeps[i].dep[j] = 0xFFFFFFFF;
        }
        for (auto jt = nodedeps[i].begin(); jt != nodedeps[i].end(); jt++)
        {
            uint64_t nd = *jt;
            int32_t a = nd >> 5;
            int32_t b = nd & 0x1F;
            //std::cout << "node %02i: dep %03i --> a=%i b=" << i,nd,a,b << "\n";
            alldeps[i].dep[a] |= (1L << b);
        }
        //std::cout << "alldeps[%02i] = ",i); for(int32_t j=0;j<MAXDEPS;j++) printf("" << alldeps[i].dep[j] << " "; std::cout << std::endl;
    }
    for (auto &i : stat.cdep.dep)
    {
        i = 0;
    }
    tcomp.mark(true);
    cutcomp(&stat, alldeps, allres);

#ifdef TANDEM_DEBUG_TIMINGS
    tcomp.mark(false);
std::cout
<< "compute_cutsets(): t_comp_new=" << std::setprecision(3) << tcomp.time()
<< "  n_sets=" << allres.size()
<< std::endl;
#endif

#ifdef TANDEM_CUTCOMP_OLDVERSION
    tcomp.mark(true);
cutcomp_slow(res, num_deps, nodedeps, allres);
#ifdef TANDEM_DEBUG_TIMINGS
tcomp.mark(false);
std::cout << "compute_cutsets(): t_comp_old=%.3lf  n_sets=" << tcomp.time(), allres.size() << "\n";
#endif
#endif

    // eliminate supersets in 'allres'
    for (auto it = allres.begin(); it != allres.end(); it++)
    {
        CutNodeList cn;
        bool found = false;
        // find any subsets of this
        for (auto kt = allres.begin(); kt != allres.end(); kt++)
        {
            if (it == kt)
            {
                continue;
            }
            if (includes(it->begin(), it->end(), kt->begin(), kt->end()))
            {
                found = true;
                break;
            }
        }
        if (found)
        {
            continue;
        }

        for (auto jt : *it)
        {
            cn.push_back(jt);
        }
        cn.push_back(num_nodes);
        result.push_back(cn);
    }

#ifdef TANDEM_DEBUG_CUTS
    std::cout << "Tandem::compute_cutsets():  " << result.size() << " sets found:" << std::endl;
for (uint64_t i = 0; i < result.size(); i++)
{
std::cout << "    set " << std::setw(2) << std::setfill('0') << i << " = { ";
for (uint64_t j = 0; j < result[i].size(); j++)
{
std::cout << result[i][j] << " ";
}
std::cout << "}" << std::endl;
}
#endif

    return result;
}

/// todo: what does this do? CHECK: filter to cuts not longer than (min_len + num_cuts)
/// \param cutsets
/// \param deltalen
void Tandem::filter_cutsets(std::vector<CutNodeList> &cutsets, int32_t deltalen)
{
    std::vector<CutNodeList> result;
    uint64_t len_min = 0xFFFFFFFF;
    uint64_t n = cutsets.size();
    for (uint64_t i = 0; i < n; i++)
    {
        uint64_t len = cutsets[i].size();
        if (len < len_min)
        {
            len_min = len;
        }
    }

    uint64_t len_max = len_min + deltalen;

    for (uint64_t i = 0; i < n; i++)
    {
        if (cutsets[i].size() <= len_max)
        {
            result.push_back(cutsets[i]);
        }
    }
    cutsets = result;
}

bool Tandem::cutcomp_slow(
    bool *res2, uint64_t ndeps, std::set<uint64_t> nodedeps[],
    std::set<std::set<uint64_t>> &allres
)
{
    if (res2 == nullptr)
    {
        return false;
    }

    // current set of satisfied dependencies
    std::vector<std::byte> cdeps(ndeps);
    std::fill(cdeps.begin(), cdeps.end(), std::byte{0x00});

    // highest order node which has been checked so far
    uint64_t currnode = 0;

    // compute cdeps
    for (uint64_t i = 1; i < num_nodes; i++)
    {
        if (!res2[i])
        {
            continue;
        }
        currnode = i;
        std::set<uint64_t> &d = nodedeps[i]; // dependencies satisfied by node *it
        for (auto jt = d.begin(); jt != d.end(); jt++)
        {
            cdeps[*jt] = std::byte{0x01};
        }
    }

    //std::cout << "res={ ");for(std::set<uint64_t>::iterator it=res.begin(); it!=res.end(); it++)std::cout << "" << *it << " ";printf("}" << std::endl;
    //std::cout << "cdeps={ ");for(int32_t i=0; i<ndeps; i++)if(cdeps[i]) std::cout << "" << i << " ";printf("}" << std::endl;
    for (uint64_t n = currnode + 1; n < num_nodes; n++)
    {
        //std::cout << "check node " << n << "\n";
        // don't consider nodes already into result
        if (res2[n])
        {
            continue;
        }
        // avoid the case of three cuts at consecutive nodes
        if ((n >= 2) && res2[n - 1] && res2[n - 2])
        {
            continue;
        }
        // check if this node satisfies new dependencies
        std::set<uint64_t> &dep = nodedeps[n];
        int32_t found = -1;
        for (auto it = dep.begin(); it != dep.end(); it++)
        {
            if (!((bool) cdeps[*it]))
            {
                found = *it;
                break;
            }
        }
        if (found < 0)
        {
            continue;
        } // no new dependencies, prune
        // add node n
        //std::cout << "add node %i (new dep=" << n,found << ") --> ";
        res2[n] = true;

        std::vector<std::byte> ddeps(ndeps);
        std::copy(cdeps.begin(), cdeps.end(), ddeps.begin());

        // update list of satisfied dependencies
        for (auto jt = dep.begin(); jt != dep.end(); jt++)
        {
            ddeps[*jt] = std::byte{0x01};
        }
        // check if all dependencies are now satisfied (termination condition)
        bool finished = true;
        for (uint64_t i = 0; i < ndeps; i++)
        {
            if (!((bool) ddeps[i]))
            {
                finished = false;
                break;
            }
        }
        //std::cout << "ddeps={ ";for(int32_t i=0; i<ndeps; i++)if(ddeps[i]) std::cout << "%i ",i);printf("} - finished=" << (finished?"true":"false") << "\n";

        if (!finished)
        {
            finished = cutcomp_slow(res2, ndeps, nodedeps, allres);
        }
        else
        {
            //std::cout << "found set={ ");for(int32_t h=0;h<num_nodes;h++)if(res2[h])std::cout << "" << h << " ";printf("}" << std::endl;
            std::set<uint64_t> tmpset;
            for (uint64_t i = 0; i < num_nodes; i++)
            {
                if (res2[i])
                {
                    tmpset.insert(i);
                }
            }
            allres.insert(tmpset);
            //std::cout << "set max_depth = " << maxdepth << "\n";
        }

        res2[n] = false;
    }
    return false;
}

void Tandem::cutcomp(compstatus *stat, std::vector<depsse> deps, std::set<std::set<uint64_t>> &allres)
{
    compstatus stat2;
    bool found;
    uint64_t tmp;
    if (stat == nullptr)
    {
        return;
    }
    for (uint64_t n = stat->cnode + 1; n < num_nodes; n++)
    {
        //std::cout << "checking node " << n << "\n";
        if (stat->cres & (1L << n))
        {
            continue;
        }
        if (n >= 2)
        {
            tmp = 3L << (n - 2);
        }
        if ((n >= 2) && ((stat->cres & tmp) == tmp))
        {
            continue;
        }
        depsse &newdep = deps[n];
        uint64_t depsok = 0;
        found = false;
        for (int32_t i = 0; i < MAXDEPS; i++)
        {
            stat2.cdep.dep[i] = stat->cdep.dep[i] | newdep.dep[i];
            //std::cout << "dep%02i: old=%08x curr=%08x  res=" << i,stat->cdep.dep[i],newdep.dep[i],stat2.cdep.dep[i] << "\n";
            if (stat2.cdep.dep[i] == 0xFFFFFFFF)
            {
                depsok++;
            }
            if (stat2.cdep.dep[i] != stat->cdep.dep[i])
            {
                found = true;
            }
        }
        if (!found)
        {
            continue;
        }
        //std::cout << "adding node %i; depsok=%i/" << n,depsok,MAXDEPS << "\n";
        stat2.cres = stat->cres | (1L << n);
        stat2.cnode = n;
        if (depsok < MAXDEPS)
        {
            cutcomp(&stat2, deps, allres);
        }
        else
        {
            std::set<uint64_t> tmpset;
            for (uint64_t i = 0; i < num_nodes; i++)
            {
                if (stat2.cres & (1L << i))
                {
                    tmpset.insert(i);
                }
            }
            //std::cout << "found set={ ");for(std::set<uint64_t>::iterator jt=tmpset.begin();jt!=tmpset.end();jt++)std::cout << "" << *jt << " ";printf("}" << std::endl;
            allres.insert(tmpset);
        }
    }
}

/* Cut the tandem into a vector of sub-tandems according to the specified list of cuts
 * Nodes in 'cutset' must appear in ascending order (sorted vector)
 * Returns a vector of sub-tandems, flows are initialized with original parameters
 */
std::vector<Tandem> Tandem::cut(CutNodeList &cutset)
{
    std::vector<Tandem> result;
    int32_t ncuts = cutset.size();
    if (ncuts == 0)
    {
        result.push_back(*this);
        return result;
    }
#ifdef TANDEM_DEBUG_CUTS
    std::cout << "Tandem::cut(): set contains the following cuts: { ";
for (int32_t i = 0; i < ncuts; i++)
std::cout << cutset[i] << " ";
std::cout << "}" << std::endl;
#endif
    // make the cuts
    uint64_t cnode = 0;
    for (int32_t i = 0; i < ncuts; i++)
    {
        Tandem t;
        t.Clear(); // this call is redundant...
        // fill-in the nodes array: determine the subset of nodes
        uint64_t cutnode = cutset[i]; // zero-based numbering
        if (cutnode > num_nodes)
        {
            continue;
        } // little safety check
        t.num_nodes = cutnode - cnode;
        t.nodes = std::vector<Node>(num_nodes);
#ifdef TANDEM_DEBUG_CUTS
        std::cout
<< "Sub-tandem " << i
<< ": nodes " << cnode << "-" << cutnode - 1
<< " (len=" << t.num_nodes << ")"
<< std::endl;
#endif
        // copy the subset of nodes
        for (uint64_t n = 0; n < t.num_nodes; n++)
        {
            t.nodes[n] = nodes[cnode + n];
        }
        // fill-in the flows array: split the flows affected
        t.flows = std::vector<Flow>(num_flows); // we'll have at most 'num_flows' flows in the sub-tandem
        // the sub-tandem spans from node 'cnode' (included) to 'cutset[i]' (excluded)
        t.num_flows = 0;
        for (uint64_t f = 0; f < num_flows; f++)
        {
            // handle the tagged flow assignement
            if (f == tagged_flow)
            {
                t.tagged_flow = t.num_flows;
            }
            // if the flow does not cross the current sub-tandem, skip it
            if ((flows[f].exit < cnode) || (flows[f].src >= cutnode))
            {
#ifdef TANDEM_DEBUG_CUTS
                std::cout << "  skip (" << flows[f].src << "," << flows[f].exit << ")" << std::endl;
#endif
                continue;
            }
            // this flow intersects the sub-tandem: trim it (i,j) --> (i,h) or (k,j)
            t.flows[t.num_flows] = flows[f];
            t.flows[t.num_flows].src = std::max(cnode, flows[f].src) - cnode;
            t.flows[t.num_flows].exit = std::min(cutnode - 1, flows[f].exit) - cnode;
            bool cut_out = flows[f].exit > cutnode - 1;
            bool cut_in = flows[f].src < cnode;
#ifdef TANDEM_DEBUG_CUTS
            std::cout
<< "  cut  (" << flows[f].src << "," << flows[f].exit << ") "
<< "as (" << t.flows[t.num_flows].src << "," << t.flows[t.num_flows].exit << ")  "
<< "cut_in=" << cut_in << " cut_out=" << cut_out << ""
<< std::endl;
#endif
            t.num_flows++;
        }
        t.ComputeNestingLevels();
        // store the new sub-tandem
        result.push_back(t);
        // move the current node pointer after the cut
        cnode = cutnode;
    }
    return result;
}

uint64_t Tandem::aggregate_flows(bool merge_tf)
{
    //std::cout << "Tandem::aggregate_flows(%s): %i flows, tagged_flow=" << (merge_tf?"true":"false"),num_flows,tagged_flow << "\n";
    for (uint64_t i = 0; i < num_flows - 1; i++)
    {
        // check if this flow can be aggregated with others
        for (uint64_t j = i + 1; j < num_flows; j++)
        {
            if ((flows[i].src == flows[j].src) && (flows[i].exit == flows[j].exit))
            {
                if (!merge_tf && ((j == tagged_flow) || (i == tagged_flow)))
                {
                    continue;
                }

                // aggregate j into i
                flows[i].burst += flows[j].burst;
                flows[i].rate += flows[j].rate;

                if (j == tagged_flow)
                {
                    tagged_flow = i;
                }
                else if (tagged_flow > j)
                {
                    tagged_flow--;
                }

                // shift down the flows array to eliminate j
                for (uint64_t k = j; k < num_flows - 1; k++)
                {
                    flows[k] = flows[k + 1];
                }

                // write junk into the left space, for safety...
                flows[i].uid = FLOW_INVALID_ID;
                flows[num_flows - 1].src = flows[num_flows - 1].exit = FLOW_INVALID_ID;
                num_flows--;

                j--; // compensate for shifted indexes
            }
        }
    }

    ComputeNestingLevels();

    // correct the nesting level of the tagged flow (when ar != 0)
    flows[tagged_flow].nesting_level = 1;

    //std::cout << "Tandem::aggregate_flows(): updated tagged_flow = %i, new num flows = " << tagged_flow,num_flows << "\n";
    return tagged_flow;
}

/*	Compute the output arrival curve for flow 'fid' at the exit of the current tandem
 * 	and store the parameters (sigma, rho) into 'fdest'
 */
bool Tandem::compute_output_arrival_curve(uint64_t fid, Flow *fdest, int32_t mode)
{
    if ((fid >= num_flows) || (fdest == nullptr))
    {
        std::cout
            << "Tandem::compute_output_arrival_curve(fid=" << fid << ", fdest=" << fdest
            << "): invalid parameter(s), bug?"
            << std::endl;

        return false;
    }
    Flow &f = flows[fid];
    if (f.exit != num_nodes - 1)
    {
        std::cout << "Tandem::compute_output_arrival_curve(fid=" << fid << "): flow does not exit the tandem, bug?"
                  << std::endl;
        return false;
    }

    double Dmin = 0.0;

    // check for special case: single-flow tandem
    if (num_flows == 1)
    {
        // do the easy math
        Dmin = LUDB_trivial();
    }
    else
    {
        // compute Dmin according to equation (13), section 4.1
        uint64_t nchildren = nesting_tree->num_children();
        if (nchildren == 0)
        {
            std::cout << "Tandem::compute_output_arrival_curve(fid=" << fid
                      << "): error, root node has no direct children\n";
            std::cout << "Affected tandem is:" << std::endl;
            Print();
            return false;
        }
        // compute the sum of LUDBs of the first-order child tnodes
        for (uint64_t i = 0; i < nchildren; i++)
        {
            uint64_t nsimpl;
            Dmin += LUDB_quiet(nesting_tree->children[i], mode, nsimpl);
        }
        // add the latency of C(1,N) of the possible leaves of the root node
        for (uint64_t i = 0; i < nesting_tree->num_leaves; i++)
        {
            uint64_t n = nesting_tree->leaves[i];
            Dmin += nodes[n].latency;
        }
    }

    //Compute output burst based on computed Dmin
    fdest->rate = f.rate;
    fdest->burst = f.burst + f.rate * Dmin;

    return true;
}

/* Compute all the necessary output arrival curves at the current tandem
 * Input: 't_next' is a reference to the subsequent tandem in the chain
 */
void Tandem::output_arrival_curves(Tandem &t_next, int32_t mode)
{
    for (uint64_t f = 0; f < num_flows; f++)
    {
        Tandem tt = *this; // a copy of the original sub-tandem
        Flow *tf = tt.getFlow(f);
        if (tf == nullptr)
        {
            continue;
        }
        // determine if the arrival curve of flow f is needed at the next tandem
        Flow *fnext = t_next.getFlowId(tf->uid);
        if (fnext == nullptr)
        {
            continue;
        }
        // more checks
        if (tf->exit != tt.NumNodes() - 1)
        {
            continue;
        } // flow exits midway, skip
        // this output curve is actually needed
        if (tf->src > 0)
        {
            // flow enters halfway, must split this sub-tandem recursively
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
            std::cout << "    !!!  flow=#%i uid=%i enters at node " << f, tf->uid, tf->src << ", must cut this sub-tandem there\n";
#endif
            // perform the cut
            CutNodeList cs;
            cs.push_back(tf->src);
            cs.push_back(tt.NumNodes());
            std::vector<Tandem> subt = tt.cut(cs);
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
            std::cout << "    splitted the subtandem into " << subt.size() << " sub-subtandems (should be 2)\n";
for (uint64_t i = 0; i < subt.size(); i++)
subt[i].Print();
#endif
            if (subt.size() != 2)
            {
                std::cout << "    wrong size of sub-splitted tandem" << std::endl;
                continue;
            }
            // recursively compute the OAC's required for the 2nd sub-tandem
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
            std::cout << "    recursion enter" << std::endl;
#endif
            subt[0].output_arrival_curves(subt[1], mode);
            // now proceed with the computation of the original OAC (2nd sub-tandem only is involved)
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
            std::cout << "    recursion exit" << std::endl;
#endif
            uint64_t tagged = subt[1].getFlowIdx(
                tf->uid
            ); // remap the tagged flow in original tandem into 2nd sub-subtandem
            if (tagged == FLOW_INVALID_ID)
            {
                std::cout << "BUG: unable to map tagged flow in sub-subtandem!" << std::endl;
                continue;
            }
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
            std::cout << "    mapped tagged flow (uid=%i) in sub-subtandem found at fid=" << tf->uid, tagged << "\n";
#endif
            subt[1].SetTaggedFlow(tagged);
            uint64_t tf2 = subt[1].aggregate_flows(false);
            subt[1].BuildNestingTree();
            bool ok = subt[1].compute_output_arrival_curve(tf2, fnext, mode);
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
            if (ok)
std::cout
<< "\t[OK] output arrival curve for tmp_tf="
<< tf2
<< " (uid=" << tf->uid
<< "): burst=" << std::fixed << std::setprecision(3) << fnext->burst
<< ", rate=" << std::fixed << std::setprecision(3) << fnext->rate
<< std::endl;
#endif
            continue;
        }
        tt.SetTaggedFlow(f);
        uint64_t tmptf = tt.aggregate_flows(false);
        tt.BuildNestingTree();
        // compute the output arrival curve
        bool res;
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
        std::cout << "    output_arrival_curve(): compute_output_arrival_curve(tmp_tf=" << tmptf << " uid=" << tf->uid << ")" << std::endl;
#endif
        res = tt.compute_output_arrival_curve(tmptf, fnext, mode);
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
        if (!res)
std::cout << "   compute_output_arrival_curve(tmp_tf=" << tmptf << ") returned error" << std::endl;
else
std::cout
<< "    [OK] output arrival curve for tmp_tf=" << tmptf << " (uid=" << tf->uid
<< "): burst=" << std::setprecision(3) << fnext->burst
<< ", rate=" << std::setprecision(3) << fnext->rate
<< std::endl;
#endif
    }
}

/// Handle LUDB computation for trivial cases.
/// Actually, only one case is handled: one flow crossing the whole tandem
/// \return The computed LUDB
double Tandem::LUDB_trivial()
{
    // case 1): one flow crossing the whole tandem
    if (num_flows == 1)
    {
        double latency = 0.0;
        double rate = TNODE_LUDB_INFINITY;
        for (uint64_t i = 0; i < num_nodes; i++)
        {
            latency += nodes[i].latency;
            rate = std::min(rate, nodes[i].rate);
        }

#ifdef TANDEM_DEBUG_LUDB
        std::cout << "Tandem::LUDB_trivial(): case 1) applied." << std::endl;
#endif
        cached_ludb = latency + flows[0].burst / rate;
        return cached_ludb;
    }

    // no supported cases detected, return error
    return -1000.0;
}

/* ************************************************************************************************************* */
/* ********************                   Lower-Bound related methods                       ******************** */
/* ************************************************************************************************************* */

// Auxiliary function which returns the n-th bit of the value x
bool inline Tandem::bit(uint64_t x, uint64_t n)
{
    uint64_t mask = (uint64_t) pow(2, n);
    return (x & mask) ? true : false;
}

// Returns true is Flow f traverses Node n
bool inline Tandem::isNodeCrossedByFlow(uint64_t n, uint64_t f)
{
    if (f >= num_flows)
    {
        return false;
    }
    return ((flows[f].src <= n) && (flows[f].exit >= n));
}

// Returns true is Flow f enters the tandem at Node n
bool inline Tandem::isNodeEnteredByFlow(uint64_t n, uint64_t f)
{
    if (f >= num_flows)
    {
        return false;
    }
    return (flows[f].src == n);
}

// Returns true is Flow f leaves the tandem at Node n
bool inline Tandem::isNodeLeftByFlow(uint64_t n, uint64_t f)
{
    if (f >= num_flows)
    {
        return false;
    }
    return (flows[f].exit == n);
}

//	Check if the given flow interferes with the tagged flow at the specified node
bool inline Tandem::isInterferingFlow(uint64_t n, uint64_t f, uint64_t tagged)
{
    uint64_t n0 = flows[tagged].src;
    if (!isNodeCrossedByFlow(n, f) || isNodeCrossedByFlow(n0, f))
    {
        return false;
    }
    // more conditions
    if (flows[f].exit < n0)
    {
        return false;
    }
    if (flows[f].src > flows[tagged].exit)
    {
        return false;
    }
    return true;
}

/*	This function performs FIFO de-multiplexing at the node exit point.
	The CDF_i contribution of the specified flow is devised from the overall CDF
	considering the CAF_i and CAF.
*/
bool Tandem::compute_CDF_i(uint64_t node_id, uint64_t flow_id)
{
    Curve caf, cdf, caf_i, cdf_i;
    double caf_latency = 0.0, cdf_latency = 0.0, caf_i_latency = 0.0;
    int32_t s_caf, s_caf_i, s_cdf, s_cdf_i;
    double x_caf, y_caf, x_caf_i, x_cdf, x_cdf_i;
    LinearSegment segcdf_i, segcaf, segcdf, segcaf_i;
    LinearSegment segcdf_next, segcaf_next;
    bool res;

    if (node_id >= num_nodes)
    {
#ifdef TANDEM_DEBUG_CDF_I
        std::cout << "Tandem::compute_CDF_i() invoked with invalid node ID." << std::endl;
#endif
        return false;
    }

    Node &node = nodes[node_id];

    caf_i = node.CAF_i[flow_id];
    caf = node.CAF;
    cdf = node.CDF;

#ifdef TANDEM_DEBUG_CDF_I
    std::cout << "Computing Node%i.CDF_" << node_id, flow_id << ":\n";
std::cout << "CAF-original = ";
caf.Print();
std::cout << "CAF_" << flow_id << "-original = ";
caf_i.Print();
std::cout << "CDF-original = ";
cdf.Print();
#endif

    int32_t nsegs_caf = caf.num_segments();
    int32_t nsegs_cdf = cdf.num_segments();
    int32_t nsegs_caf_i = caf_i.num_segments();
    int32_t nsegs_cdf_i = 0;

    if ((nsegs_caf == 0) || (nsegs_cdf == 0))
    {
        std::cout << "Tandem::compute_CDF_i(): empty CAF or CDF at node " << node_id << "\n";
        return false;
    }

    // CAF and CAF_i should have the same latency, so as CDF and CDF_i
    // Compared to them, CDF is retarded by the delay of the node's service curve (beta)
    cdf_latency = node.CDF.getLatency();
    caf_latency = node.CAF.getLatency();
    caf_i_latency = node.CAF_i[flow_id].getLatency();

    caf.removeLatency();
    if (caf_latency > 0.0)
    {
#ifdef TANDEM_DEBUG_CDF_I
        std::cout << "shifting CAF_i by " << -caf_latency << " horizontally\n";
#endif
        caf_i.shiftRight(-caf_latency);
    }
    cdf.removeLatency();

    // smooth CAF (remove colinear segments and recompute rates)
    caf.remove_colinear_segments();

    // update segments counts
    nsegs_caf = caf.num_segments();
    nsegs_cdf = cdf.num_segments();
    nsegs_caf_i = caf_i.num_segments();

#ifdef TANDEM_DEBUG_CDF_I
    std::cout << "Original latencies: CAF=%lf, CAF_i=%lf, CDF=" << caf_latency, caf_i_latency, cdf_latency << "\n";
std::cout << "Latencies removed, all curves now centered in origin" << std::endl;
std::cout << "CAF = ";
caf.Print();
std::cout << "CAF_" << flow_id << " = ";
caf_i.Print();
std::cout << "CDF = ";
cdf.Print();
#endif

    // pre-process CAF to add any missing y-breakpoints which are present in CDF instead
    res = false;
    for (int32_t i = 0; i < nsegs_cdf; i++)
    {
        segcdf = cdf.getSeg(i);
        if (caf.has_BPY_at(segcdf.y) < 0)
        {
#ifdef TANDEM_DEBUG_CDF_I
            std::cout << "CAF pre-processing: CDF has a BP at Y=%lf (X=" << segcdf.y, segcdf.x << ") which is not present in CAF, added\n";
#endif
            caf.add_BPY(segcdf.y);
            res = true;
        }
    }
    if (res)
    {
        nsegs_caf = caf.num_segments();
#ifdef TANDEM_DEBUG_CDF_I
        std::cout << "CAF-updated = ";
caf.Print();
#endif
    }

    // pre-process CAF to add any missing y-breakpoints which are present in CAF_i instead
    res = false;
    for (int32_t i = 0; i < nsegs_caf_i; i++)
    {
        segcaf_i = caf_i.getSeg(i);
        if (caf.has_BPX_at(segcaf_i.x) < 0)
        {
#ifdef TANDEM_DEBUG_CDF_I
            std::cout << "CAF pre-processing: CAF_i has a BP at X=%lf (Y=" << segcaf_i.x, segcaf_i.y << ") which is not present in CAF, added\n";
#endif
            caf.add_BPX(segcaf_i.x);
            res = true;
        }
    }
    if (res)
    {
        nsegs_caf = caf.num_segments();
#ifdef TANDEM_DEBUG_CDF_I
        std::cout << "CAF-updated = ";
caf.Print();
#endif
    }

    cdf_i.Zero();

    s_caf = 0;
    s_cdf = 0;
    x_caf = 0.0;
    y_caf = 0.0;
    x_cdf = 0.0;
    x_cdf_i = 0.0;
    double txbits, mytxbits;

    segcdf_i.Zero();
    while (s_caf < nsegs_caf)
    {
        // fetch current and next CAF segments
        segcaf = caf.getSeg(s_caf);
        if (s_caf < nsegs_caf - 1)
        {
            segcaf_next = caf.getSeg(s_caf + 1);
        }
        else
        {
            segcaf_next = segcaf;
        }
        // determine the number of bits transmitted in the current segment
        txbits = segcaf_next.y - segcaf.y;
        y_caf = segcaf.y + txbits / 2.0; // mid-point y in CAF segment of interest
        x_cdf = cdf.f_inv(y_caf, false);
        s_cdf = cdf.getSegmentDefining(x_cdf);
        segcdf = cdf.getSeg(s_cdf);

#ifdef TANDEM_DEBUG_CDF_I
        std::cout << "CAF segment: (%lf, %lf) <--> (%lf, %lf)   total bits transmitted = " << segcaf.x, segcaf.y, segcaf_next.x, segcaf_next.y, txbits << "\n";
#endif

        if ((Curve::IsZero(segcaf.rate)) || (txbits < Curve::Epsilon))
        {
            s_caf++;
            continue;
        }

        if (segcaf.id >= 0)
        {
            // someone's burst is being transmitted
            if (segcaf.id == flow_id)
            {
                // this flow is bursting, assign all CDF's rate to it
                segcdf_i.rate = segcdf.rate;
                mytxbits = txbits;
#ifdef TANDEM_DEBUG_CDF_I
                std::cout << "my burst: " << std::fixed << std::setprecision(3) << mytxbits << " bits transmitted by me\n";
#endif
            }
            else
            {
                // a different flow is bursting, this flow's rate is zero
                segcdf_i.rate = 0.0;
                mytxbits = 0.0;
#ifdef TANDEM_DEBUG_CDF_I
                std::cout << "alien burst (flow %02i), " << segcaf.id, mytxbits << " bits transmitted by me\n";
#endif
            }
            x_caf = segcaf.x;
            x_cdf = cdf.f_inv(segcaf_next.y, false);
        }
        else
        {
            // no burst is being transmitted
            x_caf = (segcaf_next.x + segcaf.x) / 2.0; // mid-point x in CAF interval of interest
            double caf_i_slope = caf_i.getSlope(x_caf, false);
            //x_cdf = cdf.f_inv(y_caf, false); // left or right limit???
            double cdf_slope = cdf.getSlope(x_cdf, false);
            if (segcaf.rate == 0.0)
            {
                // nobody was transmitting in the considered interval
                segcdf_i.rate = 0.0;
                mytxbits = 0.0;
                std::cout << "Tandem::compute_CDF_i(): should never get here, inconsistency detected!" << std::endl;
            }
            else
            {
                // compute the rate quota of this flow
                x_cdf = cdf.f_inv2(segcaf_next.y, false);
                segcdf_i.rate = cdf_slope * (caf_i_slope / segcaf.rate);
                mytxbits = segcdf_i.rate * (x_cdf - segcdf_i.x);
            }
#ifdef TANDEM_DEBUG_CDF_I
            std::cout << "Midpoint for CDF: Y=%lf   Midpoint for CAF/CAF_i: X=" << y_caf, x_caf << "\n";
std::cout << "affected CDF interval: %lf <--> %lf   slope = " << cdf.f_inv(segcaf.y, false), x_cdf, cdf_slope << " ";
std::cout << "  (segment: ";
segcdf.Print();
std::cout << ")" << std::endl;
std::cout << "affected CAF_i interval: %lf <--> %lf   slope = %lf   my bits transmitted = " << segcaf.x, segcaf_next.x, caf_i_slope, mytxbits << "\n";
#endif
        }
        //std::cout << "cdf_i.add_segment(): ";segcdf_i.Print();std::cout << std::endl;
        cdf_i.add_segment(segcdf_i);
        segcdf_i.y += mytxbits;
        segcdf_i.x = x_cdf;
        s_caf++;
    }
    segcdf_i.rate = 0.0;
    //std::cout << "add-last-segment: ";segcdf_i.Print();std::cout << std::endl;
    cdf_i.add_segment(segcdf_i);

#ifdef TANDEM_DEBUG_CDF_I
    std::cout << "CDF_" << flow_id << " = ";
cdf_i.Print();
#endif

    // sanity checks
    res = true;
#ifdef TANDEM_DEBUG_CHECKS
    if (std::fabs(cdf_i.getTotalTraffic() - caf_i.getTotalTraffic()) > LINEARSEGMENT_EPSILON)
    {
        std::cout << "Tandem::compute_CDF_" << flow_id
                  << "(): sanity check error: bits count mismatch at Node" << node_id
                  << " between CAF_" << flow_id
                  << " (" << caf_i.getTotalTraffic()
                  << ") and CDF_" << flow_id
                  << " (" << cdf_i.getTotalTraffic()
                  << "), bug somewhere!\n";
        std::cout << "difference: " << caf_i.getTotalTraffic() - cdf_i.getTotalTraffic() << "\n";
        res = false;
    }
#endif

    // finally re-add the CDF's latency to the CDF_i
    cdf_i.shiftRight(cdf_latency);
    // cleanup the resulting curve by removing colinear segments and recomputing rates (smoothing)
    cdf_i.remove_colinear_segments();

    // post-process the curve segments in order to reduce nasty rounding errors, etc.
#ifdef TANDEM_DEBUG_CHECKS
    cdf_i.SanityCheck();
#endif

#ifdef TANDEM_DEBUG_CDF_I
    std::cout << "Tandem::compute_CDF_" << flow_id << "() finished\n";
#endif
    node.CDF_i[flow_id] = cdf_i;
    return res;
}

/*	This method processes all nodes in the tandem (starting from the node where the tagged flow enters)
	and computes the delay experienced by the last bit of the tagged flow at the exit of the tandem
	Each interfering flow which enters the path along the way can transmit a burst at the beginning or at the end of the interval;
	the particular choice made for this run is encoded by the i-th bit of the input parameter
	Returns a negative number in case of errors.
*/
double Tandem::compute_delay(std::vector<bool> &combo)
{
    double delay = 0.0;
    double x1, x2;
    double tcaf, tcdf;
    bool res;

#ifdef TANDEM_DEBUG_LOWERBOUND
    std::cout << "Tandem::compute_delay(): combo ";
for (int32_t i = 0; i < num_flows; i++)
std::cout << "" << (int32_t)combo[i] << "";
std::cout << std::endl;
#endif

    // todo: this check may be unnecessary after conversion [] -> std::vector
    // if (combo == nullptr)
    // 	return TANDEM_LB_ERROR_INVALID_PARAM;

    // let's start from the node where the tagged flow enters
    uint64_t current_node = flows[tagged_flow].src;
    uint64_t last_node = flows[tagged_flow].exit;
    // the active interval for incoming flows is initially null
    x1 = 0.0;
    x2 = 0.0;

    while (current_node <= last_node)
    {
#ifdef TANDEM_DEBUG_LOWERBOUND
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Tandem::compute_delay(): processing Node " << current_node << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
#endif
        Node &node = nodes[current_node];

        // first, we must determine the interval (x1,x2) for the CAF_i of new flows entering the tandem at this node
        if (current_node > flows[tagged_flow].src)
        {
            double xx1 = Curve::PlusInfinity;
            double xx2 = 0.0;
            for (int32_t i = 0; i < num_flows; i++)
            {
                if (isNodeCrossedByFlow(current_node, i) && !isNodeEnteredByFlow(current_node, i))
                {
                    double x = nodes[current_node - 1].CDF_i[i].getLastBit(tagged_flow);
                    if (x > xx2)
                    {
                        xx2 = x;
                    }
                    x = nodes[current_node - 1].CDF_i[i].getFirstBit(tagged_flow);
                    if (x < xx1)
                    {
                        xx1 = x;
                    }
                }
            }
            x1 = nodes[current_node - 1].CDF.getFirstBit(tagged_flow);
            x2 = nodes[current_node - 1].CDF.getLastBit(tagged_flow);
            x1 = xx1;
            x2 = xx2;
#ifdef TANDEM_DEBUG_LOWERBOUND
            std::cout << "Node%i.CAF is active between x1=%lf and x2=" << current_node, xx1, xx2 << "\n";
std::cout << "New flows will transmit between x1=%lf and x2=" << x1, x2 << "\n";
#endif
        }

        // let's setup CAF and the various CAF_i in input at this node
        node.CAF.Zero();
        //Curve caftmp;
        for (int32_t i = 0; i < num_flows; i++)
        {
            bool entering = isNodeEnteredByFlow(current_node, i);
            bool crossing = isNodeCrossedByFlow(current_node, i);

            // check if this flow is entering the tandem at this node
            if (entering)
            {
                // create its brand-new CAF_i
                flows[i].create_caf(x1, x2, combo[i], i);
                //TODO: maybe we should check if this is exactly the tagged flow and do something special? (probably not)
                node.CAF_i[i] = flows[i].caf;
#ifdef TANDEM_DEBUG_LOWERBOUND
                std::cout << "Node%i.CAF" << current_node, i << " (new flow) = ";
node.CAF_i[i].Print();
#endif
            }
                // if not, check whether this flow has entered the tandem at an upstream node
            else if (crossing)
            {
                // copy its CAF_i from the CDF_i at the previous node
                if (current_node < 1)
                {
                    std::cout << "compute_delay(): SANITY CHECK 1 FAILED - CRASH!!!" << std::endl;
                    break;
                }
                node.CAF_i[i] = nodes[current_node - 1].CDF_i[i];
#ifdef TANDEM_DEBUG_LOWERBOUND
                std::cout << "Node%i.CAF%i = Node%i.CDF" << current_node, i, current_node - 1, i << " = ";
node.CAF_i[i].Print();
#endif
            }
                // this flow is not traversing the node at all, zero() the corresponding CAF_i entry
            else
            {
                node.CAF_i[i].Zero();
            }

            // sum the CAF_i just computed to the node's global CAF
            res = node.CAF.add(node.CAF_i[i]);
        }
#ifdef TANDEM_DEBUG_LOWERBOUND
        std::cout << "Node" << current_node << ".CAF = ";
node.CAF.Print();
std::cout << "Node" << current_node << ".Beta = ";
node.beta.Print();
#endif
        // now we can compute CDF = CAF (x) beta
        node.CDF = node.CAF;
        res = node.CDF.ConvolveWithLatencyRate(node.beta);

#ifdef TANDEM_DEBUG_CHECKS
        // perform some sanity checks about the computed CDF
        tcaf = node.CAF.getTotalTraffic();
        tcdf = node.CDF.getTotalTraffic();
        if (!Curve::AreEqual(tcaf, tcdf))
        {
            std::cout
                << "Tandem::compute_delay(): bitcount mismatch at Node" << current_node
                << " between CAF (" << tcaf
                << ") and CDF (" << tcdf
                << ") - bug in convolution?\n";
            std::cout << "Computed CDF is:\nNode" << current_node << ".CDF = ";
            node.CDF.Print();
            return TANDEM_LB_ERROR_INTERNAL;
        }
#endif

#ifdef TANDEM_DEBUG_LOWERBOUND
        std::cout << "Node" << current_node << ".CDF = ";
node.CDF.Print();
#endif

        // and finally we compute the various CDF_i... note that this step is not necessary at the last node!
        for (int32_t i = 0; i < num_flows; i++)
        {
            bool crossing = isNodeCrossedByFlow(current_node, i);
            bool leaving = isNodeLeftByFlow(current_node, i);

            // check if the flow has just crossed the node and it's not leaving the tandem here
            if (crossing && !leaving)
            {
                // we must compute its CDF_i
                res = compute_CDF_i(current_node, i);
                //res = compute_CDF_i_new(current_node, i);
#ifdef TANDEM_DEBUG_LOWERBOUND
                std::cout << "Node%i.CDF_" << current_node, i << " = ";
node.CDF_i[i].Print();
#endif
            }
                // this flow is leaving or it's not traversing the node at all: we don't need to compute the CDF_i
            else
            {
                node.CDF_i[i].Zero();
            }
        }

        // and that's it, we have finished to process this node
        current_node++;
    }

    delay = nodes[last_node].CDF.getLastBit(tagged_flow);

    return delay;
}

/*	This is the main function which performs the LowerBound algorithm on the tandem
*/
double
Tandem::LowerBound(bool quiet, uint64_t start_combo, uint64_t combo_range_len, float percentage, bool deterministic)
{
    uint64_t num_combos, eff_num_combos;
    double lowerbound = 0.0;
    clock_t t_start, t_end;
    double t_elapsed;
    RNG rng(1400);
    int32_t max_flows = 8 * sizeof(uint64_t);
    std::vector<bool> b_combo(num_flows);
    std::vector<uint64_t> lb_combos;

    if (!quiet)
    {
        std::cout << "Tandem::LowerBound() algorithm started." << std::endl;
        std::cout << "\n----------------------------------------------------" << std::endl;
    }

    //compute_output_rates(tagged_flow);

    if (num_flows > max_flows)
    {
        // replace uint64_t with an IndexSet to overcome the limitation
        std::cout << "ERROR: the current implementation does not support more than " << max_flows << " flows.\n";
        return -1.0;
    }

    // clear the curves at each node: beta, CAF, CAF_i, CDF and CDF_i
    for (uint64_t i = 0; i < num_nodes; i++)
    {
        nodes[i].clear_curves(num_flows);
    }

    num_combos = std::pow(2, num_flows) - 1; // we actually count <= 2^N - 1
    if (!quiet)
    {
        std::cout
            << "Number of combinations: " << num_combos + 1
            << "  (2 ^ " << num_flows << ")"
            << std::endl;
    }

    // to save time, we will only try "greedy" for flows which enter the tandem at node 0
    // in fact, "greedy" or "late" is the same at node 0 because flows only transmit the burst
    uint64_t combo_mask = 0ULL;
    for (int32_t i = 0; i < num_flows; i++)
    {
        if (isNodeEnteredByFlow(0, i))
        {
            combo_mask |= (1ULL << i);
        }
    }
    if (!quiet)
    {
        std::cout << "These flows transmit only the burst and will be considered greedy only:  ";
    }
    eff_num_combos = 1ULL;
    for (int32_t i = 0; i < num_flows; i++)
    {
        if (isNodeEnteredByFlow(0, i))
        {
            if (!quiet)
            {
                std::cout << std::setw(2) << std::setfill('0') << i << " ";
            }
            eff_num_combos <<= 1;
        }
    }
    eff_num_combos = (num_combos + 1) / eff_num_combos;

    //start_combo = 1536;  start_combo = 98344;
    //excel2 #32000 problemi cdf_i al nodo 6
    //nest8-20 #156 problemi convoluzione nodo 6
    if (combo_range_len != 0)
    {
        num_combos = start_combo + combo_range_len - 1;
    }

    uint64_t combo_base = start_combo;
    double combo_fstep = 100.0 / percentage;
    double combo_accum = 0.0, combo_istep = 0.0;
    uint64_t combo_step = (uint64_t) (combo_fstep < 1.01) ? 1 + floor(combo_fstep) : 1;
    uint64_t combo = combo_base;

    if (!quiet)
    {
        std::cout << "\nRandom combinations percentage: ";
        if (percentage >= 0.00001)
        {
            std::cout << std::setprecision(3) << percentage << " %" << std::endl;
        }
        else
        {
            std::cout << "" << percentage << " %%\n";
        }
        if (combo_fstep < 1.01)
        {
            std::cout << "Effective number of combinations that will be computed: " << eff_num_combos << std::endl;
        }
        else
        {
            eff_num_combos = (uint64_t) floor(num_combos * (percentage / 100.0));
            std::cout << "Approximate number of combinations that will be computed: " << eff_num_combos
                      << std::endl;
        }
    }
    eff_num_combos = 0;

    // start the benchmark timer
    t_start = clock();
    // start the lower-bound computations!
    while (combo_base <= num_combos)
    {
#ifdef TANDEM_LB_RANDOM_COMBOS
        if (combo_fstep > 1.01)
        {
            combo_accum += combo_fstep;
            double r = modf(combo_accum, &combo_istep);
            if (combo_istep < 1.0)
            {
                combo_step = 1;
            }
            else
            {
                combo_step = (uint64_t) combo_istep;
                combo_accum = r;
            }
            combo = combo_base + (uint64_t) floor(combo_step * rng.uniform(0.0, 1.0));
            if (combo > num_combos)
            {
                continue;
            }
        }
        else
        {
            // optimization: compute one combination only involving the flows which enter the tandem at node 0
            if (combo_base & combo_mask)
            {
                combo_base++;
                continue;
            }
            combo = combo_base;
        }
#else
        if (combo_base & combo_mask)
{
combo_base++;
continue;
}
combo = combo_base;
#endif
        // call the function which computes the delay for this combination and update the maximum
        double delay;
        for (int32_t i = 0; i < num_flows; i++)
        {
            b_combo[i] = bit(combo, i);
        }
#ifdef TANDEM_CATCH_EXCEPTIONS
        try
        {
            delay = compute_delay(b_combo);
        }
        catch (std::exception &ex)
        {
            std::cout << "Tandem::LowerBound(combo=" << combo << "): exception caught inside compute_delay(): "
                      << ex.what() << std::endl;
            delay = -1.0;
        }
#else
        delay = compute_delay(b_combo);
#endif
        if (delay < 0.0)
        {
            std::cout << "Tandem::LowerBound(combo=" << combo << "): computation failed." << std::endl;
            combo_base += combo_step;
            continue;
        }
        if (delay > lowerbound + LINEARSEGMENT_EPSILON)
        {
            lb_combos.clear();
            lb_combos.push_back(combo);
            lowerbound = delay;
            if (!quiet)
            {
                std::cout << "#" << std::setfill('0') << std::setw(8) << combo << ": new LowerBound = "
                          << std::fixed << std::setprecision(6) << lowerbound << std::endl;
            }
        }
        else if (std::fabs(delay - lowerbound) < LINEARSEGMENT_EPSILON)
        {
            lb_combos.push_back(combo);
        }
        //lowerbound = std::max(lowerbound, delay);
        //std::cout << "Combo %06i: delay computed = " << combo,delay,lowerbound << "   curr_lowerbound = %.3lf\n";
        combo_base += combo_step;
        eff_num_combos++;
    }
    t_end = clock();
    t_elapsed = ((double) t_end - (double) t_start) / (double) CLOCKS_PER_SEC;

    if (!quiet)
    {
        std::cout << "Lower Delay Bound = " << lowerbound << " obtained with the following combinations:\n";
        for (int32_t i = 0; i < lb_combos.size(); i++)
        {
            std::cout << "  " << std::setw(3) << i << ": #" << std::setw(8) << lb_combos[i] << " = ";
            for (int32_t j = 0; j < num_flows; j++)
            {
                std::cout << "" << (int32_t) bit(lb_combos[i], j) << "";
            }
            std::cout << std::endl;
        }
        std::cout << "Effective combinations computed: " << eff_num_combos << std::endl;
        if (!deterministic)
        {
            std::cout << "Computation time: " << std::fixed << std::setprecision(3) << t_elapsed << " seconds";
            if (t_elapsed > 0.1)
            {
                std::cout << "  (" << std::fixed << std::setprecision(2) << (double) (eff_num_combos) / t_elapsed
                          << " combos/s)";
            }
        }
        std::cout << std::endl;
    }

    return lowerbound;
}

/*	This function compute the max output rates of flows at each node
*/
bool Tandem::compute_output_rates(uint64_t tagged)
{
    int32_t node_start = flows[tagged].src;
    int32_t node_end = flows[tagged].exit;
    std::vector<double> rates(num_flows);
    double R;
    bool rate_reduction = false;

    if (tagged >= num_flows)
    {
        return false;
    }

    // handle the first node(s) separately
    for (int32_t i = 0; i <= node_start; i++)
    {
        nodes[i].set_output_rates(num_flows, {}, false);
    }
    nodes[node_start].set_output_rate(tagged, nodes[node_start].rate);

    // handle remaining nodes in the path of the tagged flow
    for (int32_t n = node_start + 1; n <= node_end; n++)
    {
#ifdef TANDEM_DEBUG_LBOPT
        std::cout << "Node " << n + 1 << ":";
#endif
        R = 0.0;
        for (int32_t f = 0; f < num_flows; f++)
        {
#ifdef TANDEM_DEBUG_LBOPT
            std::cout << " ";
flows[f].Print();
#endif
            if (f == tagged)
            {
                rates[f] = nodes[n - 1].get_output_rate(f);
                R += rates[f];
#ifdef TANDEM_DEBUG_LBOPT
                std::cout << "/tf=" << std::fixed << std::setprecision(3) << rates[f] << " ";
#endif
                continue;
            }
            if (!isInterferingFlow(n, f, tagged))
            {
                rates[f] = -1.0;
#ifdef TANDEM_DEBUG_LBOPT
                std::cout << "/NI ";
#endif
                continue;
            }
            if (isInterferingFlow(n - 1, f, tagged))
            {
                rates[f] = nodes[n - 1].get_output_rate(f);
                R += rates[f];
#ifdef TANDEM_DEBUG_LBOPT
                std::cout << "/r" << n, rates[f] << "=%.3lf ";
#endif
            }
            else if (isNodeEnteredByFlow(n, f))
            {
                rates[f] = flows[f].rate;
                R += rates[f];
#ifdef TANDEM_DEBUG_LBOPT
                std::cout << "/rho=" << std::fixed << std::setprecision(3) << rates[f] << " ";
#endif
            }
        }
#ifdef TANDEM_DEBUG_LBOPT
        std::cout << " - R=" << std::fixed << std::setprecision(3) << R << "";
#endif
        if (R > nodes[n].rate + LINEARSEGMENT_EPSILON)
        {
            // downscale rates
            double scale = nodes[n].rate / R;
            for (int32_t f = 0; f < num_flows; f++)
            {
                if (rates[f] >= 0.0)
                {
                    rates[f] *= scale;
                }
            }
            rate_reduction = true;
#ifdef TANDEM_DEBUG_LBOPT
            std::cout << ", scale by " << scale << "";
#endif
        }
        else
        {
            rate_reduction = false;
        }
#ifdef TANDEM_DEBUG_LBOPT
        std::cout << std::endl;
#endif
        nodes[n].set_output_rates(num_flows, rates, rate_reduction);
    }

    return true;
}

const std::byte RULE_BURST_EARLY{0x00};
const std::byte RULE_BURST_LATE{0x01};
const std::byte RULE_BURST_BOTH{0x02};
const std::byte RULE_ORDER_NONE{0x00};
const std::byte RULE_ORDER_BEFORE{0x01};
const std::byte RULE_ORDER_AFTER{0x02};
const std::byte RULE_ORDER_CLASH{0x03};

uint64_t Tandem::compute_lowerbound_rules(uint64_t tagged, std::vector<std::byte> &burst_time)
{
    int32_t node;
    std::vector<std::vector<std::byte> > order_rules(num_flows, std::vector<std::byte>(num_flows, RULE_ORDER_NONE));

    std::fill(burst_time.begin(), burst_time.end(), RULE_BURST_BOTH);
    // memset(burst_time, RULE_BURST_BOTH, num_flows);
    // memset(order_rules, RULE_ORDER_NONE, num_flows * num_flows);

    for (int32_t f = 0; f < num_flows; f++)
    {
#ifdef TANDEM_DEBUG_LBOPT
        std::cout << "Tagged flow = ";
flows[f].Print();
std::cout << std::endl;
#endif
        compute_output_rates(f);
#ifdef TANDEM_DEBUG_LBOPT
        std::cout << std::endl;
#endif
    }

    for (int32_t f = 0; f < num_flows; f++)
    {
        compute_output_rates(f);
        for (int32_t n = flows[f].src; n <= flows[f].exit; n++)
        {
            if (nodes[n].is_rate_reducing())
            {
                // all interfering flows which leave the tandem after node <n> must trasmit their burst before
                // the first bit of flow <f> enters the node
                for (int32_t l = 0; l < num_flows; l++)
                {
                    if ((l == f) || (flows[l].exit >= n))
                    {
                        continue;
                    }
                    for (int32_t m = flows[l].src; m <= flows[l].exit; m++)
                    {
                        if (isInterferingFlow(m, l, f))
                        {
                            // <l> transmits before <f>
                            order_rules[l][f] |= RULE_ORDER_BEFORE;
                        }
                    }
                }
            }
            else
            {
                // all interfering flows which enter the tandem at node <n> must transmit their burst at the end,
                // just before the last bit of flow <f> enters the node
                for (int32_t l = 0; l < num_flows; l++)
                {
                    if (l == f)
                    {
                        continue;
                    }
                    if (isNodeEnteredByFlow(n, l) && isInterferingFlow(n, l, f))
                    {
                        // <l> transmits at the end of <f>
                        order_rules[l][f] |= RULE_ORDER_AFTER;
                    }
                }
            }
        }
    }
    // now we have filled in the order_rules matrix
#ifdef TANDEM_DEBUG_LBOPT
    std::cout << "Dump of rules matrix:" << std::endl;
std::cout << "   ";
for (int32_t i = 0; i < num_flows; i++)
std::cout << "" << i << "|";
std::cout << std::endl;
for (int32_t i = 0; i < num_flows; i++)
{
std::cout << "" << i << "|";
for (int32_t j = 0; j < num_flows; j++)
{
switch (order_rules[i][j])
{
case RULE_ORDER_NONE:
std::cout << "-";
break;
case RULE_ORDER_BEFORE:
std::cout << "<";
break;
case RULE_ORDER_AFTER:
std::cout << ">";
break;
case RULE_ORDER_CLASH:
std::cout << "X";
break;
default:
break;
}
std::cout << " |";
}
std::cout << std::endl;
}
#endif
    // now determine the actual combinations
    uint64_t num_combos = 1;
    for (int32_t i = 0; i < num_flows; i++)
    {
        if (order_rules[i][tagged] == RULE_ORDER_BEFORE)
        {
            burst_time[i] = RULE_BURST_EARLY;
        }
        else if (order_rules[i][tagged] == RULE_ORDER_CLASH)
        {
            burst_time[i] = RULE_BURST_BOTH;
            num_combos *= 2;
        }
        else
        {
            burst_time[i] = RULE_BURST_LATE;
        }
    }
    // extra reduction of combinations for flows which enter the tandem at node zero
    for (int32_t i = 0; i < num_flows; i++)
    {
        if (burst_time[i] != RULE_BURST_BOTH)
        {
            continue;
        }
        if (isNodeEnteredByFlow(flows[tagged].src, i))
        {
            burst_time[i] = RULE_BURST_EARLY;
            num_combos /= 2;
        }
    }

    return num_combos;
}

/*	This method applies the exclusion criteria to determine
*/
double Tandem::LowerBound_experimental(float percentage)
{
    std::vector<std::byte> burst_time(num_flows);
    std::vector<bool> b_combo(num_flows);
    double lower_bound = 0.0;
    double delay = 0.0;
    uint64_t num_combos;
    uint64_t bitpos;
    RNG rng(1400);

    std::cout << "Tandem::LowerBound_experimental() algorithm started." << std::endl;
    std::cout << "\n----------------------------------------------------" << std::endl;

    // clear the curves at each node: beta, CAF, CAF_i, CDF and CDF_i
    for (int32_t i = 0; i < num_nodes; i++)
    {
        nodes[i].clear_curves(num_flows);
    }

    num_combos = compute_lowerbound_rules(tagged_flow, burst_time);
    for (int32_t i = 0; i < num_flows; i++)
    {
        b_combo[i] = (burst_time[i] == RULE_BURST_LATE);
    }

    std::cout << std::endl
              << "Effective number of combinations: " << num_combos + 1
              << ", total " << (uint64_t) (pow(2.0, num_flows))
              << " (2 ^ " << num_flows << " )"
              << std::endl;

    for (uint64_t combo = 0; combo < num_combos; combo++)
    {
#ifdef TANDEM_LB_RANDOM_COMBOS
        if (percentage < 99.999)
        {
            double p = rng.uniform(0.0, 100.0);
            if (p > percentage)
            {
                continue;
            }
        }
#endif

        // set b_combo
        bitpos = 0;
#ifdef TANDEM_DEBUG_LBOPT
        std::cout
<< "combo " << std::setw(6) << std::setfill('0') << combo
<< ": bitfield ";
#endif
        for (int32_t i = 0; i < num_flows; i++)
        {
            if (burst_time[i] == RULE_BURST_BOTH)
            {
                b_combo[i] = bit(combo, bitpos);
                bitpos++;
#ifdef TANDEM_DEBUG_LBOPT
                std::cout << " " << (int32_t)b_combo[i] << " ";
#endif
            }
#ifdef TANDEM_DEBUG_LBOPT
            else
std::cout << "(" << (int32_t)b_combo[i] << ")";
#endif
        }

        // compute lower bound for this combination
        delay = compute_delay(b_combo);
        if (delay < 0.0)
        {
            std::cout << "Combination " << combo << ": LowerBound computation failed." << std::endl;
            continue;
        }
        if (delay > lower_bound + LINEARSEGMENT_EPSILON)
        {
            lower_bound = delay;
        }
#ifdef TANDEM_DEBUG_LBOPT
        std::cout << "  -->  delay = %.3lf,  LowerBound = " << delay, lower_bound << "\n";
#endif
    }

    return lower_bound;
}

double Tandem::LowerBound_experimental_quiet(uint64_t &num_combos, float percentage)
{
    std::vector<std::byte> burst_time(num_flows);
    std::vector<bool> b_combo(num_flows);
    double lower_bound = 0.0;
    double delay = 0.0;
    uint64_t bitpos;
    RNG rng(1400);

    // clear the curves at each node: beta, CAF, CAF_i, CDF and CDF_i
    for (int32_t i = 0; i < num_nodes; i++)
    {
        nodes[i].clear_curves(num_flows);
    }

    num_combos = compute_lowerbound_rules(tagged_flow, burst_time);
    for (int32_t i = 0; i < num_flows; i++)
    {
        b_combo[i] = (burst_time[i] == RULE_BURST_LATE);
    }

    for (uint64_t combo = 0; combo < num_combos; combo++)
    {
#ifdef TANDEM_LB_RANDOM_COMBOS
        if (percentage < 99.999)
        {
            double p = rng.uniform(0.0, 100.0);
            if (p > percentage)
            {
                continue;
            }
        }
#endif

        // set b_combo
        bitpos = 0;
        for (int32_t i = 0; i < num_flows; i++)
        {
            if (burst_time[i] == RULE_BURST_BOTH)
            {
                b_combo[i] = bit(combo, bitpos);
                bitpos++;
            }
        }

        // compute lower bound for this combination
        delay = compute_delay(b_combo);
        if (delay < 0.0)
        {
            std::cout << "Combination " << combo << ": LowerBound computation failed." << std::endl;
            continue;
        }
        if (delay > lower_bound + LINEARSEGMENT_EPSILON)
        {
            lower_bound = delay;
        }
    }

    return lower_bound;
}

} // namespace deborah
