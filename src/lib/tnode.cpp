#include "tnode.h"
#include <iostream>
#include "simplex/simplex.h"
#include <math.h>
#include <time.h>
#include "rng.h"
#include <stdio.h>
#include <iomanip>

// Random Number Generator used in LUDB routines
static RNG rng_ludb(9981);

namespace deborah
{

TNode::TNode()
{
    num_tnodes = 0;
    num_leaves = 0;
    children = {};
    leaves = nullptr;
    flow = 0;
    pflow = nullptr;
    pi = nullptr;
    pi_eq = nullptr;
}

TNode::TNode(uint64_t flow_id, uint64_t nt, uint64_t nl)
{
    flow = flow_id;
    num_tnodes = nt;
    num_leaves = nl;
    if (num_tnodes > 0)
    {
        children = std::vector<std::shared_ptr<TNode>>(num_tnodes);
    }
    else
    {
        children = {};
    }
    if (num_leaves > 0)
    {
        leaves = new uint64_t[num_leaves];
    }
    else
    {
        leaves = nullptr;
    }
    pflow = nullptr;
    pi = nullptr;
    index_set.clear();
    index_schema.clear();
    index_nodes.clear();
}

TNode::~TNode()
{
    //Warning: don't use pflow, since the pointed Flow object could have been already destroyed

    //delete leaves array
    if (leaves != nullptr)
    {
        delete[] leaves;
        leaves = nullptr;
    }

    if (pi != nullptr)
    {
        delete pi;
        pi = nullptr;
    }

    if (pi_eq != nullptr)
    {
        delete pi_eq;
        pi_eq = nullptr;
    }

    flow = 0;
    pflow = nullptr;
}

void TNode::Print(bool subtree)
{
    if (pflow == nullptr)
    {
        std::cout
            << "TNode: flow id #"
            << std::setw(2) << std::setfill('0') << flow
            << " has NULL pointer to Flow, internal error!"
            << std::endl;
        return;
    }

    std::cout
        << "TNode: flow (" << pflow->src + 1 << "," << pflow->exit + 1 << ")"
        << "   flow_id=#" << std::setw(2) << std::setfill('0') << flow
        << std::endl;

    //Print leaves
    std::cout << "       leaf nodes:";
    if (num_leaves == 0)
    {
        std::cout << " none";
    }
    else
    {
        for (int32_t i = 0; i < num_leaves; i++)
        {
            std::cout << " " << std::setw(2) << std::setfill('0') << leaves[i] + 1;
        }
    }
    std::cout << std::endl;

    //Print directly nested flows: set S(h,k)
    std::cout << "       child flows:";
    if (num_tnodes == 0)
    {
        std::cout << " none";
    }
    else
    {
        for (int32_t i = 0; i < num_tnodes; i++)
        {
            auto src = children[i]->pflow->src + 1;
            auto exit = children[i]->pflow->exit + 1;
            std::cout << " (" << src << "," << exit << ")";
        }
    }
    std::cout << std::endl;

    //Print PI_C service curve associated to this node
    std::cout << "       PI_C: ";
    if (num_leaves == 0)
    {
        std::cout << "none" << std::endl;
    }
    else
    {
        std::cout
            << "latency=" << std::fixed << std::setprecision(2) << pi_c.x
            << " rate=" << std::fixed << std::setprecision(2) << pi_c.rate
            << std::endl;
    }

    //std::cout << "       PI: "); pi->Print();

    //Print set of indexes available at this node
    std::cout << "       Indexes: [";
    uint64_t idx_pos = 0;
    uint32_t idx_combos = 1;    //This overflows, but current tests expect it to ¯\_(ツ)_/¯
    for (int32_t i = 0; i < index_schema.size(); i++)
    {
        uint64_t setsize = index_schema[i];
        std::cout << " [";
        for (int32_t j = idx_pos; j < idx_pos + setsize; j++)
        {
            std::cout << " " << index_set[j];
            idx_combos *= index_set[j];
        }
        std::cout << " ]";
        idx_pos += setsize;
    }
    std::cout << " ] --> combinations = " << idx_combos << std::endl;

    fflush(stdout);

    //Recurse printing the subtree, if requested
    if (subtree)
    {
        for (int32_t i = 0; i < num_tnodes; i++)
        {
            children[i]->Print(true);
        }
    }
}

void TNode::PrintIndexSet(const IndexSet &set, bool compact)
{
    std::cout << "[";
    uint64_t idx_pos = 0;
    if (compact)
    {
        for (int32_t i = 0; i < index_schema.size(); i++)
        {
            uint64_t setsize = index_schema[i];
            std::cout << " [";
            for (int32_t j = idx_pos; j < idx_pos + setsize; j++)
            {
                std::cout << " " << set[j];
            }
            std::cout << " ]";
            idx_pos += setsize;
        }
    }
    else
    {
        for (int32_t i = 0; i < index_schema.size(); i++)
        {
            uint64_t setsize = index_schema[i];
            std::cout << " [";
            for (int32_t j = idx_pos; j < idx_pos + setsize; j++)
            {
                std::cout << " " << set[j] << "/" << index_set[j];
            }
            std::cout << " ]";
            idx_pos += setsize;
        }
    }
    std::cout << " ]";
}

bool TNode::has_leaf()
{
    if (num_leaves > 0)
    {
        return true;
    }
    return false;
}

uint64_t TNode::num_children()
{
    return num_tnodes;
}

/*	This method computes the PI() pseudo-affine service curve
	for the current t-node.
	1. Retrieve the PI curves of all the child flows
	2. Compute their equivalent service curves (--> at the child node!)
	3. Convolve them all, convolve also with the local PI and return
*/
ParametrizedPseudoAffine &TNode::getPI(const IndexSet &idx_set, std::vector<simplex::Constraint> &c_set)
{
#ifdef TNODE_VERBOSE_LUDB
    std::cout
        << ">>>>>> TNode #" << std::setw(2) << std::setfill('0') << flow
        << ": getPI() --> IndexSet = ";
    PrintIndexSet(idx_set, false);
    std::cout << std::endl;
#endif
    //if this is a leaf flow, return the node's PI directly

    if (num_tnodes == 0)
    {
#ifdef TNODE_VERBOSE_LUDB
        std::cout
            << "<<<<<< TNode #" << std::setw(2) << std::setfill('0') << flow
            << ": getPI() finished (no children)"
            << std::endl;
#endif
        return *pi;
    }

    // initialize the PI with the local pi_c curve
    if (has_leaf())
    {
        pi->MakeLatencyRate(pi_c);
    }
    else
    {
        pi->Zero();
    }

    // now convolve it with the equivalent service curve of each child flow
    for (int32_t i = 0; i < num_tnodes; i++)
    {
        // get index set for this child
        IndexSet child_set;
        if (getSubset(idx_set, i, child_set) == false)
        {
            continue;
        }
        // get PI curve of the child flow
        ParametrizedPseudoAffine &child_pi = children[i]->getPI(child_set, c_set);
        // compute equivalent service curve at the child node
        ParametrizedPseudoAffine &child_eq = children[i]->computeEquivalent(child_pi, child_set, c_set);
        // convolve the equivalent service curve into the local PI curve
        pi->Convolve(child_eq);
    }

#ifdef TNODE_VERBOSE_LUDB
    std::cout
        << "<<<<<< TNode #"  << std::setw(2) << std::setfill('0') << flow
        << ": getPI() finished"
        << std::endl;
#endif
    return *pi;
}

/*	This method computes an equivalent service queue according to corollary 2.7
	Input parameters:
		- the pseudo-affine curve that will be combined with the
	  	flow corresponding to the current node in order to determine the equivalent
		service curve experienced by the flow
		- the index set for the current node, specifying the index value selected for the max() expression
		- the constraint set, which will be filled in accordingly
*/
ParametrizedPseudoAffine &TNode::computeEquivalent(ParametrizedPseudoAffine &p, const IndexSet &set, std::vector<simplex::Constraint> &c_set)
{
    if (pi_eq == nullptr)
    {
        std::cout
            << "TNode::computeEquivalent(tnode flow #" << std::setw(2) << std::setfill('0') << flow
            << "): pi_eq == NULL!!!"
            << std::endl;
        return p;
    }

#ifdef TNODE_VERBOSE_LUDB
    std::cout
        << ">>> TNode #" << std::setw(2) << std::setfill('0') << flow
        << ": computeEquivalent(";
    PrintIndexSet(set, false);
    std::cout << ")" << std::endl;
#endif

    pi_eq->Zero();

    //determine which is the ith element of the max() expression that we have been asked for
    uint64_t max_index = set[0]; // max_index is the index of the chosen 'ith' element
    uint64_t max_index_limit = index_set[0];
    //range-check (for bugs)
    if (max_index >= max_index_limit)
    {
        std::cout
            << "!!! BUG !!! Requested Index=" << max_index
            << " but highest available is " << max_index_limit
            << std::endl;
        return *pi_eq;
    }

    /*****  now compute the value of the chosen element of the max() expression  *****/

    uint64_t nvars = p.delay.getNumVariables();
    Polynomial max_term(nvars);     // max_term = Zero
    simplex::Constraint diseq(nvars, flow); // will store the constraints in the form "ax + b >= 0"
    Polynomial ptemp(nvars);

    // an index value exceeding the number of stages of the PI corresponds to the case max()={0}
    if (max_index < p.stages.size())
    {
        // this is the case where max(0, X1,..Xn) == Xi
        max_term.SetConstant(pflow->burst); //max_term = sigma (of the flow represented by the current t-node)
        // compute the term [ (sigma - sigma_i) / rho_i ]+  (where i == max_index)
        max_term = (max_term - p.stages[max_index].sigma) / p.stages[max_index].rho;
#ifdef TNODE_VERBOSE_LUDB
        std::cout << "MAX_TERM = ";
        max_term.Print();
        std::cout << std::endl;
#endif
        //add the constraint (sigma - sigma_i) > 0 to implement the + apice, i.e. max(0, max_term)
        diseq = max_term;
        // store the constraint only if it's not a tautology
        if (!diseq.isTautology())
        {
            c_set.push_back(diseq);
        }

        // now add the constraints (sigma - sigma_i)/rho_i >= all the others stages of the local PI
        for (int32_t i = 0; i < p.stages.size(); i++)
        {
            if (i == max_index)
            {
                continue;
            }
            ptemp.Zero();
            ptemp.SetConstant(pflow->burst);
            ptemp = (ptemp - p.stages[i].sigma) / p.stages[i].rho;
            // new constraint: max_term >= ptemp
            diseq = max_term - ptemp;
            // store the constraint only if it's not a tautology
            if (!diseq.isTautology())
            {
                c_set.push_back(diseq);
            }
        }
    }
    else
    {
        // this is the case where max(0, X1,..Xn) == 0
#ifdef TNODE_VERBOSE_LUDB
        std::cout << "               case MAX()=0 hit" << std::endl;
#endif
        max_term.Zero();
        // add the constraints sigma - sigma_i < 0 for all the stages of the local PI
        for (int32_t i = 0; i < p.stages.size(); i++)
        {
            ptemp.Zero();
            ptemp.SetConstant(pflow->burst);
            ptemp = (ptemp - p.stages[i].sigma) / p.stages[i].rho;
            diseq = max_term - ptemp;
            // store the constraint only if it's not a tautology
            if (!diseq.isTautology())
            {
                c_set.push_back(diseq);
            }
        }
    }

    /*****  now max_term contains the chosen element of the max(0, X1,..,Xn) term  *****/

    // now make max_term = S(i,j) + (sigma - sigma_i) / rho_i
    max_term.SumVariableCoefficient(flow, 1.0);

    // compute the delay term of the equivalent service curve: D = PI.D + max() + S(i,j)
    pi_eq->delay = p.delay + max_term; // D = PI.D + max() +S(i,j)

    // compute the leaky-bucket stages of the equivalent service curve
    ParametrizedLeakyBucket newstage(nvars);
    for (int32_t i = 0; i < p.stages.size(); i++) //the number of stages is the same as the PI
    {
        // new_sigma = rho_x * [max() + S(i,j)] - (sigma - sigma_x)
        newstage.sigma = (max_term * p.stages[i].rho) - pflow->burst + p.stages[i].sigma;
        // new_rho = rho_x - rho
        newstage.rho = p.stages[i].rho - pflow->rate;
        // store computed stage
        pi_eq->AddStage(newstage);
    }

#ifdef TNODE_VERBOSE_LUDB
    std::cout
        << "<<< TNode #" << std::setw(2) << std::setfill('0') << flow
        << ": Equivalent service curve: ";
    pi_eq->Print();
#endif
    return *pi_eq;
}

/*	This method computes the equivalent service curve
	for the current t-node.
	1. Retrieve all the equivalent service curves of the child nodes
	2. Convolve them all and convolve also with the node's PI_C, obtaining the local PI
	3. Compute the equivalent service curve of the local PI and return it
*/
ParametrizedPseudoAffine TNode::computeEquivalentServiceCurve(const IndexSet &idx_set)
{
    //ignored
    std::vector<simplex::Constraint> constraints;

    // copied from getDelayBound
    // todo: what is a PI? why is it needed?
    auto pi = getPI(idx_set, constraints);
    auto eq = computeEquivalent(pi, idx_set, constraints);

    return eq;
}

/*	This method returns a reference to a Polynomial which contains the delay bound
	computed at the current t-node for the given set of indexes
	Input:
		idx_set is an instance of the index set corresponding to the desired combination
		constraints will be filled in with the constraints resulting from the desired combination
	Output:
		a Polynomial reference to the delay bound expression: h(PI,a)
*/
simplex::Objective& TNode::getDelayBound(const IndexSet &idx_set, std::vector<simplex::Constraint> &constraints)
{
    ParametrizedPseudoAffine &p = getPI(idx_set, constraints);
    ParametrizedPseudoAffine &pe = computeEquivalent(p, idx_set, constraints);
    //subtract the S(h,k) term at the root node
    pe.delay.SumVariableCoefficient(flow, -1.0);
    return pe.delay;
}

// Convert a number into the corresponding IndexSet combination
bool TNode::getCombo(IndexSet &set, uint64_t n)
{
    uint64_t q, r = n;
    set.clear();
    int32_t size = index_set.size();
    for (int32_t i = size - 1; i >= 0; i--)
    {
        q = r % index_set[i];
        set.insert(set.begin(), q);
        r = r / index_set[i];
    }
    return true;
}

/*	Increments the IndexSet and return the position of the highest digit that has been incremented
	'set' is the IndexSet to be incremented
	'index' is the first index position to be incremented (higher positions will be unchanged)
	If 'index' is negative, the increment will be applied to the whole set of indexes normally
	The method returns the position of the highest index incremented, or -1 if a wrap-around has occurred
*/
int32_t TNode::incCombo(IndexSet &set, int32_t index)
{
    int32_t size = index_set.size();
    // if we are not starting from the bottom, zero the last index positions
    if (index < 0)
    {
        index = size - 1;
    }
    else
    {
        for (int32_t i = index + 1; i < size; i++)
        {
            set[i] = 0;
        }
    }
    // perform the increment
    while (index >= 0)
    {
        if ((++set[index]) < index_set[index])
        {
            return index;
        }
        // carry set: increment also the next index
        set[index] = 0;
        index--;
    }
    // signal wrap-around
    return -1;
}

// Returns the total number of combinations corresponding to the given set of indexes
uint64_t TNode::getNumCombos()
{
    uint64_t num = 1;
    for (int32_t i = 0; i < index_set.size(); i++)
    {
        num *= index_set[i];
    }
    return num;
}

/* Returns the index subset corresponding to the Ith child of this node
 * Parameters: n_child = number of set desired (0..num_tnodes-1)
 */
bool TNode::getSubset(IndexSet &subset, int32_t n_child)
{
    return getSubset(index_set, n_child, subset);
}

bool TNode::getSubset(const IndexSet &set, int32_t n_child, IndexSet &subset)
{
    subset.clear();
    uint64_t count, size;

    if ((n_child < 0) || (n_child >= num_tnodes))
    {
        return false;
    }

    // entry #0 in index_schema is the index for the current node
    count = 0;
    // now sum the number of indexes of the preceding children
    for (int32_t i = 0; i <= n_child; i++)
    {
        count += index_schema[i];
    }

    // fetch the number of indexes for the child requested
    size = index_schema[n_child + 1];
    // copy the subset of indexes of the child into the new set
    for (int32_t i = 0; i < size; i++)
    {
        subset.push_back(set[count + i]);
    }

    return true;
}

//todo: distinguish this two LUDB methods

/// Computes the LUDB at the current TNode level, rather than at the Tandem level.
/// \param delay will contain the expression of the objective function of the "optimum simplex"
/// \param solution will contain the solution which gives the optimum
/// \param nvars is the number of variables at the top level (# of flows)
/// \param good_combos will contain the IndexSets that have produces good simplexes at the lower levels.
/// Can be NULL, in this case information about good combinations will be discarded.
/// \return
double TNode::computeLUDB(simplex::Objective &delay, simplex::Solution &solution, int32_t nvars, IndexSetVector *good_combos)
{
    uint64_t num_combos = getNumCombos();
    IndexSet iset;
    int32_t iresult;
    double MinDelayBound = TNODE_LUDB_INFINITY;

    simplex::Solution best_optimum(nvars);
    simplex::Objective best_delay(nvars);

    for (uint64_t combo = 0; combo < num_combos; combo++)
    {
        // get the combination of indexes corresponding to this iteration
        getCombo(iset, combo);

        // compute the delay bound expression and constraints at the root t-node for this iteration
        std::vector<simplex::Constraint> constraints;
        simplex::Objective &delay = getDelayBound(iset, constraints);

        // compute the simplex
        simplex::Solution optimum(delay.getNumVariables());

        if (constraints.size() == 0 && delay.getNumVariables() == 1)
        {
            optimum.setConstant(delay.getConstant());
            iresult = SIMPLEX_RESULT_OK;
            //std::cout << "TNode::computeLUDB(): trivial simplex skipped.\n");
        }
        else
        {
            iresult = simplex::minimum(delay, constraints, flow, optimum);
        }

        if (iresult != SIMPLEX_RESULT_OK)
        {
            continue;
        }

        // update the least delay bound, if necessary
        if (optimum.getConstant() < MinDelayBound)
        {
            MinDelayBound = optimum.getConstant();
            best_optimum = optimum;
            best_delay = delay;
        }
        if (good_combos != nullptr)
        {
            good_combos->push_back(iset);
        }
    }

    solution = best_optimum;
    delay = best_delay;

    return MinDelayBound;
}

/// Computes the LUDB at the current t-node by recursively computing it at the children t-nodes and trying
/// only the combinations reported by them.
/// \param delay will contain the expression of the objective function of the "optimum simplex"
/// \param solution will contain the solution which gives the optimum
/// \param nvars is the number of variables at the top level (# of flows)
/// \param good_combos will contain the IndexSets that have produces good simplexes at the lower levels
/// \param max_good_combos
/// \param eta
/// \return the number of simplexes computed in this sub-tree
uint64_t TNode::computeLUDB(
    simplex::Objective &delay, simplex::Solution &solution, int32_t nvars, IndexSetVector &good_combos,
    int32_t max_good_combos, bool eta
)
{
    simplex::Objective d(nvars);
    simplex::Solution s(nvars);
    simplex::Solution sol(nvars);

    simplex::Solution optimum(nvars);
    simplex::Objective delay_opt(nvars);

    uint64_t num_combos = getNumCombos();
    IndexSet set;
    std::vector<simplex::Constraint> constraints;
    double MinDelayBound = TNODE_LUDB_INFINITY;
    uint64_t num_simplexes = 0;
    uint64_t num_good_combos = 0;

    std::vector<IndexSetVector> children_combos(num_tnodes);
    std::vector<int32_t> children_num_combos(num_tnodes);
    std::vector<uint64_t> children_combo(num_tnodes);

    clock_t t_start = clock();

    // terminate recursion when a t-node has no children
    if (num_tnodes == 0)
    {
        //std::cout << "TNode::computeLUDB_fast(t-node=%i): no child t-nodes, compute local simplex\n",flow);
        computeLUDB(delay, solution, nvars, &good_combos);
        return 1;
    }

    // compute LUDB at all child nodes
    // each child tells the father the set of good combinations to be tried on him
    // the father will try all the indexes obtained by all the possible juxtapositions of the children's good combos,
    // plus its own index. In turn, the father will tell its own grandfather only the combinations that were feasible
#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
    std::cout
        << "TNode::computeLUDB_experimental(t-node=" << flow << "):"
        << " computing simplexes at " << num_tnodes << " child t-nodes,"
        << " max_good_combos = " << max_good_combos
        << std::endl;
#endif
    num_combos = index_set[0];
    for (int32_t i = 0; i < num_tnodes; i++)
    {
        d.Zero();
        s.Zero();
        num_simplexes += children[i]->computeLUDB(
            d, s, nvars, children_combos[i], max_good_combos,
            eta
        );
        children_num_combos[i] = children_combos[i].size();
        num_combos *= children_num_combos[i];
#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
        std::cout
            << "TNode" << std::setw(2) << std::setfill('0') << flow
            << ".child" << std::setw(2) << std::setfill('0') << i
            << " returned " << children_num_combos[i] << " good combo(s):"
            << std::endl;

        for (int32_t j = 0; j < children_num_combos[i]; j++)
        {
            std::cout
                << "\t" << std::setw(2) << std::setfill('0') << j
                <<":  [ ";

            int32_t size = children_combos[i][j].size();
            for (int32_t k = 0; k < size; k++)
                std::cout << children_combos[i][j][k] << " ";
            std::cout << "]\n";
        }
#endif
    }

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
    std::cout
        << "TNode::computeLUDB_experimental(t-node=" << flow
        << "): computing " << num_combos << " simplexes:"
        << std::endl;
#endif
    if (num_combos > 1000000)
    {
        eta = true;
    }
    for (uint64_t combo = 0; combo < num_combos; combo++)
    {
        // check if it's time to update ETA
        if (eta && ((combo & 0x1FFF) == 0x1000))
        {
            clock_t t_now = clock();
            double t_elapsed = ((double) t_now - (double) t_start) / (double) CLOCKS_PER_SEC;
            std::cout
                << std::setw(6) << std::setfill('0') << combo << "/" << std::setw(6) << std::setfill('0')
                << num_combos
                << ": ETA = " << ((num_combos - combo) * t_elapsed) / combo << " s "
                << std::setprecision(6) << "LUDB = " << MinDelayBound
                << std::endl;
        }

        // compute indexes for children
        uint64_t q, r;
        r = combo;
        for (int32_t i = num_tnodes - 1; i >= 0; i--)
        {
            q = r % children_num_combos[i];
            children_combo[i] = q;
            r = r / children_num_combos[i];
        }
        q = r % index_set[0];

        // assemble top-level combination
        set.clear();
        set.push_back(q);
        for (int32_t i = 0; i < num_tnodes; i++)
        {
            auto v = children_combo[i];
            IndexSet cs = (children_combos[i])[v];
            int32_t size = cs.size();
            for (int32_t j = 0; j < size; j++)
            {
                set.push_back(cs[j]);
            }
        }

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
        std::cout << "TNode " << std::setw(2) << std::setfill('0') << flow << ":  ";
        PrintIndexSet(set, true);
#endif

        // compute top-level simplex
        constraints.clear();
        simplex::Objective &del = getDelayBound(set, constraints);
        int32_t iresult = simplex::minimum(del, constraints, flow, sol);
#ifdef TNODE_DEBUG_LUDB_CONSTRAINTS
        std::cout << " --> " << constraints.size() << " constraints:" << std::endl;
        for (int32_t con = 0; con < constraints.size(); con++)
        {
            std::cout << "C" << std::setw(3) << std::setfill('0') << con << ": ";
            constraints[con].Print();
            std::cout << std::endl;
        }
        std::cout << "Obj.Fun: ";
        del.Print();
        std::cout << std::endl;
        std::cout << "Result";
#endif
        if (iresult != SIMPLEX_RESULT_OK)
        {
            // simplex was unfeasible
#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
            std::cout << " --> N/F" << std::endl;
#endif
            continue;
        }
        if (max_good_combos != LUDB_MODE_EXACT)
        {
            // heuristic mode
            if (sol.getConstant() < (MinDelayBound - TNODE_LUDB_EPSILON))
            {
                MinDelayBound = sol.getConstant();
                optimum = sol;
                delay_opt = del;
                good_combos.clear();
                good_combos.push_back(set);
                num_good_combos = 1;
#ifndef TNODE_DEBUG_LUDB_EXPERIMENTAL
                //if(eta) std::cout << "TNode %02i:  ",flow); PrintIndexSet(set,true); std::cout << " --> new LUDB = %lf\n",MinDelayBound);
#endif
            }
            else if (std::fabs(sol.getConstant() - MinDelayBound) <= TNODE_LUDB_EPSILON)
            {
                good_combos.push_back(set);
                num_good_combos++;
            }
        }
        else
        {
            // exact mode
            if (sol.getConstant() < MinDelayBound)
            {
                MinDelayBound = sol.getConstant();
                optimum = sol;
                delay_opt = del;
#ifndef TNODE_DEBUG_LUDB_EXPERIMENTAL
                //if(eta) std::cout << "TNode %02i:  ",flow); PrintIndexSet(set,true); std::cout << " --> new LUDB = %lf\n",MinDelayBound);
#endif
            }
            good_combos.push_back(set);
        }

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
        std::cout << " --> " << sol.GetConstant() << std::endl;
#endif
    }

    solution = optimum;
    delay = delay_opt;
    num_simplexes += num_combos;

    // in heuristic mode, trim the set of good combinations to the desired size
    if (max_good_combos != LUDB_MODE_EXACT)
    {
        // delete exceeding combinations randomly
        while (num_good_combos > max_good_combos)
        {
            uint64_t rnd_combo = floor(rng_ludb.uniform(0.0, num_good_combos));
            if (rnd_combo >= num_good_combos)
            {
                std::cout << "TNode::computeLUDB(): RNG right limit trespassed!" << std::endl;
                rnd_combo = num_good_combos - 1;
            }
            good_combos.erase(good_combos.begin() + rnd_combo);
            num_good_combos--;
        }
    }

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
    std::cout << "TNode::computeLUDB_experimental(): finished, " << num_simplexes << " simplexes computed in this sub-tree" << std::endl;
#endif
    return num_simplexes;
}

} // namespace deborah
