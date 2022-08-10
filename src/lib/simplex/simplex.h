#include "../tnode.h"

#include "objective.h"
#include "constraint.h"
#include "solution.h"

namespace simplex
{

#define SIMPLEX_RESULT_OK             0
#define SIMPLEX_RESULT_UNFEASIBLE    -1
#define SIMPLEX_RESULT_UNBOUNDED    -2
#define SIMPLEX_RESULT_UNK_ERROR    -999

#define SIMPLEX_ERROR_OK             0
#define SIMPLEX_ERROR_MEMORY        -1


//MMAX = Max number of constraints
#define  MMAX  320

//NMAX = Max number of variables
#define  NMAX  63

[[maybe_unused]] extern int32_t interactive_test_simplex();

extern int32_t minimum(simplex::Objective &objective, std::vector<simplex::Constraint> &constraints, uint64_t dead_var, simplex::Solution &solution);

extern int32_t reset();

extern void print_lp(simplex::Objective &EcFun, std::vector<simplex::Constraint> &constraints, uint64_t dead_var, simplex::Solution &Optimum);

extern void print_solution(simplex::Solution &solution, int result);
}
