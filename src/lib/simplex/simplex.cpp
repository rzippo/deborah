#include "simplex.h"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>

/***************************************************************
*           LINEAR PROGRAMMING: THE SIMPLEX METHOD             *
* ------------------------------------------------------------ *
* SAMPLE RUN:                                                  *
* Maximize z = x1 + x2 + 3x3 -0.5x4 with conditions:           * 
*          x1  + 2x3 <= 740                                    *
*          2x2 - 7x4 <= 0                                      *
*          x2  - x3 + 2x4 >= 0.5                               *
*          x1 + x2 + x3 +x4 = 9                                *
*          and all x's >=0.                                    *
*                                                              *
* Number of variables in E.F.: 4                               *
* Number of <= inequalities..: 2                               *
* Number of >= inequalities..: 1                               *
* Number of = equalities.....: 1                               *
* Input Economic Function:                                     *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: 3                                           *
* Coefficient # 4: -0.5                                        *
* Constant term..: 0                                           *
* Input constraint # 1:                                        *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 0                                           *
* Coefficient # 3: 2                                           *
* Coefficient # 4: 0                                           *
* Constant term..: 740                                         *
* Input constraint # 2:                                        *
* Coefficient # 1: 0                                           *
* Coefficient # 2: 2                                           *
* Coefficient # 3: 0                                           *
* Coefficient # 4: -7                                          *
* Constant term..: 0                                           *
* Input constraint # 3:                                        *
* Coefficient # 1: 0                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: -1                                          *
* Coefficient # 4: 2                                           *
* Constant term..: 0.5                                         *
* Input constraint # 4:                                        *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: 1                                           *
* Coefficient # 4: 1                                           *
* Constant term..: 9                                           *
*                                                              *
* Input Table:                                                 *
*    0.00    1.00    1.00    3.00   -0.50                      *
*  740.00   -1.00    0.00   -2.00    0.00                      *
*    0.00    0.00   -2.00    0.00    7.00                      *
*    0.50    0.00   -1.00    1.00   -2.00                      *
*    9.00   -1.00   -1.00   -1.00   -1.00                      *
*                                                              *
* Maximum of E.F. =    17.02500                                *
*  X1 =    0.000000                                            *
*  X2 =    3.325000                                            *
*  X3 =    4.725000                                            *
*  X4 =    0.950000                                            *
*                                                              *
* ------------------------------------------------------------ *
* Reference: "Numerical Recipes By W.H. Press, B. P. Flannery, *
*             S.A. Teukolsky and W.T. Vetterling, Cambridge    *
*             University Press, 1986" [BIBLI 08].              *
*                                                              *
*                       C++ Release 1.0 By J-P Moreau, Paris   *
***************************************************************/

namespace simplex
{

#define  REAL  double

typedef REAL MAT[MMAX][NMAX];

static MAT A;
static int32_t IPOSV[MMAX], IZROV[NMAX];
static int32_t i, j, ICASE, N, M, M1, M2, M3;
static REAL R;

void simp1(MAT, int32_t, int32_t *, int32_t, int32_t, int32_t *, REAL *);

void simp2(MAT, int32_t, int32_t, int32_t *, int32_t, int32_t *, int32_t, REAL *);

void simp3(MAT, int32_t, int32_t, int32_t, int32_t);

void simplx(
    MAT a, int32_t m, int32_t n, int32_t m1, int32_t m2, int32_t m3, int32_t *icase, int32_t *izrov,
    int32_t *iposv
)
{
    /*----------------------------------------------------------------------------------------
         USES simp1,simp2,simp3
         Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3,
         and output parameters a, icase, izrov, and iposv are described above (see reference).
         Parameters: MMAX is the maximum number of constraints expected; NMAX is the maximum number
         of variables expected; EPS is the absolute precision, which should be adjusted to the
         scale of your variables.
         -----------------------------------------------------------------------------------------*/
    int32_t i, ip, ir, is, k, kh, kp, m12, nl1, nl2, l1[NMAX], l2[MMAX], l3[MMAX];
    REAL bmax, q1, EPS = 1e-6;
    if (m != m1 + m2 + m3)
    {
        std::cout << " Bad input constraint counts in simplx." << std::endl;
        return;
    }
    nl1 = n;
    for (k = 1; k <= n; k++)
    {
        l1[k] = k;     //Initialize index list of columns admissible for exchange.
        izrov[k] = k;  //Initially make all variables right-hand.
    }
    nl2 = m;
    for (i = 1; i <= m; i++)
    {
        if (a[i + 1][1] < 0.0)
        {
            std::cout << " Bad input tableau in simplx, Constants bi must be nonnegative." << std::endl;
            return;
        }
        l2[i] = i;
        iposv[i] = n + i;
/*------------------------------------------------------------------------------------------------
		  Initial left-hand variables. m1 type constraints are represented by having their slackv ariable 
		  initially left-hand, with no artificial variable. m2 type constraints have their slack 
		  variable initially left-hand, with a minus sign, and their artificial variable handled implicitly 
		  during their first exchange. m3 type constraints have their artificial variable initially 
		  left-hand.     
		  ------------------------------------------------------------------------------------------------*/
    }
    for (i = 1; i <= m2; i++)
    {
        l3[i] = 1;
    }
    ir = 0;
    if (m2 + m3 == 0)
    {
        goto e30;
    } //The origin is a feasible starting solution. Go to phase two.
    ir = 1;
    for (k = 1; k <= n + 1; k++)
    { //Compute the auxiliary objective function.
        q1 = 0.0;
        for (i = m1 + 1; i <= m; i++)
        {
            q1 += a[i + 1][k];
        }
        a[m + 2][k] = -q1;
    }
    e10:
    simp1(a, m + 1, l1, nl1, 0, &kp, &bmax);    //Find max. coeff. of auxiliary objective fn
    if (bmax <= EPS && a[m + 2][1] < -EPS)
    {
        *icase = -1;    //Auxiliary objective function is still negative and can’t be improved,
        return;       //hence no feasible solution exists.
    }
    else if (bmax <= EPS && a[m + 2][1] <= EPS)
    {
//Auxiliary objective function is zero and can’t be improved; we have a feasible starting vector. 
//Clean out the artificial variables corresponding to any remaining equality constraints by 
//goto 1’s and then move on to phase two by goto 30.
        m12 = m1 + m2 + 1;
        if (m12 <= m)
        {
            for (ip = m12; ip <= m; ip++)
            {
                if (iposv[ip] == ip + n)
                {   //Found an artificial variable for an equalityconstraint.
                    simp1(a, ip, l1, nl1, 1, &kp, &bmax);
                    if (bmax > EPS)
                    {
                        goto e1;
                    } //Exchange with column corresponding to maximum
                }
            }
        }                         //pivot element in row.
        ir = 0;
        m12 = m12 - 1;
        if (m1 + 1 > m12)
        {
            goto e30;
        }
        for (i = m1 + 1; i <= m1 + m2; i++)
        {     //Change sign of row for any m2 constraints
            if (l3[i - m1] == 1)
            {             //still present from the initial basis.
                for (k = 1; k <= n + 1; k++)
                {
                    a[i + 1][k] *= -1.0;
                }
            }
        }
        goto e30;                       //Go to phase two.
    }

    simp2(a, m, n, l2, nl2, &ip, kp, &q1); //Locate a pivot element (phase one).

    if (ip == 0)
    {                         //Maximum of auxiliary objective function is
        *icase = -1;                          //unbounded, so no feasible solution exists.
        return;
    }
    e1:
    simp3(a, m + 1, n, ip, kp);
//Exchange a left- and a right-hand variable (phase one), then update lists. 
    if (iposv[ip] >= n + m1 + m2 + 1)
    { //Exchanged out an artificial variable for an
        //equality constraint. Make sure it stays
        //out by removing it from the l1 list.
        for (k = 1; k <= nl1; k++)
        {
            if (l1[k] == kp)
            {
                goto e2;
            }
        }
        e2:
        nl1 = nl1 - 1;
        for (is = k; is <= nl1; is++)
        {
            l1[is] = l1[is + 1];
        }
    }
    else
    {
        if (iposv[ip] < n + m1 + 1)
        {
            goto e20;
        }
        kh = iposv[ip] - m1 - n;
        if (l3[kh] == 0)
        {
            goto e20;
        }  //Exchanged out an m2 type constraint.
        l3[kh] = 0;                  //If it’s the first time, correct the pivot column
        //or the minus sign and the implicit
        //artificial variable.
    }
    a[m + 2][kp + 1] += 1.0;
    for (i = 1; i <= m + 2; i++)
    {
        a[i][kp + 1] *= -1.0;
    }
    e20:
    is = izrov[kp];             //Update lists of left- and right-hand variables.
    izrov[kp] = iposv[ip];
    iposv[ip] = is;
    if (ir != 0)
    {
        goto e10;
    }       //if still in phase one, go back to 10.
//End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it. 
    e30:
    simp1(a, 0, l1, nl1, 0, &kp, &bmax); //Test the z-row for doneness.
    if (bmax <= EPS)
    {          //Done. Solution found. Return with the good news.
        *icase = 0;
        return;
    }
    simp2(a, m, n, l2, nl2, &ip, kp, &q1);   //Locate a pivot element (phase two).
    if (ip == 0)
    {             //Objective function is unbounded. Report and return.
        *icase = 1;
        return;
    }
    simp3(a, m, n, ip, kp);       //Exchange a left- and a right-hand variable (phase two),
    goto e20;                 //update lists of left- and right-hand variables and
}                           //return for another iteration.

// The preceding routine makes use of the following utility subroutines: 

void simp1(MAT a, int32_t mm, int32_t *ll, int32_t nll, int32_t iabf, int32_t *kp, REAL *bmax)
{
//Determines the maximum of those elements whose index is contained in the supplied list 
//ll, either with or without taking the absolute value, as flagged by iabf. 
    int32_t k;
    REAL test;
    *kp = ll[1];
    *bmax = a[mm + 1][*kp + 1];
    if (nll < 2)
    {
        return;
    }
    for (k = 2; k <= nll; k++)
    {
        if (iabf == 0)
        {
            test = a[mm + 1][ll[k] + 1] - (*bmax);
        }
        else
        {
            test = std::fabs(a[mm + 1][ll[k] + 1]) - std::fabs(*bmax);
        }
        if (test > 0.0)
        {
            *bmax = a[mm + 1][ll[k] + 1];
            *kp = ll[k];
        }
    }
    return;
}

void simp2(MAT a, int32_t m, int32_t n, int32_t *l2, int32_t nl2, int32_t *ip, int32_t kp, REAL *q1)
{
    REAL EPS = 1e-6;
//Locate a pivot element, taking degeneracy into account. 
    int32_t i, ii, k;
    REAL q, q0, qp;
    *ip = 0;
    if (nl2 < 1)
    {
        return;
    }
    for (i = 1; i <= nl2; i++)
    {
        if (a[i + 1][kp + 1] < -EPS)
        {
            goto e2;
        }
    }
    return;  //No possible pivots. Return with message.
    e2:
    *q1 = -a[l2[i] + 1][1] / a[l2[i] + 1][kp + 1];
    *ip = l2[i];
    if (i + 1 > nl2)
    {
        return;
    }
    for (i = i + 1; i <= nl2; i++)
    {
        ii = l2[i];
        if (a[ii + 1][kp + 1] < -EPS)
        {
            q = -a[ii + 1][1] / a[ii + 1][kp + 1];
            if (q < *q1)
            {
                *ip = ii;
                *q1 = q;
            }
            else if (q == *q1)
            {  //We have a degeneracy.
                for (k = 1; k <= n; k++)
                {
                    qp = -a[*ip + 1][k + 1] / a[*ip + 1][kp + 1];
                    q0 = -a[ii + 1][k + 1] / a[ii + 1][kp + 1];
                    if (q0 != qp)
                    {
                        goto e6;
                    }
                }
                e6:
                if (q0 < qp)
                {
                    *ip = ii;
                }
            }
        }
    }
    return;
}

void simp3(MAT a, int32_t i1, int32_t k1, int32_t ip, int32_t kp)
{
//Matrix operations to exchange a left-hand and right-hand variable (see text). 
    int32_t ii, kk;
    REAL piv;
    piv = 1.0 / a[ip + 1][kp + 1];
    if (i1 >= 0)
    {
        for (ii = 1; ii <= i1 + 1; ii++)
        {
            if (ii - 1 != ip)
            {
                a[ii][kp + 1] *= piv;
                for (kk = 1; kk <= k1 + 1; kk++)
                {
                    if (kk - 1 != kp)
                    {
                        a[ii][kk] -= a[ip + 1][kk] * a[ii][kp + 1];
                    }
                }
            }
        }
    }
    for (kk = 1; kk <= k1 + 1; kk++)
    {
        if (kk - 1 != kp)
        {
            a[ip + 1][kk] = -a[ip + 1][kk] * piv;
        }
    }
    a[ip + 1][kp + 1] = piv;
    return;
}

/***************************************************************************************************************/

[[maybe_unused]] extern void dump_table()
{
    std::cout << "Input Table:" << std::endl;
    std::cout << "------------------------------" << std::endl;
    for (i = 1; i <= M + 1; i++)
    {
        for (j = 1; j <= N + 1; j++)
        {
            std::cout << std::setw(8) << std::setprecision(2) << A[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << "------------------------------" << std::endl;
}

/*	This is a test function which asks for problem data interactively (keyboard input)
        It computes the MAX of the economic function
*/
[[maybe_unused]] int32_t interactive_test_simplex()
{
    std::cout << std::endl;
    std::cout << " Number of variables in E.F.: ";
    std::cin >> N;
    std::cout << " Number of <= inequalities..: ";
    std::cin >> M1;
    std::cout << " Number of >= inequalities..: ";
    std::cin >> M2;
    std::cout << " Number of = equalities.....: ";
    std::cin >> M3;

    M = M1 + M2 + M3;  // Total number of constraints

    for (i = 1; i <= M + 2; i++)
    {
        for (j = 1; j <= N + 1; j++)
        {
            A[i][j] = 0.0;
        }
    }

    std::cout << " Input Economic Function to MAXIMIZE:" << std::endl;
    for (i = 2; i <= N + 1; i++)
    {
        std::cout << " Coefficient #" << i - 1 << ": ";
        std::cin >> A[1][i];
    }

    std::cout << " Constant term : ";
    std::cin >> A[1][1];
// input constraints    
    for (i = 1; i <= M; i++)
    {
        std::cout << " Input constraint #" << i << ": " << std::endl;
        for (j = 2; j <= N + 1; j++)
        {
            std::cout << " Coefficient #" << j - 1 << ": ";
            std::cin >> R;
            A[i + 1][j] = -R;
        }

        std::cout << " Constant term : ";
        std::cin >> A[i + 1][1];
    }

    std::cout << "\n Input Table:" << std::endl;
    for (i = 1; i <= M + 1; i++) \

    {
        for (j = 1; j <= N + 1; j++)
        {
            std::cout << std::setprecision(2) << std::setw(8) << A[i][j] << " ";
        }

        std::cout << std::endl;
    }

    simplx(A, M, N, M1, M2, M3, &ICASE, IZROV, IPOSV);

    if (ICASE == 0)
    {  //result ok.
        std::cout << "\n Maximum of E.F. = " << A[1][1] << std::endl;

        for (i = 1; i <= N; i++)
        {
            bool print0 = true;
            for (j = 1; j <= M; j++)
            {
                if (IPOSV[j] == i)
                {
                    std::cout << "  X" << i << " = " << A[j + 1][1] << std::endl;
                    print0 = false;
                    break;
                }
            }

            if (print0)
            {
                std::cout << "  X" << i << " = " << 0.0 << std::endl;
            }
        }
    }
    else
    {
        std::cout << " No solution (error code = " << ICASE << ")." << std::endl;
    }

    std::cout << std::endl;
    return 0;
}

//#define DEBUG_SIMPLEX

/// This method runs the simplex algorithm to minimize the passed objective function
/// \param objective Economic Function as returned by the TNode class (the first variable will be skipped)
/// \param constraints set of constraints in the form ax+b<=0 or ax+b<=0 (x>=0 are assumed implicitly)
/// \param dead_var the variable at this index must be excluded from the simplex
/// \param solution optimal solution
/// \return Error code, see SIMPLEX_RESULT macros
int32_t minimum(simplex::Objective &objective, std::vector<simplex::Constraint> &constraints, const uint64_t dead_var, simplex::Solution &solution)
{
#ifdef DEBUG_SIMPLEX
    print_lp(objective, constraints, dead_var, solution);
#endif

    int32_t result = 0;
    int32_t nvars = objective.getNumVariables();
    uint64_t idx;
    int32_t num_constraints;
    std::vector<simplex::Constraint> gte;
    std::vector<simplex::Constraint> lte;
    std::vector<double> solution_values(nvars);

    if (nvars > NMAX)
    {
        std::cout << "Simplex error: max number of variables exceeded (" << nvars << ", max " << NMAX << ")"
                  << std::endl;
        return SIMPLEX_RESULT_UNFEASIBLE;
    }

    //solution.Zero();
    num_constraints = constraints.size();
    if (num_constraints > MMAX)
    {
        std::cout << "Simplex error: max number of constraints exceeded (N=" << num_constraints << ")" << std::endl;
        return SIMPLEX_RESULT_UNFEASIBLE;
    }

    // split the constraint set into two distinct sets (for <= and >= types)
    for (i = 0; i < constraints.size(); i++)
    {
        Constraint c1(constraints[i]);
        if (-c1.Inequality.getConstant() < 0.0)
        {
            c1.InvertDirection();
        } // negate B because the library expects ax<=b or ax>=b
        if (c1.greater_than_zero)
        {
            gte.push_back(c1);
        }
        else
        {
            lte.push_back(c1);
        }
    }

    N = nvars - 1; // discard the dead variable
    M1 = lte.size(); // number of <= constraints
    M2 = gte.size(); // number of >= constraints
    M3 = 0; // number of == constraints
    M = M1 + M2 + M3;  // total number of constraints

    // reset problem data
    for (i = 1; i <= M + 2; i++)
    {
        for (j = 1; j <= N + 1; j++)
        {
            A[i][j] = 0.0;
        }
    }

    // store economic function data (negate all coefficients, since min(cx) = -max(-cx)
    A[1][1] = -objective.getConstant(); // constant term
    idx = 2;
    for (i = 0; i < N + 1; i++)
    {
        if (i == dead_var)
        {
            continue;
        }
        A[1][idx] = -objective.getCoefficient(i);
        idx++;
    }

    // store constraints data
    for (i = 1; i <= M; i++)
    {
        Constraint &cp = (i <= M1) ? lte[i - 1] : gte[i - M1 - 1];
        //ith constraint
        A[i + 1][1] = -cp.Inequality.getConstant(); // constant term: negate it because the library expects ax<=b or ax>=b
        idx = 2;
        for (j = 0; j < N + 1; j++)
        {
            if (j == dead_var)
            {
                continue;
            }
            R = cp.Inequality.getCoefficient(j);
            A[i + 1][idx] = -R;
            idx++;
        }
    }

    //dump_table();

    // run the simplex!
    simplx(A, M, N, M1, M2, M3, &ICASE, IZROV, IPOSV);

    if (ICASE == 0)
    {
        // simplex was feasible
        solution.setConstant(-A[1][1]);

        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= M; j++)
            {
                if (IPOSV[j] == i)
                {
                    solution_values[i - 1] = A[j + 1][1];
                    goto e3;
                }
            }
            solution_values[i - 1] = 0.0;
            e3:;
        }

        j = 0;
        for (i = 0; i < solution.getNumVariables(); i++)
        {
            if (i == dead_var)
            {
                continue;
            }
            solution.setCoefficient(i, solution_values[j]);
            j++;
        }
        //std::cout << "SOLUTION = "); solution.Print(); std::cout << std::endl;
        result = SIMPLEX_RESULT_OK;
    }
    else
    {
        // simplex was unfeasible
        for (i = 0; i < N; i++)
        {
            solution_values[i] = 0.0;
        }
        if (ICASE == 1)
        {
            result = SIMPLEX_RESULT_UNBOUNDED;
        }
        else if (ICASE == -1)
        {
            result = SIMPLEX_RESULT_UNFEASIBLE;
        }
        else
        {
            result = SIMPLEX_RESULT_UNK_ERROR;
        }
    }

#ifdef DEBUG_SIMPLEX
    print_solution(solution, result);
#endif

    return result;
}

std::string format_value(double value)
{
    std::stringstream ss;
    if(value >= 0)
    {
        ss << "+" << value;
    }
    else
    {
        ss << value;
    }
    return ss.str();
}

void print_lp(simplex::Objective &EcFun, std::vector<simplex::Constraint> &constraints, uint64_t dead_var, simplex::Solution &Optimum)
{
    std::cout << "Minimize" << std::endl;
    for(int var_idx = 0; var_idx < EcFun.getNumVariables(); var_idx++)
    {
        if(var_idx == dead_var)
                continue;
        std::cout << format_value(EcFun.getCoefficient(var_idx)) << " x_" << var_idx << " ";
    }
    std::cout << format_value(EcFun.getConstant()) << std::endl;

    std::cout << std::endl;
    std::cout << "Subject To" << std::endl;
    for(auto constraint : constraints)
    {
        for (int var_idx = 0; var_idx < constraint.Inequality.getNumVariables(); ++var_idx)
        {
            if(var_idx == dead_var)
                continue;
            std::cout << format_value(constraint.Inequality.getCoefficient(var_idx)) << " x_" << var_idx << " ";
        }
        std::cout << " >= " << format_value(-constraint.Inequality.getConstant()) << std::endl;
    }
    std::cout << "End" << std::endl;
}

void print_solution(simplex::Solution &solution, int dead_var, int result)
{
    std::cout << "Solution found by integrated algorithm:" << std::endl;
    switch (result)
    {
        case SIMPLEX_RESULT_OK:
        {
            std::cout << "Problem is FEASIBLE" << std::endl;
            std::cout << "Objective value: " << solution.getConstant() << std::endl;
            for (int var_idx = 0; var_idx < solution.getNumVariables(); ++var_idx)
            {
                if(var_idx == dead_var)
                    continue;
                std::cout << "x_" << var_idx << " = " << solution.getCoefficient(var_idx) << std::endl;
            }
            break;
        }

        case SIMPLEX_RESULT_UNFEASIBLE:
        {
            std::cout << "Problem is UNFEASIBLE" << std::endl;
            break;
        }

        case SIMPLEX_RESULT_UNBOUNDED:
        {
            std::cout << "Problem is UNFEASIBLE" << std::endl;
            break;
        }

        case SIMPLEX_RESULT_UNK_ERROR:
        {
            std::cout << "UNKNOWN ERROR during simplex algorithm" << std::endl;
            break;
        }

        default:
            std::cout << "!! Unrecognized return code !!" << std::endl;
            break;
    }
}

} //end of namespace simplex
