#include "constraint.h"
#include <iostream>
#include <stdio.h>

namespace simplex
{

Constraint::Constraint(uint64_t num_variables) : Inequality(num_variables)
{
    greater_than_zero = true;
    nvars = num_variables;
    tnode_id = CPOLY_TNODEID_NULL;
}

Constraint::Constraint(uint64_t num_variables, uint64_t tnode) : Inequality(num_variables)
{
    greater_than_zero = true;
    nvars = num_variables;
    tnode_id = tnode;
}

Constraint::Constraint(const Constraint &cp) : Inequality(cp.nvars)
{
    Inequality = cp.Inequality;
    greater_than_zero = cp.greater_than_zero;
    nvars = Inequality.getNumVariables();
    tnode_id = cp.tnode_id;
}

Constraint::~Constraint()
{
}

Constraint &Constraint::operator=(const Constraint &cp)
{
    if (&cp != this)
    {
        Inequality = cp.Inequality;
        greater_than_zero = cp.greater_than_zero;
        tnode_id = cp.tnode_id;
    }
    return *this;
}

Constraint &Constraint::operator=(const deborah::Polynomial &p)
{
    if (&p != &Inequality)
    {
        Inequality = p;
        greater_than_zero = true;
    }
    return *this;
}

bool Constraint::isTautology()
{
    // only one type of check: const >= 0   or  const <= 0
    if (!Inequality.isConstant())
    {
        return false;
    }
    if (greater_than_zero)
    {
        return (Inequality.getConstant() >= 0);
    }
    else
    {
        return (Inequality.getConstant() <= 0);
    }
}

bool Constraint::isOxymoron()
{
    // only one type of check: const >= 0   or  const <= 0
    if (!Inequality.isConstant())
    {
        return false;
    }
    if (greater_than_zero)
    {
        return (Inequality.getConstant() < 0);
    }
    else
    {
        return (Inequality.getConstant() > 0);
    }
}

// This method flips the direction of the Inequality by negating the coefficients
void Constraint::InvertDirection()
{
    Inequality.Negate();
    greater_than_zero = !greater_than_zero;
}

/*	Checks whether the constraint is respected by the given point
*/
bool Constraint::isSatisfied(deborah::Polynomial &p)
{
    if (nvars != p.getNumVariables())
    {
        std::cout << "Constraint::isSatisfied(): num_variables mismatch" << std::endl;
        return false;
    }
    double val = Inequality.getConstant();
    for (uint64_t i = 0; i < nvars; i++)
    {
        val += Inequality.getCoefficient(i) * p.getCoefficient(i);
    }
    if (greater_than_zero)
    {
        return (val >= 0.0);
    }
    else
    {
        return (val <= 0.0);
    }
}

void Constraint::Print()
{
    Inequality.Print();
    if (greater_than_zero)
    {
        std::cout << " >= 0";
    }
    else
    {
        std::cout << " <= 0";
    }

    std::cout << "  (id=" << tnode_id << ")";
}

}
