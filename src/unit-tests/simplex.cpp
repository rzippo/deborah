#include <iostream>
#include "catch.hpp"

#include "../lib/simplex/constraint.h"
#include "../lib/simplex/simplex.h"
#include "../lib/tnode.h"

TEST_CASE ("Simplex: test_es4")
{
    uint64_t nvars = 5;
    int32_t dead_var = 0;
    simplex::Objective fobj(nvars);
    simplex::Constraint tmp(nvars);
    simplex::Solution opt(nvars);
    std::vector<simplex::Constraint> cset;
    bool result;
    double minvalue = 0.0;

    fobj.setCoefficient(0, 0.0); //dead variable!
    fobj.setCoefficient(1, 0.1);
    fobj.setCoefficient(2, 0.9);
    fobj.setCoefficient(3, 0.08);
    fobj.setCoefficient(4, 1.2);
    fobj.setConstant(0.0);

    tmp.Inequality.SetCoefficient(1, 9.9);
    tmp.Inequality.SetCoefficient(2, 20.0);
    tmp.Inequality.SetCoefficient(3, 2.0);
    tmp.Inequality.SetCoefficient(4, 22.0);
    tmp.Inequality.SetConstant(70.0);
    cset.push_back(tmp);

    tmp.Inequality.SetCoefficient(1, 1.2);
    tmp.Inequality.SetCoefficient(2, 1.1);
    tmp.Inequality.SetCoefficient(3, 0.2);
    tmp.Inequality.SetCoefficient(4, 26.0);
    tmp.Inequality.SetConstant(50.0);
    cset.push_back(tmp);

    tmp.Inequality.SetCoefficient(1, 75.3);
    tmp.Inequality.SetCoefficient(2, 0.0);
    tmp.Inequality.SetCoefficient(3, 4.0);
    tmp.Inequality.SetCoefficient(4, 0.0);
    tmp.Inequality.SetConstant(250.0);
    cset.push_back(tmp);

    tmp.Inequality.SetCoefficient(1, 0.0);
    tmp.Inequality.SetCoefficient(2, 0.0);
    tmp.Inequality.SetCoefficient(3, 1.0);
    tmp.Inequality.SetCoefficient(4, 0.0);
    tmp.Inequality.SetConstant(5.0);
    tmp.greater_than_zero = false;
    cset.push_back(tmp);

    std::cout << "*** TEST: ES4 PROBLEM" << std::endl;

    std::cout << "Problem in LP syntax:" << std::endl;
    simplex::print_lp(fobj, cset, dead_var, opt);
    std::cout << std::endl;

    result = simplex::minimum(fobj, cset, dead_var, opt);

    if (result)
    {
        std::cout << "Simplex has computed the optimum: FOBJ = " << opt.getConstant() << std::endl;
        std::cout << "Optimum variable values: " << std::endl;
        for (int i = 0; i < opt.getNumVariables(); i++)
            std::cout << "x_" << i << ": " << opt.getCoefficient(i) << std::endl;
    }
    else
    {
        std::cout << "Simplex did not find a solution" << std::endl;
    }
}
