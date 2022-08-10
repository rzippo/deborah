#include "catch.hpp"

#include <iostream>

#include "../lib/tandem.h"

/*	*** TEST code for LowerBound-related classes and operations ***
*/
TEST_CASE( "LowerBound-related classes and operations", "[tandem]" )
{
    Curve beta, caf1, caf2, cdf;
    LinearSegment ls;
    bool res;

    auto flows = std::vector<Flow>(2);
    flows[0].rate = 40.0;
    flows[0].burst = 4.0;
    flows[1].rate = 60.0;
    flows[1].burst = 6.0;

    SECTION( "Two final bursts")
    {
        flows[0].create_caf(10.0, 11.0, true, 0);
        caf1 = flows[0].caf;
        std::cout << "CAF1 ";
        caf1.Print();
        caf2.create_token_bucket(0.0, 10.0); // sigma, rho
        flows[1].create_caf(10.0, 11.0, true, 1);
        caf2 = flows[1].caf;
        std::cout << "CAF2 ";
        caf2.Print();
        caf2.add(caf1);
        std::cout << "CAF total ";
        caf2.Print();
        beta.create_latency_rate(10.0, 10.0); // theta, R
        std::cout << "beta ";
        beta.Print();
        cdf = caf2;
        res = cdf.ConvolveWithLatencyRate(beta);

        REQUIRE (res);

        std::cout << "CDF ";
        cdf.Print();
    }

    SECTION( "TWO INITIAL BURSTS")
    {
        flows[0].create_caf(10.0, 11.0, false, 0);
        caf1 = flows[0].caf;
        std::cout << "CAF1 ";
        caf1.Print();
        caf2.create_token_bucket(0.0, 10.0); // sigma, rho
        flows[1].create_caf(10.0, 11.0, false, 1);
        caf2 = flows[1].caf;
        std::cout << "CAF2 ";
        caf2.Print();
        caf2.add(caf1);
        std::cout << "CAF total ";
        caf2.Print();
        beta.create_latency_rate(10.0, 10.0); // theta, R
        std::cout << "beta ";
        beta.Print();
        cdf = caf2;
        res = cdf.ConvolveWithLatencyRate(beta);

        REQUIRE (res);

        std::cout << "CDF ";
        cdf.Print();
    }

    SECTION( "TWO AT BEGINNING, TWO AT END")
    {
        flows[0].create_caf(10.0, 11.0, true, 0);
        caf1 = flows[0].caf;
        std::cout << "CAF1 ";
        caf1.Print();
        caf2.create_token_bucket(0.0, 10.0); // sigma, rho
        flows[1].create_caf(10.0, 11.0, false, 1);
        caf2 = flows[1].caf;
        std::cout << "CAF2 ";
        caf2.Print();
        caf2.add(caf1);
        caf1 = caf2;
        caf2.add(caf1);
        std::cout << "CAF total ";
        caf2.Print();
        beta.create_latency_rate(10.0, 10.0); // theta, R
        std::cout << "beta ";
        beta.Print();
        cdf = caf2;
        res = cdf.ConvolveWithLatencyRate(beta);

        REQUIRE (res);

        std::cout << "CDF ";
        cdf.Print();
    }

    SECTION( "ONE AT BEGINNING, TAGGED FLOW")
    {
        flows[0].create_caf(0.0, 0.0, true, 0);
        caf1 = flows[0].caf;
        std::cout << "CAF1 ";
        caf1.Print();
        caf2.create_token_bucket(0.0, 10.0); // sigma, rho
        flows[1].create_caf(0.0, 0.0, true, 1);
        caf2 = flows[1].caf;
        std::cout << "CAF2 ";
        caf2.Print();
        caf2.add(caf1);
        std::cout << "CAF total ";
        caf2.Print();
        beta.create_latency_rate(10.0, 10.0); // theta, R
        std::cout << "beta ";
        beta.Print();
        cdf = caf2;
        res = cdf.ConvolveWithLatencyRate(beta);

        REQUIRE (res);

        std::cout << "CDF ";
        cdf.Print();
    }

    SECTION( "ONE CONTINUOS, ONE GAP")
    {
        flows[0].burst = 0.0;
        flows[0].create_caf(10.0, 12.0, true, 0);
        caf1 = flows[0].caf;
        std::cout << "CAF1 ";
        caf1.Print();
        caf2.create_token_bucket(0.0, 10.0); // sigma, rho
        flows[1].create_caf(10.0, 11.0, true, 1);
        caf2 = flows[1].caf;
        std::cout << "CAF2 ";
        caf2.Print();
        caf2.add(caf1);
        std::cout << "CAF total ";
        caf2.Print();
        beta.create_latency_rate(10.0, 10.0); // theta, R
        std::cout << "beta ";
        beta.Print();
        cdf = caf2;
        res = cdf.ConvolveWithLatencyRate(beta);

        REQUIRE (res);

        std::cout << "CDF ";
        cdf.Print();
    }

    SECTION( "TWO CONTINUOS")
    {
        flows[0].burst = 0.0;
        flows[1].burst = 0.0;
        flows[0].create_caf(10.0, 12.0, true, 0);
        caf1 = flows[0].caf;
        std::cout << "CAF1 ";
        caf1.Print();
        caf2.create_token_bucket(0.0, 10.0); // sigma, rho
        flows[1].create_caf(11.0, 13.0, true, 1);
        caf2 = flows[1].caf;
        std::cout << "CAF2 ";
        caf2.Print();
        caf2.add(caf1);
        std::cout << "CAF total ";
        caf2.Print();
        beta.create_latency_rate(10.0, 10.0); // theta, R
        std::cout << "beta ";
        beta.Print();
        cdf = caf2;
        res = cdf.ConvolveWithLatencyRate(beta);

        REQUIRE (res);

        std::cout << "CDF ";
        cdf.Print();
    }
}
