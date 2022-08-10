#include "node.h"

namespace deborah
{

Node::Node()
{
    latency = 0.0;
    rate = 0.0;
    CDF_i = nullptr;
    CAF_i = nullptr;
    max_output_rates = nullptr;
}

Node::Node(uint64_t src, uint64_t exit, double theta, double R)
{
    latency = theta;
    rate = R;
    max_output_rates = nullptr;
}

Node::~Node()
{
    if (CDF_i != nullptr)
    {
        delete[] CDF_i;
    }
    CDF_i = nullptr;
    if (CAF_i != nullptr)
    {
        delete[] CAF_i;
    }
    CAF_i = nullptr;
    if (max_output_rates != nullptr)
    {
        delete[] max_output_rates;
    }
    max_output_rates = nullptr;
}

void Node::clear_curves(uint64_t n_flows)
{
    CDF.Zero();
    CAF.Zero();
    beta.create_latency_rate(latency, rate);
    if (CDF_i != nullptr)
    {
        for (int32_t i = 0; i < n_flows; i++)
        {
            CDF_i[i].Zero();
        }
    }
    else
    {
        CDF_i = new Curve[n_flows];
    }
    if (CAF_i != nullptr)
    {
        for (int32_t i = 0; i < n_flows; i++)
        {
            CAF_i[i].Zero();
        }
    }
    else
    {
        CAF_i = new Curve[n_flows];
    }
}

bool Node::set_output_rates(uint64_t n_flows, std::vector<double> out_rates, bool rate_red)
{
    if (max_output_rates == nullptr)
    {
        max_output_rates = new double[n_flows];
    }
    if (max_output_rates == nullptr)
    {
        return false;
    }

    if (out_rates.empty())
    {
        for (int32_t i = 0; i < n_flows; i++)
        {
            max_output_rates[i] = -1.0;
        }
    }

    else
    {
        for (int32_t i = 0; i < n_flows; i++)
        {
            max_output_rates[i] = out_rates[i];
        }
    }
    rate_reduction = rate_red;
    return true;
}

/*	Return the max output rate of the specified flow at this node
	A negative value means that the flow does not interfere with the tagged flow at this node
*/
double Node::get_output_rate(uint64_t flow)
{
    if (max_output_rates == nullptr)
    {
        return 0.0;
    }
    return max_output_rates[flow];
}

bool Node::set_output_rate(uint64_t flow, double rate)
{
    if (max_output_rates == nullptr)
    {
        return false;
    }

    max_output_rates[flow] = rate;
    return true;
}

bool Node::is_rate_reducing()
{
    return rate_reduction;
}

} // namespace deborah
