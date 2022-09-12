#include "linearsegment.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>

namespace deborah
{

LinearSegment::LinearSegment()
{
    x = 0.0;
    y = 0.0;
    rate = LINEARSEGMENT_POSITIVE_INFINITE;
    leftopen = false;
    id = -1;
    gap = 0.0;
}

LinearSegment::LinearSegment(double x0, double r)
{
    x = x0;
    y = 0.0;
    rate = r;
    leftopen = false;
    id = -1;
    gap = 0.0;
}

LinearSegment::LinearSegment(double x0, double y0, double r)
{
    x = x0;
    y = y0;
    rate = r;
    leftopen = false;
    id = -1;
    gap = 0.0;
}

LinearSegment &LinearSegment::operator=(const LinearSegment &s2)
{
    if (this != &s2)
    {
        x = s2.x;
        y = s2.y;
        rate = s2.rate;
        leftopen = s2.leftopen;
        id = s2.id;
        gap = s2.gap;
    }
    return *this;
}

LinearSegment::LinearSegment(const LinearSegment &ls)
{
    x = ls.x;
    y = ls.y;
    rate = ls.rate;
    leftopen = ls.leftopen;
    id = ls.id;
    gap = ls.gap;
}

LinearSegment::~LinearSegment()
{
}

void LinearSegment::setLatencyRate(double l, double r)
{
    x = l;
    y = 0.0;
    rate = r;
}

void LinearSegment::setAffine(double b, double r)
{
    x = 0.0;
    y = b;
    rate = r;
}

double LinearSegment::value(double x1)
{
    if (std::fabs(x1 - x) < LINEARSEGMENT_EPSILON)
    {
        return y;
    }
    return (x1 - x) * rate + y;
}

double LinearSegment::get_X_intersection(LinearSegment &other)
{
    if (other.rate == LINEARSEGMENT_POSITIVE_INFINITE)
    {
        return LINEARSEGMENT_POSITIVE_INFINITE;
    }
    if (std::fabs(rate - other.rate) < LINEARSEGMENT_EPSILON)
    {
        return LINEARSEGMENT_POSITIVE_INFINITE;
    }
    double y1 = y - x * rate;
    double y2 = other.y - other.x * other.rate;
    return (y2 - y1) / (rate - other.rate);
}

bool LinearSegment::convolve(double x1, double rate1)
{
    if (y != 0.0)
    {
        return false;
    }

    x += x1;
    if ((rate1 < rate) || (rate == LINEARSEGMENT_POSITIVE_INFINITE))
    {
        rate = rate1;
    }

    return true;
}

bool LinearSegment::convolve(LinearSegment &other)
{
    //supports only latency-rate segments
    if ((y != 0.0) || (other.y != 0.0))
    {
        return false;
    }
    return convolve(other.x, other.rate);
}

void LinearSegment::add(LinearSegment &s2, double xpos, bool lopen)
{
    x = xpos;
    y = value(xpos) + s2.value(xpos);
    rate = rate + s2.rate;
    leftopen = lopen;
}

void LinearSegment::sub(LinearSegment &s2, double xpos, bool lopen)
{
    x = xpos;
    y = value(xpos) - s2.value(xpos);
    rate = rate - s2.rate;
    leftopen = lopen;
}

void LinearSegment::Zero()
{
    x = y = rate = 0.0;
    leftopen = false;
    id = -1;
    gap = 0.0;
}

void LinearSegment::Print()
{
    std::cout
        << "(" << std::setprecision(6) << x
        << ", " << std::setprecision(6) << y
        << ")  R = ";

    if (rate == LINEARSEGMENT_POSITIVE_INFINITE)
    {
        std::cout << "+INF";
    }
    else if (rate == LINEARSEGMENT_NEGATIVE_INFINITE)
    {
        std::cout << "-INF";
    }
    else
    {
        std::cout << std::setprecision(6) << rate;
    }

    //std::cout << "  leftopen=%s  flow_burst_id=" << (leftopen?"true":"false"),id,gap << "  gap=%.3lf";
    if (id >= 0)
    {
        std::cout
            << "   bursting_flow = %02i" << std::setw(2) << std::setfill('0') << id
            << "   burst_size = %.3lf" << std::setprecision(3) << gap;
    }
}

}
