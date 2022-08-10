/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   								   *
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
#include "flow.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

namespace deborah
{

Flow::Flow()
{
    src = 0;
    exit = 0;
    nesting_level = 0;
    burst = 0.0;
    rate = 0.0;

    uid = FLOW_INVALID_ID;
}

Flow::~Flow()
{
}

/* Check if this flow is nested into the passed one
 *
 */
bool Flow::is_nested(Flow &f)
{
    if ((src >= f.src) && (exit <= f.exit))
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Check if this flow intersects the passed one
bool Flow::intersects(Flow &f)
{
    // check if (i,j) and (h,k) are such that: i < h <= j < k
    return ((src < f.src) && (f.src <= exit) && (exit < f.exit));
}

bool Flow::contains(uint64_t node_id)
{
    if ((node_id >= src) && (node_id <= exit))
    {
        return true;
    }
    return false;
}

std::string Flow::toString()
{
    std::ostringstream ss;
    ss << "(" << src + 1 << "," << exit + 1 << ")";
    return ss.str();
}

bool Flow::create_caf(double x1, double x2, bool ending_burst, uint64_t flow_id)
{
    bool res;

    if ((x1 == 0.0) && (x2 == 0.0))
    {
        // it's the tagged flow: single burst at zero
        res = caf.create_burst_line(0.0, 0.0, 0.0, burst, flow_id);
        return res;
    }

    if (ending_burst)
    {
        res = caf.create_line_burst(x1, x2, rate, burst, flow_id);
    }
    else
    {
        res = caf.create_burst_line(x1, x2, rate, burst, flow_id);
    }
    return res;
}

void Flow::Print()
{
    std::cout << toString();
}

}
