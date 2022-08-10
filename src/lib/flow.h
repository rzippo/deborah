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
#ifndef DEBORAHFLOW_H
#define DEBORAHFLOW_H

#include <string>
#include "curve.h"

namespace deborah
{

#define FLOW_INVALID_ID    0xFFFFFFFF

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class Flow
{
 public:
  Flow();

  ~Flow();

  bool is_nested(Flow &f);

  bool intersects(Flow &f);

  bool contains(uint64_t node_id);

  uint64_t src; // source node of the flow (zero-based)
  uint64_t exit; // exit node of the flow (zero-based)
  uint64_t nesting_level;
  double burst;
  double rate;

  std::string toString();

  void Print();

  Curve caf;

  bool create_caf(double x1, double x2, bool ending_burst, uint64_t flow_id);

  uint64_t uid; // unique flow ID (for tracking flows across cuts)
};

}

#endif
