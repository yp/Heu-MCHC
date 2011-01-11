/**
 *
 *
 *                               Heu-MCHC
 *
 * A fast and accurate heuristic algorithm for the haplotype inference
 * problem on pedigree data with recombinations and mutations
 *
 * Copyright (C) 2009,2010,2011  Yuri PIROLA
 *
 * Distributed under the terms of the GNU General Public License (GPL)
 *
 *
 * This file is part of Heu-MCHC.
 *
 * Heu-MCHC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Heu-MCHC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
#ifndef _LOCUS_GRAPH_HPP_
#define _LOCUS_GRAPH_HPP_

#include "data.hpp"
#include <set>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "gen-ped-IO.h"
#include "util.h"


typedef boost::adjacency_list_traits < boost::listS,
													boost::vecS,
													boost::undirectedS >::vertex_descriptor LG_Vertex;

typedef boost::adjacency_list_traits < boost::listS,
													boost::vecS,
													boost::undirectedS >::edge_descriptor LG_Edge;

class LG_VertexProperty {
public:
  LG_VertexProperty(void)
		:individ(NULL), ls(ped::LS_UNDETERMINED), ptilde(-1)
  {}
  pindiv individ;
  ped::locus_status ls;
  int ptilde;

  LG_Vertex pred_vertex;
  LG_Edge pred_edge;
};

class LG_EdgeProperty {
public:
  LG_EdgeProperty(void)
		:d(0), in_spanning_forest(false)
  {}
  int d;
  bool in_spanning_forest;
};

typedef boost::property < boost::graph_name_t, int > Locus;

typedef boost::adjacency_list < boost::listS,
										  boost::vecS,
										  boost::undirectedS,
										  LG_VertexProperty,
										  LG_EdgeProperty,
										  Locus > LocusGraph;

void
build_constraints_from_pedigree(pgenped gp,
										  ped::cycle_constraints_t& cycle_constraints,
										  ped::e_vars_t& m_vars_univ);

#endif
