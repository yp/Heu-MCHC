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
 * along with Heu-MCHC.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
#include <string>
#include <iostream>
#include <fstream>

#include <map>

#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graphviz.hpp>

#include <cstdio>

#include "locus-graph.hpp"
#include "irregular_mat.hpp"
#include "belief-propagation.hpp"
#include "util.h"
#include "log.h"
#include <m4ri/m4ri.h>

using namespace boost;
using namespace ped;

#define MZD_INIT( r, c ) mzd_init( r<32?32:r, c<32?32:c )

template <typename T>
class TLN_t {
public:
  TLN_t<T>* const pred;
  T data;
  TLN_t(TLN_t<T>* const _pred, T _data)
		:pred(_pred), data(_data) {}

};

typedef TLN_t<ped::h_variable_t> tln_h_t;
typedef TLN_t<ped::e_variable_t> tln_e_t;


static bool
is_indiv_homoz_at_locus(const pindiv i, const int l) {
  return i!=NULL && i->g[l]!=GEN2;
}

class ped_lg_label_writer {
private:
  LocusGraph& g;
  int l;
public:
  ped_lg_label_writer(LocusGraph& _lg, int _l)
		:g(_lg), l(_l)
  {}

  void operator()(std::ostream& out, const LG_Vertex& v) {
	 out << "[label=\"I" << g[v].individ->id <<
		"^"<< l <<" (ls=" << ped::locus_abbrv_status_names[g[v].ls] << ", p~=" << g[v].ptilde <<")\"]";
  }

  void operator()(std::ostream& out, const LG_Edge& e) {
	 out << "[label=\"d=" << g[e].d << "\"]";
  }
};

class general_locus_graph_visitor: public boost::default_dfs_visitor {
private:
  LocusGraph& lg;

public:
  general_locus_graph_visitor(LocusGraph& _lg)
		:lg(_lg)
  {}

  void initialize_vertex(LG_Vertex v, const LocusGraph& g) {
	 lg[v].pred_vertex= v;
  }

  void tree_edge(LG_Edge e, const LocusGraph& g) {
	 DEBUG("Edge " ST_FMT "-" ST_FMT " is in the spanning forest.",
			 source(e, lg), target(e, lg));
	 lg[e].in_spanning_forest= true;
	 lg[target(e, lg)].pred_vertex= source(e, lg);
	 lg[target(e, lg)].pred_edge= e;
  }
};


static LocusGraph *
build_locus_graph_from_pedigree(pgenped gp,
										  const unsigned int l) {
  my_assert(gp!=NULL);
  assert_ulim(l, gp->n_loci);
  DEBUG("Building locus graph of locus %d", l);

  LocusGraph * plg= new LocusGraph(gp->n_indiv, l);
  LocusGraph& lg= *plg;

  int index= 0;
  BGL_FORALL_VERTICES(v, lg, LocusGraph) {
	 pindiv i= gp->individuals[index];
	 ++index;
	 lg[v].individ= i;
	 lg[v].ptilde= -1;
	 lg[v].ls= ped::LS_UNDETERMINED;
// Set the status of the locus
	 if (is_indiv_homoz_at_locus(i, l)) {
		lg[v].ls= ped::LS_IMMUTABLE;
		lg[v].ptilde= i->g[l];
	 } else {
// Individual heterozygous
		if (is_indiv_homoz_at_locus(i->f, l)) {
// Father homozygous
		  if (is_indiv_homoz_at_locus(i->m, l)) {
// Both parents homozygous
			 if (i->f->g[l]==i->m->g[l]) {
// Same homozygousity -> Ambiguous
				lg[v].ls= ped::LS_PRED_AMBIGUOUS;
				lg[v].ptilde= i->f->g[l];
			 } else {
// Different homozygousity -> Doubly predetermined
				lg[v].ls= ped::LS_PRED_DOUBLE;
				lg[v].ptilde= i->f->g[l];
			 }
		  } else {
// Only father homozygous
			 lg[v].ls= ped::LS_PRED_SINGLE_FATHER;
			 lg[v].ptilde= i->f->g[l];
		  }
		} else if (is_indiv_homoz_at_locus(i->m, l)) {
// Only mother homozygous
		  lg[v].ls= ped::LS_PRED_SINGLE_MOTHER;
		  lg[v].ptilde= (i->m->g[l]+1)%2;
		}
	 }

	 my_assert(lg[v].ptilde!=-1 || lg[v].ls==ped::LS_UNDETERMINED);

  }

  index= 0;
  BGL_FORALL_VERTICES(v, lg, LocusGraph) {
	 pindiv i= gp->individuals[index];
	 if (i->f!=NULL && i->f->g[l]==GEN2) {
		DEBUG("Adding edge %d %d", i->id, i->fi);
		LocusGraph::edge_descriptor ed= add_edge(i->id, i->fi, lg).first;
		lg[ed].d= 0;
	 }
	 if (i->m!=NULL && i->m->g[l]==GEN2) {
		DEBUG("Adding edge %d %d", i->id, i->mi);
		LocusGraph::edge_descriptor ed= add_edge(i->id, i->mi, lg).first;
		lg[ed].d= (gp->individuals[index]->g[l]==GEN2)?1:0;
	 }
	 ++index;
  }
  general_locus_graph_visitor vis(lg);
  depth_first_search(lg, visitor(vis));


#ifdef LOG_GRAPHS
  write_graphviz(std::cout, lg,
					  ped_lg_label_writer(lg, l),
					  ped_lg_label_writer(lg, l));
#endif

  return plg;
}

static LocusGraph**
build_all_locus_graphs_from_pedigree(pgenped gp) {
  my_assert(gp!=NULL);
  LocusGraph** lgs= NPALLOC(LocusGraph*, gp->n_loci);

  for(unsigned int i= 0; i<gp->n_loci; ++i) {
	 lgs[i]= build_locus_graph_from_pedigree(gp, i);
  }

  return lgs;
}

static LocusGraph*
build_general_locus_graph_from_pedigree(pgenped gp) {
// Build a general locus graph (w/o checking homozygousity)
  my_assert(gp!=NULL);
  DEBUG("Building general locus graph");
  LocusGraph* const plg= new LocusGraph(gp->n_indiv, -1);
  LocusGraph& lg= *plg;
  int index= 0;
  BGL_FORALL_VERTICES(v, lg, LocusGraph) {
	 pindiv i= gp->individuals[index];
	 ++index;
	 lg[v].individ= i;
  }

  index= 0;
  BGL_FORALL_VERTICES(v, lg, LocusGraph) {
	 pindiv i= gp->individuals[index];
	 if (i->f!=NULL) {
		LocusGraph::edge_descriptor ed= add_edge(i->id, i->fi, lg).first;
	 }
	 if (i->m!=NULL) {
		LocusGraph::edge_descriptor ed= add_edge(i->id, i->mi, lg).first;
	 }
	 ++index;
  }

  general_locus_graph_visitor vis(lg);
  depth_first_search(lg, visitor(vis));

#ifdef LOG_GRAPHS
  write_graphviz(std::cout, lg,
					  ped_lg_label_writer(lg, -1),
					  ped_lg_label_writer(lg, -1));
#endif

  return plg;
}

static ped::e_variable_t
get_m_variable_from_edge(const LocusGraph& lg,
								 const int u,
								 const int v,
								 const int l) {
// Determino chi e' il genitore
  bool u_child_of_v= (lg[u].individ->fi==v) ||
	 (lg[u].individ->mi==v);

  my_assert(u_child_of_v || (lg[v].individ->fi==u) ||
				(lg[v].individ->mi==u));

  if (u_child_of_v)
	 return ped::e_variable_t(ped::e_variable_t::MUT,
									  u, l, lg[u].individ->fi==v?0:1);
  else
	 return ped::e_variable_t(ped::e_variable_t::MUT,
									  v, l, lg[v].individ->fi==u?0:1);
}

static ped::h_variable_t
get_h_variable_from_edge(const LocusGraph& lg,
								 const int u,
								 const int v,
								 const int l) {
// Determino chi e' il genitore
  bool u_child_of_v= (lg[u].individ->fi==v) ||
	 (lg[u].individ->mi==v);

  my_assert(u_child_of_v || (lg[v].individ->fi==u) ||
				(lg[v].individ->mi==u));

  if (u_child_of_v)
	 return ped::h_variable_t(u, lg[u].individ->fi==v?0:1);
  else
	 return ped::h_variable_t(v, lg[v].individ->fi==u?0:1);
}

static void
get_rec_variables_from_edge_to_tln(const LocusGraph& lg,
											  const int u,
											  const int v,
											  const int l,
											  tln_e_t*& evars,
											  int & cont) {
#ifndef NO_RECOMBINATIONS
// Determino chi e' il genitore
  bool u_child_of_v= (lg[u].individ->fi==v) ||
	 (lg[u].individ->mi==v);

  my_assert(u_child_of_v || (lg[v].individ->fi==u) ||
				(lg[v].individ->mi==u));
  const int p= (u_child_of_v)?v:u;
  const int c= (u_child_of_v)?u:v;
  const int parenthood= lg[c].individ->fi==p?0:1;

  for (int i= 1; i<=l; ++i) {
	 if (lg[p].individ->g[i]==GEN2) {
		evars= new tln_e_t(evars,
								 ped::e_variable_t(ped::e_variable_t::REC,
														 c, i, parenthood));
		++cont;
	 }
  }
#endif
}

static void
get_rec_variables_from_edge_to_set(const LocusGraph& lg,
											  const int u,
											  const int v,
											  const int l,
											  ped::e_vars_t& evars) {
#ifndef NO_RECOMBINATIONS
// Determino chi e' il genitore
  bool u_child_of_v= (lg[u].individ->fi==v) ||
	 (lg[u].individ->mi==v);

  my_assert(u_child_of_v || (lg[v].individ->fi==u) ||
				(lg[v].individ->mi==u));
  const int p= (u_child_of_v)?v:u;
  const int c= (u_child_of_v)?u:v;
  const int parenthood= lg[c].individ->fi==p?0:1;

  for (int i= 1; i<=l; ++i) {
	 if (lg[p].individ->g[i]==GEN2) {
		evars.insert(ped::e_variable_t(ped::e_variable_t::REC,
												 c, i, parenthood));
	 }
  }
#endif
}

static ped::e_variable_t
get_extremal_m_var(const ped::locus_status ls,
						 const int v,
						 const int locus) {
  switch (ls) {
	 case ped::LS_UNDETERMINED:
		fail();
		break;
	 case ped::LS_IMMUTABLE:
// none to add
		return ped::e_variable_t();
		break;
	 case ped::LS_PRED_SINGLE_FATHER:
		return ped::e_variable_t(ped::e_variable_t::MUT, v, locus, 0);
		break;
	 case ped::LS_PRED_SINGLE_MOTHER:
		return ped::e_variable_t(ped::e_variable_t::MUT, v, locus, 1);
		break;
	 case ped::LS_PRED_DOUBLE:
		return ped::e_variable_t(ped::e_variable_t::MUT, v, locus, 0);
		break;
	 case ped::LS_PRED_AMBIGUOUS:
		return ped::e_variable_t(ped::e_variable_t::MUT, v, locus, 0);
		break;
  }
  fail();
  return ped::e_variable_t();
}

template <typename T>
class to_visit_t {
private:
  std::list<T> datas;

public:
  void push(const T& el) {
	 datas.push_front(el);
  }

  T pop(void) {
	 T el= datas.front();
	 datas.pop_front();
	 return el;
  }

  bool empty(void) const {
	 return datas.empty();
  }
};

template <typename SET>
static void
read_tree_path_events(tln_e_t* const n, SET& set) {
  tln_e_t* tmp= n;
  while (tmp->pred!=NULL) {
#ifdef NO_RECOMBINATIONS
	 if (tmp->data.kind() != ped::e_variable_t::REC)
#endif
#ifdef NO_MUTATIONS
	 if (tmp->data.kind() != ped::e_variable_t::MUT)
#endif
		set.insert(tmp->data);
	 tmp= tmp->pred;
  }
}

template <typename SET>
static void
read_tree_path_hvars(tln_h_t* const n, SET& set) {
  tln_h_t* tmp= n;
  while (tmp->pred!=NULL) {
	 set.insert(tmp->data);
	 tmp= tmp->pred;
  }
}

static void
add_constraints_from_v_in_locus(const LocusGraph& lg,
										  const LG_Vertex root,
										  cycle_constraints_t& ccv,
										  const int locus,
										  const int n_indiv,
										  std::vector<bool>& visited) {
  my_assert(lg[root].ls != ped::LS_UNDETERMINED);
  FINETRACE("Adding constraints from predetermined vertex %d.", root);

  to_visit_t< LG_Vertex > to_visit;
  to_visit_t< tln_h_t* > tv_h;
  to_visit_t< tln_e_t* > tv_e;
  to_visit_t< int > tv_ne;
  std::vector< int > sum(n_indiv, -1);

  std::list< tln_h_t* > td_h;
  std::list< tln_e_t* > td_e;
  std::list< int > td_ne;

  ped::e_vars_t::iterator it;
  ped::e_variable_t m;
  ped::h_variable_t h;

  tln_e_t* const root_e= new tln_e_t(NULL, e_variable_t());
  tln_h_t* const root_h= new tln_h_t(NULL, h_variable_t());

  sum[root]= lg[root].ptilde;

  int cont= 1;
  m= get_extremal_m_var(lg[root].ls, root, locus);
  if (m.i>=0) {
	 tln_e_t* n= new tln_e_t(root_e, m);
	 tv_e.push(n);
	 cont= 2;
  } else {
	 tv_e.push(root_e);
  }
  tv_ne.push(cont);
  tv_h.push(root_h);
  to_visit.push(root);

  while (!to_visit.empty()) {
	 LG_Vertex v= to_visit.pop();
	 tln_h_t* hvars= tv_h.pop();
	 tln_e_t* evars= tv_e.pop();
	 int ne= tv_ne.pop();
	 td_h.push_front(hvars);
	 td_e.push_front(evars);
	 td_ne.push_front(ne);
// Iterates over its out edges
	 graph_traits<LocusGraph>::out_edge_iterator e, e_end, next;
	 tie(e, e_end) = out_edges(v, lg);
	 for (; e != e_end; ++e) {
		if (!lg[*e].in_spanning_forest)
		  continue;
		const LG_Vertex u= target(*e, lg);
		if (sum[u]==-1) {
		  FINETRACE("Considering edge %d-%d", v, u);
		  sum[u]= sum[v]+lg[*e].d;
		  m= get_m_variable_from_edge(lg, u, v, locus);
		  tln_e_t* new_evars= new tln_e_t(evars, m);
		  cont= 1;
		  get_rec_variables_from_edge_to_tln(lg, u, v, locus, new_evars, cont);
		  h= get_h_variable_from_edge(lg, u, v, locus);
		  tln_h_t* new_hvars= new tln_h_t(hvars, h);
		  to_visit.push(u);
		  tv_e.push(new_evars);
		  tv_ne.push(cont);
		  tv_h.push(new_hvars);
		  if (lg[u].ls != ped::LS_UNDETERMINED) {
			 DEBUG("Constraint detected between vertex " ST_FMT " and vertex " ST_FMT ".", root, u);
			 visited[u]= true;
			 cycle_constraint_t* cc= new cycle_constraint_t;
			 ccv.push_back(cc);
			 cc->constant= (sum[u]+lg[u].ptilde)%2==1;
			 m= get_extremal_m_var(lg[u].ls, u, locus);
			 if (m.i>=0) {
				new_evars= new tln_e_t(new_evars, m);
				++cont;
			 }
			 read_tree_path_hvars(new_hvars, cc->h_vars);
			 read_tree_path_events(new_evars, cc->events);
		  }
// Free memory
//			 delete new_hvars;
//			 for (int i= 0; i<cont; ++i) {
//				tln_e_t* tmp= new_evars;
//				new_evars= new_evars->pred;
//				delete tmp;
//			 }
//		  }
		}
	 }
  }

  while (!td_h.empty()) {
	 delete td_h.front();
	 td_h.pop_front();
  }
  while (!td_e.empty()) {
	 int ne= td_ne.front();
	 tln_e_t* evars= td_e.front();
	 for (int i= 0; i<ne; ++i) {
		tln_e_t* tmp= evars;
		evars= evars->pred;
		delete tmp;
	 }
	 td_ne.pop_front();
	 td_e.pop_front();
  }
}

static void
add_immutable_immutable_constraints
  (pgenped pg,
	ped::e_vars_t& m_vars_univ,
	ped::cycle_constraints_t& cycle_constraints)
{
  for (unsigned int i= 0; i<pg->n_indiv; ++i) {
	 if (pg->individuals[i]->fi==-1) {
		my_assert(pg->individuals[i]->mi==-1);
	 } else {
		pindiv const ii= pg->individuals[i];
		for (unsigned int j= 0; j<pg->n_loci; ++j) {
		  if (ii->g[j]!=GEN2 && ii->f->g[j]!=GEN2) {
			 if (ii->g[j]!=ii->f->g[j]) {
				ped::e_variable_t mv(ped::e_variable_t::MUT, i, j, 0);
				m_vars_univ.insert(mv);
				ped::cycle_constraint_t* cc= new ped::cycle_constraint_t;
				cc->events.insert(mv);
				cc->constant= true;
				cycle_constraints.push_back(cc);
				DEBUG("Added immutable constraint for indiv. %d at locus %d"
						" because he and his father are homozygous and they have"
						" different alleles.",
						mv.i, mv.l);
			 }
		  }
		  if (ii->g[j]!=GEN2 && ii->m->g[j]!=GEN2) {
			 if (ii->g[j]!=ii->m->g[j]) {
				ped::e_variable_t mv(ped::e_variable_t::MUT, i, j, 1);
				m_vars_univ.insert(mv);
				ped::cycle_constraint_t* cc= new ped::cycle_constraint_t;
				cc->events.insert(mv);
				cc->constant= true;
				cycle_constraints.push_back(cc);
				DEBUG("Added immutable constraint for indiv. %d at locus %d"
						" because he and his mother are homozygous and they have"
						" different alleles.",
						mv.i, mv.l);
			 }
		  }
		  if ((ii->g[j]==GEN2) &&
				(ii->f->g[j] == ii->m->g[j]) &&
				(ii->f->g[j]!=GEN2)) {
			 ped::e_variable_t mv0(ped::e_variable_t::MUT, i, j, 0);
			 if (m_vars_univ.find(mv0)==m_vars_univ.end()) {
				m_vars_univ.insert(mv0);
				ped::e_variable_t mv1(ped::e_variable_t::MUT, i, j, 1);
				m_vars_univ.insert(mv1);
				ped::cycle_constraint_t* cc= new ped::cycle_constraint_t;
				cc->events.insert(mv0);
				cc->events.insert(mv1);
				cc->constant= true;
				cycle_constraints.push_back(cc);
				DEBUG("Added special constraint for indiv. %d at locus %d"
						" because he is ambiguously predetermined.",
						i, j);
			 }
		  }
		}
	 }
  }
}


// Compute the tree path between u and v on lg
static void
tree_path(const LocusGraph& lg,
			 const LG_Vertex& u,
			 const LG_Vertex& v,
			 std::list< LG_Vertex >& path_u_v) {
  list<LG_Edge> path_v_root;
  LG_Vertex tmp= v;
  while (lg[tmp].pred_vertex != tmp) {
	 path_v_root.push_front(lg[tmp].pred_edge);
	 tmp= lg[tmp].pred_vertex;
  }
  list<LG_Edge> path_u_root;
  tmp= u;
  while (lg[tmp].pred_vertex != tmp) {
	 path_u_root.push_front(lg[tmp].pred_edge);
	 tmp= lg[tmp].pred_vertex;
  }
  while (!path_u_root.empty() && !path_v_root.empty()
			&& (path_u_root.front()==path_v_root.front())) {
	 path_u_root.pop_front();
	 path_v_root.pop_front();
  }
  tmp= u;
  path_u_v.push_front(u);
  while (!path_u_root.empty()) {
	 LG_Edge e= path_u_root.back();
	 if (source(e, lg)==tmp) {
		tmp= target(e, lg);
	 } else if (target(e, lg)==tmp) {
		tmp= source(e, lg);
	 } else {
		ERROR("Error in determining the tree path from " ST_FMT " to " ST_FMT ".", u, v);
		fail();
	 }
	 path_u_v.push_back(tmp);
	 path_u_root.pop_back();
  }
  tmp= path_u_v.back();
  while (!path_v_root.empty()) {
	 LG_Edge e= path_v_root.front();
	 if (source(e, lg)==tmp) {
		tmp= target(e, lg);
	 } else if (target(e, lg)==tmp) {
		tmp= source(e, lg);
	 } else {
		ERROR("Error in determining the tree path from " ST_FMT " to " ST_FMT ".", u, v);
		fail();
	 }
	 path_u_v.push_back(tmp);
	 path_v_root.pop_front();
  }
#ifdef LOG_DEBUG_ENABLED
  DEBUG_NOTN("The foundamental cycle of edge (" ST_FMTL(4) " -" ST_FMTL(4) " ) "
				 "contains the path {", u, v);
  for (std::list< LG_Vertex >::const_iterator it= path_u_v.begin();
		 it != path_u_v.end(); ++it) {
	 if (it!=path_u_v.begin())
		DEBUG_MSG("%s",", ");
	 DEBUG_MSG("" ST_FMTL(4) " ", *it);
  }
  DEBUG_MSG("%s","}\n");
#endif
}

static void
add_non_tree_constraints(const LocusGraph& glg,
								 const LocusGraph& lg,
								 cycle_constraints_t& ccv,
								 const int locus) {
  BGL_FORALL_EDGES(e, lg, LocusGraph) {
	 if (lg[e].in_spanning_forest)
		continue;
	 const LG_Vertex u= source(e, lg);
	 const LG_Vertex v= target(e, lg);
	 DEBUG("Adding constraints for non-tree edge " ST_FMT "-" ST_FMT ".", u, v);
	 std::list< LG_Vertex > cycle_e;
	 tree_path(lg, u, v, cycle_e);

	 unsigned int predet_vertices= 0;
	 for (std::list< LG_Vertex >::const_iterator it= cycle_e.begin();
			it != cycle_e.end(); ++it) {
		if (lg[*it].ls!=ped::LS_UNDETERMINED)
		  ++predet_vertices;
	 }
	 DEBUG("The foundamental cycle of (" ST_FMT "-" ST_FMT ") has %u predetermined vertices.",
			 u, v, predet_vertices);
	 if (predet_vertices<=1) {
		DEBUG("The foundamental cycle exists in the locus graph %d.", locus);
		DEBUG("Adding NON-tree constraint.");
		ped::cycle_constraint_t* c= new cycle_constraint_t;
		ccv.push_back(c);
		LG_Vertex prev= v;
		c->constant= false;
		for (std::list< LG_Vertex >::const_iterator it= cycle_e.begin();
			  it != cycle_e.end(); ++it) {
		  bool found= false;
		  LG_Edge lge;
		  tie(lge, found)= edge(prev, *it, lg);
		  my_assert(found);
		  c->constant ^= (lg[lge].d==1);
#ifndef NO_MUTATIONS
		  c->events.insert(get_m_variable_from_edge(lg, prev, *it, locus));
#endif
#ifndef NO_RECOMBINATIONS
		  get_rec_variables_from_edge_to_set(lg, prev, *it, locus, c->events);
#endif
		  c->h_vars.insert(get_h_variable_from_edge(lg, prev, *it, locus));
		  prev= *it;
		}
	 } else {
// Compute the predetermined endpoints of the path
		DEBUG("The foundamental cycle does not exist in the locus graph %d.",
				locus);
		DEBUG("Adding tree constraint.");
		LG_Vertex x= u;
		for (std::list< LG_Vertex >::const_iterator it= cycle_e.begin();
			  (it != cycle_e.end()) && (lg[x].ls==ped::LS_UNDETERMINED);
			  ++it) {
		  x= *it;
		}
		my_assert(lg[x].ls!=ped::LS_UNDETERMINED);
		LG_Vertex y= v;
		for (std::list< LG_Vertex >::const_reverse_iterator it= cycle_e.rbegin();
			  (it != cycle_e.rend()) && (lg[y].ls==ped::LS_UNDETERMINED);
			  ++it) {
		  y= *it;
		}
		my_assert(lg[y].ls!=ped::LS_UNDETERMINED);
		my_assert(x!=y);
		DEBUG("Endpoints of the constraint " ST_FMT "-" ST_FMT ".", x, y);
		ped::cycle_constraint_t* c= new cycle_constraint_t;
		ccv.push_back(c);
		LG_Vertex prev= v;
		c->constant= false;
		bool avoid= false;
		for (std::list< LG_Vertex >::const_iterator it= cycle_e.begin();
			  it != cycle_e.end(); ++it) {
		  if (!avoid) {
			 bool found= false;
			 LG_Edge lge;
			 DEBUG("" ST_FMT " " ST_FMT "", prev, *it);
			 tie(lge, found)= edge(prev, *it, lg);
			 my_assert(found);
			 c->constant ^= (lg[lge].d==1);
#ifndef NO_MUTATIONS
			 c->events.insert(get_m_variable_from_edge(lg, prev, *it, locus));
#endif
#ifndef NO_RECOMBINATIONS
			 get_rec_variables_from_edge_to_set(lg, prev, *it, locus, c->events);
#endif
			 c->h_vars.insert(get_h_variable_from_edge(lg, prev, *it, locus));
		  }
		  prev= *it;
		  if (*it==y)
			 avoid= false;
		  if (*it==x)
			 avoid= true;
		}
#ifndef NO_MUTATIONS
// Add extremal m variables
		ped::e_variable_t m= get_extremal_m_var(lg[x].ls, x, locus);
		if (m.i>=0) {
		  c->events.insert(m);
		}
		m= get_extremal_m_var(lg[y].ls, y, locus);
		if (m.i>=0) {
		  c->events.insert(m);
		}
#endif
// Add predetermined phases
		c->constant ^= ((lg[x].ptilde+lg[y].ptilde)%2)==1;
	 }
  }
}


static void
add_special_constraints(ped::e_vars_t& m_vars_univ,
								LocusGraph** lgs,
								ped::cycle_constraints_t& cycle_constraints)
{
  ped::e_vars_t new_vars;
  ped::e_vars_t::iterator iMVU= m_vars_univ.begin();
  const ped::e_vars_t::iterator& iMVUend= m_vars_univ.end();
  for (; iMVU!=iMVUend; ++iMVU) {
	 const ped::e_variable_t& mv= *iMVU;
	 if (mv.kind()==ped::e_variable_t::MUT && mv.p==0 && (*lgs[mv.l])[mv.i].ls>=ped::LS_PRED_DOUBLE) {
		ped::cycle_constraint_t* cc= new ped::cycle_constraint_t;
		cc->events.insert(mv);
		ped::e_variable_t mv1= mv;
		mv1.p= 1;
		cc->events.insert(mv1);
		cc->constant= (*lgs[mv.l])[mv.i].ls==ped::LS_PRED_AMBIGUOUS;
		cycle_constraints.push_back(cc);
		m_vars_univ.insert(iMVU, mv1);
		new_vars.insert(mv1);
		DEBUG("Added special constraint for indiv. %d at locus %d because he is %s.",
				mv.i, mv.l, ped::locus_status_names[(*lgs[mv.l])[mv.i].ls]);
	 }
  }
}



static void
build_basic_constraints(LocusGraph* const glg,
								LocusGraph* const * const lgs,
								pgenped gp,
								ped::cycle_constraints_t& constraints) {
/*********************************************************
 *
 * Per ogni locus graph Gl
 *   Per ogni vertice v
 *     if v is predetermined then
 *       visita il grafo a partire da v e aggiungi i vincoli
 *
 *********************************************************/
  for (unsigned int l= 0; l<gp->n_loci; ++l) {
	 DEBUG("Starting analysis of constraints of locus %d.", l);
	 LocusGraph& lg= *(lgs[l]);
	 std::vector<bool> visited(gp->n_indiv, false);
	 for (unsigned int v= 0; v<gp->n_indiv; ++v) {
		if (lg[v].ls!=ped::LS_UNDETERMINED && !visited[v]) {
		  visited[v]= true;
		  TRACE("Vertex %d is predetermined.", v);
		  add_constraints_from_v_in_locus(lg, v, constraints,
													 l, gp->n_indiv, visited);
		}
	 }
	 add_non_tree_constraints(*glg, lg, constraints, l);
  }
}

void
build_constraints_from_pedigree(pgenped gp,
										  ped::cycle_constraints_t& final_ct,
										  ped::e_vars_t& e_vars_univ) {

#if (defined NO_RECOMBINATIONS) && (defined NO_MUTATIONS)
#error "Both Recombinations and Mutations are excluded. Impossible to continue."
#endif
#ifdef NO_MUTATIONS
#warning "Only RECOMBINATIONS are considered."
  WARN("Only RECOMBINATIONS are considered.");
#endif
#ifdef NO_RECOMBINATIONS
#warning "Only MUTATIONS are considered."
  WARN("Only MUTATIONS are considered.");
#endif
  my_assert(gp!=NULL);
// Build the locus graphs
  LocusGraph* glg= build_general_locus_graph_from_pedigree(gp);
  LocusGraph** lgs= build_all_locus_graphs_from_pedigree(gp);

// Build the basic constraints
  ped::cycle_constraints_t basic_ct;
  build_basic_constraints(glg, lgs, gp, basic_ct);
#ifdef LOG_DEBUG_ENABLED
  {
	 DEBUG("Basic constraints:");
	 int _cont= 0;
	 for (cycle_constraints_t::iterator it= basic_ct.begin();
			it!=basic_ct.end();
			++it, ++_cont) {
		DEBUG("Constraint %4d.", _cont);
		DEBUG("h_vars: %s", to_c_str((*it)->h_vars));
		DEBUG("events: %s", to_c_str((*it)->events));
		DEBUG("constant: %d", ((*it)->constant)?1:0);
	 }
  }
#endif

// Build the basis of the null-space of the vector space spanned by
// the set of h-variable subsets

// 1- build the map h-variable <->position
  ped::h_vars_t h_vars;
  unsigned int i= 0;
  for (ped::cycle_constraints_t::const_iterator cit= basic_ct.begin();
		 cit!=basic_ct.end(); ++cit) {
	 for (ped::h_vars_t::const_iterator hit= (*cit)->h_vars.begin();
			hit!= (*cit)->h_vars.end(); ++hit){
		h_vars.insert(*hit);
	 }
	 ++i;
  }

  const size_t N_BASIC_CT= i;
  const size_t N_H_VARS= h_vars.size();
#ifdef LOG_STATS_ENABLED
  static bool first= true;
  STATS_IF(first, "# basic constraints " ST_FMT "", N_BASIC_CT);
  STATS_IF(first, "# h-variables " ST_FMT "", N_H_VARS);
#endif

  if (N_BASIC_CT>0) {
	 DEBUG("There are " ST_FMTL(4) "  basic constraints and " ST_FMTL(4) "  h-variables.",
			 N_BASIC_CT, N_H_VARS);

	 std::map<h_variable_t, unsigned int> h_to_pos;
	 std::vector<h_variable_t> pos_to_h(N_H_VARS);
	 i= 0;
	 for (ped::h_vars_t::const_iterator hit= h_vars.begin();
			hit!= h_vars.end(); ++hit, ++i) {
		h_to_pos[*hit]= i;
		pos_to_h.push_back(*hit);
	 }

// 2- build the binary matrix
	 mzd_t *bm= MZD_INIT(N_H_VARS, N_BASIC_CT);
	 i= 0;
	 for (ped::cycle_constraints_t::const_iterator cit= basic_ct.begin();
			cit!=basic_ct.end(); ++cit) {
		for (size_t j= 0; j<N_H_VARS; ++j)
		  mzd_write_bit(bm, j, i, 0);
		for (ped::h_vars_t::const_iterator hit= (*cit)->h_vars.begin();
			  hit!= (*cit)->h_vars.end(); ++hit){
		  mzd_write_bit(bm, h_to_pos[*hit], i, 1);
		}
		++i;
	 }
	 const size_t rank_bm= mzd_echelonize_m4ri(bm, 0, 0);
	 STATS_IF(first, "rank basic constraints " ST_FMT "", rank_bm);
//	 DEBUG("The rank of the basic constraint matrix is %4d.", rank_bm);
/**
 * 3- reorder column so that the first rank_bm are independent
 * */
	 std::vector<int> col_perm;
	 mzd_t *bm_re= MZD_INIT(rank_bm, N_BASIC_CT);
	 size_t pos= 0;
	 for (i= 0; i<rank_bm; ++i) {
		while (mzd_read_bit(bm, i, pos)==0) {
		  ++pos;
		}
// Independent column
		for (size_t j= 0; j<=i; ++j)
		  mzd_write_bit(bm_re, j, i, mzd_read_bit(bm, j, pos));
		col_perm.push_back(pos);
		++pos;
	 }
	 size_t second_part= rank_bm;
	 i= 0;
	 pos= 0;
	 while ((second_part<N_BASIC_CT) && (i<N_H_VARS)) {
		while ((pos<N_BASIC_CT) && (mzd_read_bit(bm, i, pos) == 0)) {
// Move column pos to second_part
		  for (size_t j= 0; j<i; ++j)
			 mzd_write_bit(bm_re, j, second_part, mzd_read_bit(bm, j, pos));
		  col_perm.push_back(pos);
		  ++second_part;
		  ++pos;
		}
		++i;
		++pos;
	 }
// Move remaining columns to second_part
	 while (second_part<N_BASIC_CT) {
		for (size_t j= 0; j<N_H_VARS; ++j)
		  mzd_write_bit(bm_re, j, second_part, mzd_read_bit(bm, j, pos));
		col_perm.push_back(pos);
		++second_part;
		++pos;
	 }

// Transform the independent columns into a identity matrix
	 for (i= rank_bm-1; i>0; --i) {
		for (size_t j= 0; j<i; ++j) {
		  if (mzd_read_bit(bm_re, j, i)==1) {
			 mzd_row_add(bm_re, i, j);
		  }
		}
	 }

// The last N_BASIC_CT - rank_bm columns contain the null-space

// Generate cycle-constraints
	 for (i= rank_bm; i<N_BASIC_CT; ++i) {
		DEBUG("Generating cycle constraint " ST_FMTL(4) " .", i+1-rank_bm);
		ped::cycle_constraint_t* cc= new ped::cycle_constraint_t;
		cc->events= basic_ct[col_perm[i]]->events;
		cc->constant= basic_ct[col_perm[i]]->constant;
		DEBUG("Dependent constraint %d", col_perm[i]);
		for (size_t j= 0; j<rank_bm; ++j) {
		  if (mzd_read_bit(bm_re, j, i)==1) {
			 DEBUG("+ constraint %d", col_perm[j]);
			 inplace_symmetric_difference(cc->events, basic_ct[col_perm[j]]->events);
			 cc->constant= cc->constant != basic_ct[col_perm[j]]->constant;
		  }
		}
		inplace_union(e_vars_univ, cc->events);
		final_ct.push_back(cc);
	 }
	 mzd_free(bm);
	 mzd_free(bm_re);
  } else {
	 INFO("There are no basic constraints.");
  }


  add_special_constraints(e_vars_univ, lgs, final_ct);

  add_immutable_immutable_constraints(gp, e_vars_univ, final_ct);

  for (ped::cycle_constraints_t::const_iterator cit= basic_ct.begin();
		 cit!=basic_ct.end(); ++cit) {
	 delete *cit;
  }
  for (unsigned int l= 0; l<gp->n_loci; ++l)
	 delete (lgs[l]);
  delete glg;
  pfree(lgs);
#ifdef LOG_STATS_ENABLED
  STATS_IF(first, "# initial cycle constraints " ST_FMT "", final_ct.size());
  STATS_IF(first, "# initial events variables " ST_FMT "", e_vars_univ.size());
  first= false;
#endif
}


