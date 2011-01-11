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
#include "solve.hpp"

#include "util.h"
#include "log.h"
#include "belief-propagation.hpp"
#include "locus-graph.hpp"

#ifndef MAX_ITER
#define MAX_ITER 100
#endif

#include "my_time.h"
extern pmytime pt_gauss;

#include <m4ri/m4ri.h>


using namespace ped;

static void
print_mapping_variables(FILE* fout, system_desc& sd) {
  fprintf(fout, "* Column headers\n");
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end();
		 cn= cn->next()) {
	 fprintf(fout, "* %s\n", (**cn).header()->to_c_str());
  }
}

static void
print_system_description(FILE* fout, system_desc& sd, bool force= false) {
  if (!force && sd.n_cols()-sd.n_rows()>300) {fprintf(fout, "* System description too big. Skipped.\n");return;}
  print_mapping_variables(fout, sd);
  fprintf(fout, "* System description\n");
  // fprintf(fout, "* ");
  // for (system_desc::cols_node* cn= sd.cols_begin();
  // 		 cn != sd.cols_end();
  // 		 cn= cn->next()) {
  // 	 fprintf(fout, "%5d ", (**cn).header()->id);}
  // fprintf(fout, "\n");
  for (system_desc::rows_node* rn= sd.rows_begin();
		 rn != sd.rows_end();
		 rn= rn->next()) {
	 fprintf(fout, "* " ST_FMTL(5) " ", (**rn).header()->id);
	 system_desc::cols_node* cn= sd.cols_begin();
	 system_desc::row_node* ri= (**rn).begin();
	 while (ri != (**rn).end()) {
		if ((**ri).col().header()->id == (**cn).header()->id) {
		  fprintf(fout, "X");
		  cn= cn->next();
		  ri= ri->next();
		} else {
		  fprintf(fout, " ");
		  cn= cn->next();
		}
	 }
	 fprintf(fout, "\n");
  }
}

static system_desc::cols_node*
apply_bp(system_desc& sd, const size_t iter,
			const BP_T* p,
			const double gamma,
			const size_t max_iter=100) {
  INFO("Iteration of belief propagation.");
  DEBUG("Set the priori probabilities.");
  for (system_desc::cols_node* cs= sd.cols_begin();
		 cs != sd.cols_end() &&
			!(**cs).header()->phantom;
		 cs= cs->next()) {
	 (**cs).header()->p1= p[(**cs).header()->mv.kind()];
  }
  struct _bp_config bp = {
	 gamma,
	 max_iter,
	 false //true //iter==10 //false
  };
  belief_propagation(sd, &bp);
// set to 1 one of the variables that has the highest q1 probability
  system_desc::cols_node* max1= NULL;
  BP_T q1max= BP0;
  size_t n_max= 0;
  for (system_desc::cols_node* cs= sd.cols_begin();
		 cs != sd.cols_end() &&
			!(**cs).header()->phantom;
		 cs= cs->next()) {
	 const BP_T q1t= (**cs).header()->q1;
//	 TRACE("Variable %3d with q1 %23.20f.", (**cs).header()->id, q1t);
	 if((max1==NULL) || ONE_MORE_PROBABLE_THAN(q1t, q1max)) {
		max1= cs;
		q1max= q1t;
		n_max= 1;
	 } else if (q1t == q1max) {
		++n_max;
	 }
  }
  INFO("The maximum posteriori probability is %.8e.", BP2D(q1max));
  INFO("Variables with posteriori probability %.8e are:", BP2D(q1max));
  max1= NULL;
  for (system_desc::cols_node* cs= sd.cols_begin();
		 cs != sd.cols_end();
		 cs= cs->next()) {
	 if (!(**cs).header()->phantom) {
		const BP_T q1t= (**cs).header()->q1;
		if (q1t==q1max) {
		  INFO("%s", (**cs).header()->to_c_str());
		  if ((max1==NULL) && (rnd_01() < (1.0/n_max))) {
			 max1= cs;
		  }
		  --n_max;
		}
	 }
  }
  my_assert(max1 != NULL);
  INFO("Setting variable " ST_FMTL(5) "  %s to 1.",
		 (**max1).header()->id, (**max1).header()->mv.to_c_str());
  return max1;
}


// check if the remaining constraints are all satisfied
static bool
is_system_feasible(const system_desc& sd) {
// if the panthom variables are not assigned, then the system is feasible
  bool feasible= true;
  for (system_desc::rows_node* rs= sd.rows_begin();
		 feasible && rs != sd.rows_end();
		 rs= rs->next()) {
// The last element of each row has to be a phantom variable
	 my_assert((**rs).end()->prev()->data().col().header()->phantom);
	 feasible= !(**(**rs).end()->prev()).col().header()->assigned;
  }
  return feasible;
}

static system_desc::cols_node*
remove_column(system_desc::cols_node* csn,
				  bool change) {
  DEBUG("Removing variable " ST_FMTL(5) " .", (**csn).header()->id);
  for (system_desc::col_node* cn= (**csn).begin();
		 cn != (**csn).end();
		 ){
// If requested, change the fixed value of the phantom variable
	 if (change) {
		DEBUG("Changing assignment of constraint " ST_FMTL(5) " ", (**cn).row().header()->id);
		col_h* const chp= (**cn).row().end()->prev()->data().col().header();
		my_assert(chp->phantom);
		chp->assigned= !chp->assigned;
		chp->p1= chp->assigned?BP1:BP0;
	 }
	 system_desc::row_node* const rn= (**cn).row_n();
	 rn->remove();
	 cn= cn->free_and_remove(sd_del());
  }
  delete (**csn).header();
  delete csn->pdata();
  return csn->remove();
}

static system_desc::rows_node*
remove_row(system_desc::rows_node* rsn) {
  DEBUG("Removing constraint " ST_FMTL(5) " .", (**rsn).header()->id);
  for (system_desc::row_node* rn= (**rsn).begin();
		 rn != (**rsn).end();
		 ){
	 system_desc::col_node* const cn= (**rn).col_n();
	 cn->remove();
	 rn= rn->free_and_remove(sd_del());
  }
  delete (**rsn).header();
  delete rsn->pdata();
  return rsn->remove();
}

static bool
simplify_empty_columns(system_desc& sd) {
  bool mod= false;
  for (system_desc::cols_node* csn= sd.cols_begin();
		 csn != sd.cols_end();
		 ) {
	 if ((**csn).size()==0) {
		DEBUG("it header %p", (void*)((**csn).header()));
		DEBUG("Variable [%s] does not appear in any constraints.",
				(**csn).header()->to_c_str());
		csn= remove_column(csn, false);
		mod= true;
	 } else {
		csn= csn->next();
	 }
  }
  return mod;
}

static inline unsigned int
count_ones(const unsigned int v) {
  unsigned int c; // store the total here
  static const int S[] = {1, 2, 4, 8, 16}; // Magic Binary Numbers
  static const int B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF, 0x0000FFFF};

  c = v - ((v >> 1) & B[0]);
  c = ((c >> S[1]) & B[1]) + (c & B[1]);
  c = ((c >> S[2]) + c) & B[2];
  c = ((c >> S[3]) + c) & B[3];
  c = ((c >> S[4]) + c) & B[4];
  return c;
}

#ifdef _SIMPLIFY_CYCLE_CONSTRAINTS
static void
simplify_cycle_constraints(cycle_constraints_t& cc,
											 cycle_constraints_t& ncc) {
// Build the variable set
  e_vars_t ev;
  size_t i= 0;
  for (cycle_constraints_t::iterator it= cc.begin();
		 it!=cc.end();
		 ++it, ++i) {
	 ev.insert(((*it)->events).begin(), ((*it)->events).end());
  }
  std::vector<e_variable_t> evv(ev.begin(), ev.end());

  const size_t n_cc= i;
  const size_t n_var= evv.size();

// Shuffle the variable set
  for (i= 0; i<n_var; ++i) {
	 const size_t newv= rand()%n_var;
	 e_variable_t _t= evv[i];
	 evv[i]= evv[newv];
	 evv[newv]= _t;
  }

  i= 0;
  std::map<e_variable_t, int> evm;
  for (std::vector<e_variable_t>::const_iterator it= evv.begin();
		 it!= evv.end(); ++it, ++i) {
	 evm[*it]= i;
  }


// Build the binary matrix
  mzd_t* bm= mzd_init(n_cc, n_var+1);
  size_t r= 0;
  for (cycle_constraints_t::iterator it= cc.begin();
		 it!=cc.end();
		 ++it, ++r) {
	 for (e_vars_t::const_iterator ite= ((*it)->events).begin();
			ite!=((*it)->events).end(); ++ite) {
		mzd_write_bit(bm, r, evm[*ite], 1);
	 }
	 mzd_write_bit(bm, r, n_var, ((*it)->constant)?1:0);
  }

  const size_t rank= mzd_echelonize_m4ri(bm, 1, 0);
  for (r= 0; r<rank; ++r) {
	 ped::cycle_constraint_t* c= new ped::cycle_constraint_t;
	 ncc.push_back(c);
	 for (i= r; i<n_var; ++i) {
		if (mzd_read_bit(bm, r, i)==1) {
		  c->events.insert(evv[i]);
		}
	 }
	 c->constant= (mzd_read_bit(bm, r, n_var)==1);
  }

  mzd_free(bm);
}
#endif

static bool
full_simplify_constraints_gauss(system_desc& sd, e_vars_t* const mmut)
  throw (ped_hi_exception) {
  MYTIME_start(pt_gauss);
  bool mod= false;
  const size_t num_constraints= sd.n_rows();
  const size_t num_vars= sd.n_cols()-sd.n_rows();
  INFO("Full Gauss elimination.");
  mzd_t *bm_c = mzd_init(num_constraints, num_vars+1);

  size_t c= 0;
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end() && !(**cn).header()->phantom;
		 cn= cn->next()){
	 (**cn).header()->pos= c;
	 ++c;
  }

  size_t i= 0;
  size_t n_ones= 0;
  for (system_desc::rows_node* rs= sd.rows_begin();
		 rs != sd.rows_end();
		 rs= rs->next(), ++i) {
	 for (system_desc::row_node* rn= (**rs).begin();
			rn != (**rs).end() && !(**rn).col().header()->phantom;
			rn= rn->next()) {
		mzd_write_bit(bm_c, i, (**rn).col().header()->pos, 1);
		++n_ones;
	 }
	 system_desc::row_node* pn= (**rs).end()->prev();
	 my_assert((**pn).col().header()->phantom);
	 mzd_write_bit(bm_c, i, num_vars, (**pn).col().header()->assigned?1:0);
  }
  for (size_t r= 0; r<num_constraints; ++r) {
	 const size_t newr= rand()%num_constraints;
	 mzd_row_swap(bm_c, r, newr);
  }
  INFO("Num. of one entries:     " ST_FMTL(11) " ", n_ones);
  const size_t rank_c= mzd_echelonize_m4ri(bm_c, 1, 2);
  INFO("Rank of the complete matrix:" ST_FMTL(8) " ", rank_c);
  size_t first_one= rank_c-1;
  while ((first_one <= num_vars) && (mzd_read_bit(bm_c, rank_c-1, first_one)==0))
	 ++first_one;
  if (first_one==num_vars)
	 throw ped_hi_exception("Unsatisfable configuration found (from FullGauss).");
// Erasing variables that must assume a given value
// Such variables have only one 1 entry on a row
#ifdef LOG_INFO_ENABLED
  size_t erased_vars= 0;
  size_t set_vars= 0;
  size_t erased_constr= 0;
#endif
  const size_t lim32= num_vars-(num_vars%32);
  size_t totones= 0;
  for (i= 0; i<rank_c; ++i) {
// Check if row i has only one 1 entry in bm
	 size_t n_ones_in_row= 0;
	 unsigned int wrd;
	 for (size_t c= (i/32)*32; (c<lim32)&&(n_ones_in_row<2); c+= 32) {
		wrd= mzd_read_bits(bm_c, i, c, 32);
		n_ones_in_row+= count_ones(wrd);
	 }
	 if (lim32<num_vars) {
		wrd= mzd_read_bits(bm_c, i, lim32, (num_vars%32));
		n_ones_in_row+= count_ones(wrd);
	 }
	 totones+= n_ones_in_row;
	 TRACE("Constraint %5d has " ST_FMTL(5) "  1-entries in the row-echelon form.",
			 bm_c->rowperm[i], n_ones_in_row);
	 my_assert(n_ones_in_row>0);
	 if (n_ones_in_row==1) {
		DEBUG("Constraint " ST_FMTL(5) "  has only one 1-entry in the row-echelon form.",
			  bm_c->rowperm[i]);
// Find the only one 1-entry
		size_t c=0;
		while  (c<num_vars && mzd_read_bit(bm_c, i, c)==0)
		  ++c;
		system_desc::cols_node* csn= sd.cols_begin();
		while ((**csn).header()->pos!=c) {
		  csn= csn->next();
		}
		bool change= mzd_read_bit(bm_c, i, num_vars)==1;
		ped::e_variable_t mv= (**csn).header()->mv;
#ifdef LOG_INFO_ENABLED
		++erased_vars;
#endif
		if (change) {
		  INFO("Adding variable %s to the result because of FullGaussSimplification.",
				 mv.to_c_str());
		  mmut->insert(mv);
#ifdef LOG_INFO_ENABLED
		  ++set_vars;
#endif
		} else {
		  INFO("Removing variable %s that must be unset because of FullGaussSimplification.",
				 mv.to_c_str());
		}
		remove_column(csn, change);
		mod= true;
	 }
  }
  INFO("Erased " ST_FMTL(4) "  variables (" ST_FMTL(4) "  of them added as predicted events)",
		 erased_vars, set_vars);
  if (rank_c<num_constraints) {
	 std::vector<size_t> toerase;
	 for (i=rank_c; i<num_constraints; ++i) {
		DEBUG("Erasing constraint " ST_FMTL(5) "  because it is linear dependent.", bm_c->rowperm[i]);
		toerase.push_back(bm_c->rowperm[i]);
	 }
	 std::sort(toerase.begin(), toerase.end());
	 i= 0;
	 std::vector<size_t>::const_iterator it= toerase.begin();
	 for (system_desc::rows_node* rsn= sd.rows_begin();
			rsn!= sd.rows_end() && it != toerase.end();
			++i) {
		if (i == *it) {
		  DEBUG("Erasing constraint " ST_FMTL(5) " .", *it);
		  system_desc::rows_node* rsn1= rsn->next();
		  mod= true;
		  remove_row(rsn);
		  rsn= rsn1;
		  ++it;
#ifdef LOG_INFO_ENABLED
		  ++erased_constr;
#endif
		} else {
		  rsn= rsn->next();
		}
	 }
	 INFO("Erased " ST_FMTL(4) "  constraints.", erased_constr);
  }
  mzd_free(bm_c);
  INFO("Gauss elimination terminated.");
  MYTIME_stop(pt_gauss);
  return mod;
}

// apply the gauss elimination, the constraint elimination,
// and the variable elimination
static void
simplify_constraints(system_desc& sd, e_vars_t* const mmut) {
  INFO("Simplification of the system.");
  const size_t num_constraints= sd.n_rows();
  const size_t num_vars= sd.n_cols()-sd.n_rows();
  DEBUG("Num. of remaining constraints: " ST_FMTL(5) " ", num_constraints);
  DEBUG("Num. of remaining variables: " ST_FMTL(7) " ", num_vars);
  full_simplify_constraints_gauss(sd, mmut);
  simplify_empty_columns(sd);
#ifdef LOG_TRACE_ENABLED
  DEBUG("Simplification terminated.");
  print_system_description(stderr, sd);
#endif
}



static void
build_system_description(system_desc& sd,
								 cycle_constraints_t& cycle_constraints) {
  for (cycle_constraints_t::const_iterator it= cycle_constraints.begin();
		 it!=cycle_constraints.end();
		 ++it) {
	 const e_vars_t& cc= (*it)->events;
	 system_desc::row_t& r= sd.append_row(new row_h(*it))->data();
	 e_vars_t::iterator ccmvit= cc.begin();
	 system_desc::cols_node* cs= sd.cols_begin();
	 while (cs != sd.cols_end() && ccmvit != cc.end()) {
		if ((**cs).header()->mv < *ccmvit) {
		  FINETRACE("Comparing variables %s to %s: SMALLER.",
				cs->data().header()->mv.to_c_str(), (*ccmvit).to_c_str());
		  cs= cs->next();
		} else if (*ccmvit < (**cs).header()->mv) {
		  FINETRACE("Comparing variables %s to %s: LARGER.",
				cs->data().header()->mv.to_c_str(), (*ccmvit).to_c_str());
		  FINETRACE("Insert variable %s as column.", (*ccmvit).to_c_str());
		  system_desc::col_t& c=
			 sd.insert_col_before(new col_h(*ccmvit), cs)->data();
		  sd.insert_at(r, c, new bp_entry);
		  ++ccmvit;
		} else if (*ccmvit == (**cs).header()->mv) {
		  FINETRACE("Comparing variables %s to %s: EQUAL.",
				cs->data().header()->mv.to_c_str(), (*ccmvit).to_c_str());
		  sd.insert_at(r, cs->data(), new bp_entry);
		  cs= cs->next();
		  ++ccmvit;
		}
	 }
	 while (ccmvit != cc.end()) {
		FINETRACE("Insert variable %s as column.", (*ccmvit).to_c_str());
		system_desc::col_t& c= sd.append_col(new col_h(*ccmvit))->data();
		sd.insert_at(r, c, new bp_entry);
		++ccmvit;
	 }
  }
// Add the phantom variables
  system_desc::rows_node* rn= sd.rows_begin();
  for (cycle_constraints_t::const_iterator it= cycle_constraints.begin();
		 it!=cycle_constraints.end();
		 ++it, rn= rn->next()) {
	 my_assert(rn != sd.rows_end());
	 system_desc::col_t& c= sd.append_col(new col_h((*it)->constant))->data();
	 sd.insert_at(rn->data(), c, new bp_entry);
  }
#ifdef LOG_TRACE_ENABLED
  print_system_description(stderr, sd);
#endif
}

e_vars_t*
calculate_minimum_solution(pgenped gp,
									const size_t max_mut,
									const double gamma,
									const BP_T* p1
									) throw (ped_hi_exception){
  cycle_constraints_t cycle_constraints;
  e_vars_t m_vars_univ;
  build_constraints_from_pedigree(gp, cycle_constraints,
											 m_vars_univ);
#ifdef LOG_DEBUG_ENABLED
  DEBUG("Constraints:");
  int i= 0;
  for (cycle_constraints_t::iterator it= cycle_constraints.begin();
		 it!=cycle_constraints.end();
		 ++it, ++i) {
	 DEBUG("Constraint %4d.", i);
	 DEBUG("events: %s", to_c_str((*it)->events));
	 DEBUG("constant: %d", ((*it)->constant)?1:0);
  }
#endif

  system_desc sd;
#ifdef _SIMPLIFY_CYCLE_CONSTRAINTS
  cycle_constraints_t new_cycle_constraints;
  simplify_cycle_constraints(cycle_constraints, new_cycle_constraints);
#ifdef LOG_DEBUG_ENABLED
  DEBUG("Simplified constraints:");
  i= 0;
  for (cycle_constraints_t::iterator it= new_cycle_constraints.begin();
		 it!=new_cycle_constraints.end();
		 ++it, ++i) {
	 DEBUG("Constraint %4d.", i);
	 DEBUG("events: %s", to_c_str((*it)->events));
	 DEBUG("constant: %d", ((*it)->constant)?1:0);
  }
#endif
// Build the system description
  build_system_description(sd, new_cycle_constraints);

#else // IF NOT _SIMPLIFY_CYCLE_CONSTRAINTS

// Build the system description
  build_system_description(sd, cycle_constraints);

#endif

  e_vars_t* ris= new e_vars_t;
  full_simplify_constraints_gauss(sd, ris);
#ifdef LOG_DEBUG_ENABLED
  print_mapping_variables(stderr, sd);
#endif

// Apply the BP algorithm
  size_t iter= 0;
  size_t max_iter= MAX_ITER;
#ifdef LOG_INFO_ENABLED
  for (unsigned int i= 0; i<ped::e_variable_t::N_KINDS; ++i) {
	 INFO("Base probability %s= %.4e.", ped::kind_names[i], BP2D(p1[i]));
  }
#endif
  INFO("Maximum no. of iterations= " ST_FMTL(4) " .", max_iter);
  INFO("Damping gamma= %.5e.", gamma);
  INFO("No. of constraints= " ST_FMTL(6) " .", sd.n_rows());
  INFO("No. of variables=  " ST_FMTL(7) " .", sd.n_cols()-sd.n_rows());
  try {
	 while (!is_system_feasible(sd)){
		INFO("Iteration " ST_FMTL(4) " .", iter+1);
		simplify_constraints(sd, ris);
		if (!is_system_feasible(sd)) {
		  system_desc::cols_node* mut_col= apply_bp(sd, iter,
																  p1,
																  gamma,
																  max_iter);
		  ris->insert((**mut_col).header()->mv);
		  remove_column(mut_col, true);
		}
		++iter;
		if (ris->size()>max_mut)
		  throw ped_hi_exception("Maximum mutation limit passed.");
	 }
//Destroy the system and the cycle constraints
#ifdef _SIMPLIFY_CYCLE_CONSTRAINTS
	 for (cycle_constraints_t::iterator it= new_cycle_constraints.begin();
			it!=new_cycle_constraints.end();
			++it) {
		delete *it;
	 }
#endif
	 for (cycle_constraints_t::iterator it= cycle_constraints.begin();
			it!=cycle_constraints.end();
			++it) {
		delete *it;
	 }
	 sd.free(sd_del());
  }
  catch (int tmp) {
	 for (cycle_constraints_t::iterator it= cycle_constraints.begin();
			it!=cycle_constraints.end();
			++it) {
		delete *it;
	 }
	 sd.free(sd_del());
	 delete ris;
	 throw;
  }

  return ris;
}
