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
#ifndef _BELIEF_PROPAGATION_HPP_
#define _BELIEF_PROPAGATION_HPP_

#include "data.hpp"
#include "irregular_mat.hpp"

#define ONE_MORE_PROBABLE_THAN( q1t, q1max) ((q1t)<(q1max))

typedef matrix<BP_T> dbl_mat;
typedef vect<BP_T> dbl_vec;
typedef irregular_mat<BP_T> dbl_irregular_mat;

class row_h {
private:
  static size_t cont;
public:
  size_t pos;
  const size_t id;
  ped::cycle_constraint_t* const cc;

  row_h(ped::cycle_constraint_t* const _cc)
		:pos(0), id(cont++), cc(_cc)
  {}
};

class col_h {
private:
  static size_t cont;
public:
  size_t pos;
  const size_t id;
  ped::e_variable_t mv;
  BP_T p1;
  BP_T q1;
  bool phantom;  // Is it a false variable?
  bool assigned; // Useful only for phantom variables

  col_h(const bool _assigned)
		:pos(0), id(cont++), p1(BP0), q1(BP0), phantom(true), assigned(_assigned)
  {
	 if (_assigned)
		p1= BP1;
	 else
		p1= BP0;
  }

  col_h(ped::e_variable_t _mv)
		:pos(0), id(cont++), mv(_mv), p1(BP0), q1(BP0), phantom(false), assigned(false)
  {
  }

  char* to_c_str(void) const {
	 char* str= str_buffer::get_instance().get_buffer();
	 if (phantom)
		sprintf(str, "" ST_FMTL(5) " ) phantom assigned?%s   p1=%.8e  q1=%.8e",
				  id, assigned?"YES":"NO ", BP2D(p1), BP2D(q1));
	 else
		sprintf(str, "" ST_FMTL(5) " ) %s   p1=%.8e  q1=%.8e",
				  id, mv.to_c_str(), BP2D(p1), BP2D(q1));
	 return str;
  }
};

class bp_entry {
public:
  BP_T qm0;
  BP_T rm0;
  bp_entry(void)
		:qm0(BP0), rm0(BP0)
  { }
};

typedef irr_mat<row_h*, col_h*, bp_entry> system_desc;

struct sd_del
  :public unary_function<bp_entry*, void>,
	public unary_function<col_h*, void>,
	public unary_function<row_h*, void>,
	public unary_function<system_desc::entry*, void> {

  void operator()(bp_entry* bp) const {
	 delete bp;
  }

  void operator()(col_h* ch) const {
	 delete ch;
  }

  void operator()(row_h* rh) const {
	 delete rh;
  }

  void operator()(system_desc::entry* e) const {
	 delete e->pdata();
	 delete e;
  }
};



struct _bp_config {

  BP_T gamma;

  size_t max_iter;

  bool ext_debug;

};

typedef struct _bp_config* pbp_config;


void belief_propagation(system_desc& sd, pbp_config bp);






#endif // _BELIEF_PROPAGATION_HPP_
