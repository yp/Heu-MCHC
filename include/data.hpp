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
#ifndef _DATA_HPP_
#define _DATA_HPP_

#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include <cstring>
#include "util.h"
#include "log.h"

// Only one inclusion is admitted at time
//#include "bp_gmp.h"
#include "bp_double.h"

class str_buffer {
private:
  const size_t size;
  const size_t len;
  size_t i;
  char* const buff;

//Other choices are: size= 61 increment= 11
//Other choices are: size= 107 increment= 14
  str_buffer(void)
		:size(11), len(1024*1024), i(0), buff(new char[size*len])
  {}

  str_buffer(const str_buffer& sb)
		:size(sb.size), len(sb.len), i(sb.i), buff(sb.buff)
  {}

  const str_buffer& operator=(const str_buffer& sb) {
	 return sb;
  }

  ~str_buffer(void) {
	 delete [] buff;
  }

public:
  static str_buffer& get_instance(void) {
	 static str_buffer instance;
	 return instance;
  }

  char* get_buffer(void) {
	 char* r= buff+(i*len);
	 i= (i+3)%size;
	 return r;
  }

};


namespace ped {
  typedef enum _locus_status {
	 LS_UNDETERMINED,
	 LS_IMMUTABLE,
	 LS_PRED_SINGLE_FATHER,
	 LS_PRED_SINGLE_MOTHER,
	 LS_PRED_DOUBLE,
	 LS_PRED_AMBIGUOUS
  } locus_status;

  extern const char* locus_status_names[];
  extern const char* locus_abbrv_status_names[];
  extern const char* kind_names[];

  class e_variable_t;

  std::ostream& operator<<(std::ostream& os, const e_variable_t& m);
  std::ostream& operator<<(std::ostream& os, const e_variable_t* m);


  class e_variable_t {
  private:
	 int kind_;
  public:
	 static const unsigned int N_KINDS= 2;

	 static const unsigned int MUT= 0;
	 static const unsigned int REC= 1;

	 int i;
	 int l;
	 int p;

	 e_variable_t(void)
		  :kind_(-1), i(-1), l(-1), p(-1)
	 {}

	 e_variable_t(unsigned int _kind, int _i, int _l, int _p)
		  :kind_(_kind), i(_i), l(_l), p(_p)
	 {
		my_assert(_kind==MUT || _kind==REC);
		my_assert(_p==0 || _p==1);
	 }

	 unsigned int kind() const {
		my_assert(kind_!=-1);
		return kind_;
	 }

	 bool operator==(const e_variable_t& m) const {
		return (kind_==m.kind_) && (i==m.i) && (l==m.l) && (p==m.p);
	 }

	 bool operator<(const e_variable_t& m) const {
		return ((kind_<m.kind_) ||
				  ((kind_==m.kind_) && (i<m.i)) ||
				  ((kind_==m.kind_) && (i == m.i) && (l < m.l)) ||
				  ((kind_==m.kind_) && (i == m.i) && (l == m.l) && (p > m.p)));
	 }

	 char* to_c_str(void) const {
		char* str= str_buffer::get_instance().get_buffer();
		sprintf(str, "([%s] %4d, %3d, %1d)", kind_names[kind_], i, l, p);
		return str;
	 }

  };

  typedef struct {

	 bool operator()(const e_variable_t& m1, const e_variable_t& m2) const {
		return m1<m2;
	 }

	 bool operator()(const e_variable_t* m1, const e_variable_t* m2) const {
		return *m1 < *m2;
	 }
  } e_variable_cmp;

  class h_variable_t {
  public:

	 int i;
	 int p;

	 h_variable_t(void)
		  :i(-1), p(-1)
	 {}

	 h_variable_t(int _i, int _p)
		  :i(_i), p(_p)
	 {
		my_assert(_p==0 || _p==1);
	 }

	 bool operator==(const h_variable_t& h) const {
		return (i==h.i) && (p==h.p);
	 }

	 bool operator<(const h_variable_t& h) const {
		return ((i<h.i) ||
				  ((i == h.i) && (p > h.p)));
	 }

	 char* to_c_str(void) const {
		char* str= str_buffer::get_instance().get_buffer();
		sprintf(str, "(H %4d, %1d)", i, p);
		return str;
	 }

  };

  typedef struct {

	 bool operator()(const h_variable_t& h1, const h_variable_t& h2) const {
		return h1<h2;
	 }

	 bool operator()(const h_variable_t* h1, const h_variable_t* h2) const {
		return *h1 < *h2;
	 }
  } h_variable_cmp;


  typedef std::set<h_variable_t, h_variable_cmp> h_vars_t;
  typedef std::set<e_variable_t, e_variable_cmp> e_vars_t;

  std::ostream& operator<<(std::ostream& os, const e_vars_t& mvars);

  template <typename SET>
  const char* to_c_str(const SET& evars) {
	 static char* str= new char[1000000];
	 bool first= true;
	 sprintf(str, "{");
	 for (typename SET::const_iterator eit= evars.begin();
			eit!= evars.end();
			++eit) {
		if (!first)
		  strcat(str, ", ");
		strcat(str, (*eit).to_c_str());
		first= false;
	 }
	 strcat(str, "}");
	 return str;
  }

  class cycle_constraint_t {
  public:
	 h_vars_t h_vars;
	 e_vars_t events;
	 bool constant;

	 cycle_constraint_t()
		  : h_vars(), events(), constant(false)
	 {}

	 cycle_constraint_t(const h_vars_t& _h_vars,
							  const e_vars_t& _events,
							  const bool _constant)
		  :h_vars(_h_vars), events(_events), constant(_constant)
	 {}

  };

  typedef std::vector< cycle_constraint_t* > cycle_constraints_t;

  template <typename set_t>
  void inplace_union(set_t& A,
							const set_t& B) {
	 typename set_t::iterator iA= A.begin();
	 typename set_t::iterator iB= B.begin();
	 const typename set_t::const_iterator& iAend= A.end();
	 const typename set_t::const_iterator& iBend= B.end();
	 while (iA!=iAend && iB!=iBend) {
		if (*iA == *iB) {
		  ++iA;
		  ++iB;
		} else {
		  if (*iA < *iB) {
			 ++iA;
		  } else if (*iB < *iA) {
			 A.insert(iA, *iB);
			 ++iB;
		  }
		}
	 }
	 while (iB!=iBend) {
		iA= A.insert(iA, *iB);
		++iB;
	 }
  }

  template <class set_t>
  void inplace_symmetric_difference(set_t& A,
							const set_t& B) {
	 typename set_t::iterator iA= A.begin();
	 typename set_t::iterator iB= B.begin();
	 const typename set_t::const_iterator& iAend= A.end();
	 const typename set_t::const_iterator& iBend= B.end();
	 while (iA!=iAend && iB!=iBend) {
		if (*iA == *iB) {
		  typename set_t::iterator iTMP= iA;
		  ++iA;
		  A.erase(iTMP);
		  ++iB;
		} else {
		  if (*iA < *iB) {
			 ++iA;
		  } else if (*iB < *iA) {
			 A.insert(iA, *iB);
			 ++iB;
		  }
		}
	 }
	 while (iB!=iBend) {
		iA= A.insert(iA, *iB);
		++iB;
	 }
  }


}

#endif // _DATA_HPP_
