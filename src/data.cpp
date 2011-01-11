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
#include "data.hpp"

namespace ped {

  const char* locus_status_names[] = {
	 "undetermined",
	 "immutable",
	 "single predetermined from father",
	 "single predetermined from mother",
	 "double predetermined",
	 "ambiguously predetermined"
  };

  const char* locus_abbrv_status_names[] = {
	 "?",
	 "I",
	 "SF",
	 "SM",
	 "D",
	 "A"
  };

  const char* kind_names[] = {
	 "MUT",
	 "REC"
  };


  std::ostream& operator<<(std::ostream& os, const e_variable_t& e) {
	 return os << &e;
  }

  std::ostream& operator<<(std::ostream& os, const e_variable_t* m) {
	 if (m==NULL)
		return os;
	 else
		return os << "(" << m->i << ","<<m->l<<","<<m->p<<")";
  }

  std::ostream& operator<<(std::ostream& os, const e_vars_t& mvars) {
	 bool first= true;
	 os << "{";
	 for (ped::e_vars_t::iterator mit= mvars.begin();
			mit!= mvars.end();
			++mit) {
		if (!first)
		  os << ", ";
		os << *mit;
		first= false;
	 }
	 os << "}";
	 return os;
  }



}
