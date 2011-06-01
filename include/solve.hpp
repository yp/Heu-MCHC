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
#ifndef _SOLVE_HPP_
#define _SOLVE_HPP_

#include "data.hpp"
#include "gen-ped-IO.h"

class ped_hi_exception {
private:
  std::string msg;
public:
  ped_hi_exception(const std::string& _msg)
		:msg(_msg)
  {}

  const std::string& get_message(void) const {
	 return msg;
  }
};

ped::e_vars_t*
calculate_minimum_solution(pgenped gp, const size_t max_mut, const double gamma,
									const BP_T* p1) throw (ped_hi_exception);





#endif //  _SOLVE_HPP_
