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
#ifndef _BP_DOUBLE_H_
#define _BP_DOUBLE_H_

/**
 * Header defining the type of the values of the
 * Belief Propagation Algorithm using the basic
 * type double (or float)
 **/

#include <cmath>

#ifdef BP_T
#error "Only a single definition for BP_T is admitted in the same program."
#endif

#define BP_T double


// Transform an element into a double
#define BP2D( v ) (v)
// Get the absolute value of an expression
#define BPABS( v ) (fabs(v))
// Condition to establish if a value is not a propability
#define BPINVALID( v ) (((v)>1.0) || ((v)<0.0) || isnan(v))

// Zero, One and Minimum constants
#define BP0 0.0
#define BP1 1.0
#define BP2 2.0
#define BP0_5 0.5
#define BP_MIN 1e-60
#define BP1THRESHOLD (BP1-BP_MIN)

#endif // _BP_DOUBLE_H_
