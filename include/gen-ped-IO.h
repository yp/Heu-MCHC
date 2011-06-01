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
#ifndef _GEN_PED_IO_H_
#define _GEN_PED_IO_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

/**
 *
 * Define the structure of a genotyped pedigree.
 *
 **/

#define GEN0 (0)
#define GEN1 (1)
#define GEN2 (2)


struct _indiv {
  int id;  // individual id
  struct _indiv * f; // father (NULL if founder)
  struct _indiv * m; // mother (NULL if founder)
  int fi; // id of the father (-1 if founder)
  int mi; // id of the mother (-1 if founder)
  int* g; // genotype vector
  int* p; // phase vector
  int* df; // mutation vector from father
  int* dm; // mutation vector from mother
  int* rf; // recombination vector from father
  int* rm; // recombination vector from mother
  int h_f; // haplotype inherithed from the father
  int h_m; // haplotype inherithed from the mother

  int* hv_f; // haplotype vector inherited from the father
  int* hv_m; // haplotype vector inherited from the mother

  unsigned int n_children;
  struct _indiv ** children;
};

typedef struct _indiv* pindiv;

pindiv indiv_create(const int n_loci, const int n_individuals);

void indiv_destroy(pindiv pi);

struct _genped {
  unsigned int n_loci;
  unsigned int n_indiv;
  pindiv* individuals;
};

typedef struct _genped * pgenped;

pgenped gp_create(const int n_loci, const int n_individuals);

void gp_destroy(pgenped gp);

void gp_add_trio(pgenped gp,
					  const unsigned int fi,
					  const unsigned int mi,
					  const unsigned int ci);

pgenped gp_read_from_file(FILE* fin);

pgenped gp_copy(pgenped src);

#ifdef __cplusplus
}
#endif

#endif // _GEN_PED_IO_H_
