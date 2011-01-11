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
/*
** conversion-2bob.c
**
** Made by (Yuri)
** Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
**
** Started on  Mon Jul 27 16:01:38 2009 Yuri
** Last update Sun May 12 01:17:25 2002 Speed Blue
*/

#include "gen-ped-IO.h"
#include "util.h"
#include "log-build-info.h"

int main() {
  PRINT_SYSTEM_INFORMATIONS;
  pgenped pg= gp_read_from_file(stdin);
  int* gender= NPALLOC(int, pg->n_indiv);
  for (size_t i= 0; i<pg->n_indiv; ++i) {
	 gender[i]= 1;
  }
  for (size_t i= 0; i<pg->n_indiv; ++i) {
	 if (pg->individuals[i]->fi>=0) {
		gender[pg->individuals[i]->fi]= 1;
	 }
	 if (pg->individuals[i]->mi>=0) {
		gender[pg->individuals[i]->mi]= 2;
	 }
  }
  FILE* fgen= fopen("genotypes.txt", "w");
  if (fgen == NULL) {
	 FATAL("Impossible to open file genotypes.txt for writing!");
	 fail();
  }
  FILE* fped= fopen("pedigree.txt", "w");
  if (fped == NULL) {
	 FATAL("Impossible to open file pedigree.txt for writing!");
	 fail();
  }

  fprintf(fped, "%d\n", pg->n_indiv);
  fprintf(fgen, "%d %d\n", pg->n_indiv, pg->n_loci);
  for (size_t i= 0; i<pg->n_indiv; ++i) {
	 DEBUG("Saving individual %d", i+1);
	 pindiv in= pg->individuals[i];
	 fprintf(fped, "%d\t%c",
				in->id+1, (gender[i]==1)?'M':'F');
	 if (in->fi >= 0) {
		fprintf(fped, "\t%d\t%d\n",
				  in->fi+1, in->mi+1);
	 } else {
		fprintf(fped, "\tN\tN\n");
	 }
	 for (size_t j= 0; j<pg->n_loci-1; ++j) {
		fprintf(fgen, "%d ", in->g[j]);
	 }
	 fprintf(fgen, "%d\n", in->g[pg->n_loci-1]);
  }
  fclose(fgen);
  fclose(fped);
}

