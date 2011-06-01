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
/*
** conversion-2simwalk.c
**
** Made by (Yuri Pirola)
** Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
**
** Started on  Wed Mar 17 12:18:27 2010 Yuri Pirola
** Last update Sun May 12 01:17:25 2002 Speed Blue
*/

#include "gen-ped-IO.h"
#include "util.h"
#include "log-build-info.h"

int main() {
  PRINT_SYSTEM_INFORMATIONS;
  pgenped pg= gp_read_from_file(stdin);
  char* gender= NPALLOC(char, pg->n_indiv);
  for (size_t i= 0; i<pg->n_indiv; ++i) {
	 gender[i]= 'M';
  }
  for (size_t i= 0; i<pg->n_indiv; ++i) {
	 if (pg->individuals[i]->fi>=0) {
		gender[pg->individuals[i]->fi]= 'M';
	 }
	 if (pg->individuals[i]->mi>=0) {
		gender[pg->individuals[i]->mi]= 'F';
	 }
  }
// Save MAP file
  INFO("Saving file MAP.DAT...");
  FILE* MAP= fopen("MAP.DAT", "w");
  fail_if_msg(MAP == NULL,
				  "Impossible to write to file MAP.DAT");
  for (unsigned int j= 0; j<pg->n_loci; ++j) {
	 fprintf(MAP, "L%u\n        0.01000\n", j);
  }
  fclose(MAP);
// Save LOCUS file
  INFO("Saving file LOCUS.DAT...");
  FILE* LOCUS= fopen("LOCUS.DAT", "w");
  fail_if_msg(LOCUS == NULL,
				  "Impossible to write to file LOCUS.DAT");
  for (unsigned int j= 0; j<pg->n_loci; ++j) {
	 fprintf(LOCUS,
				"L%-7uAUTOSOME 2 0\n"
				"       1 0.50000\n"
				"       2 0.50000\n", j);
  }
  fclose(LOCUS);

// Save PEDIGREE file
  INFO("Saving file PEDIGREE.DAT...");
  FILE* PEDIGREE= fopen("PEDIGREE.DAT", "w");
  fail_if_msg(PEDIGREE == NULL,
				  "Impossible to write to file PEDIGREE.DAT");
  fprintf(PEDIGREE,
			 "(I5,1X,A8)\n"
			 "(3A8,2A1,1X,%u(A3,1X))\n",pg->n_loci);
  fprintf(PEDIGREE,
			 "%-5u PEDIGREE\n",pg->n_indiv);
  for (unsigned int i= 0; i<pg->n_indiv; ++i) {
	 pindiv in= pg->individuals[i];
	 fprintf(PEDIGREE,
				"I%-7u", in->id+1);
	 if (in->fi<0)
		fprintf(PEDIGREE, "        ");
	 else
		fprintf(PEDIGREE, "I%-7u", in->fi+1);
	 if (in->mi<0)
		fprintf(PEDIGREE, "        ");
	 else
		fprintf(PEDIGREE, "I%-7u", in->mi+1);
	 fprintf(PEDIGREE, "%c  ", gender[i]);
	 for (unsigned int j= 0; j<pg->n_loci; ++j) {
		if (in->g[j]==2) {
		  fprintf(PEDIGREE, "1/2 ");
		} else {
		  fprintf(PEDIGREE, "%d/%d ", in->g[j]+1, in->g[j]+1);
		}
	 }
	 fprintf(PEDIGREE, "\n");
  }
  fclose(PEDIGREE);

  INFO("Saving file BATCH2.DAT...");
  FILE* BATCH2= fopen("BATCH2.DAT", "w");
  fail_if_msg(BATCH2 == NULL,
				  "Impossible to write to file BATCH2.DAT");
  fprintf(BATCH2,
			 "\n"
			 "01                       ! batch item number\n"
			 "1                        ! analysis: 1=Haplo; 2=LOD; 3=NPL; 4=IBD; 5=Mistyping\n"
			 "\n"
			 "02                       ! batch item number\n"
			 "01                       ! integer label for this run of the program\n"
			 "\n"
			 "03                       ! batch item number\n"
			 "General test                            \n"
			 "\n"
			 "09                       ! batch item number\n"
			 "MAP.DAT                  ! name of map file\n"
			 "\n"
			 "10                       ! batch item number\n"
			 "LOCUS.DAT                ! name of locus file\n"
			 "\n"
			 "11                       ! batch item number\n"
			 "PEDIGREE.DAT             ! name of pedigree file\n"
			 "\n"
			 "12                       ! batch item number\n"
			 "F                        ! symbol for female (case insensitive)\n"
			 "M                        ! symbol for   male (case insensitive)\n"
			 "\n"
			 "13                       ! batch item number\n"
			 "N                        ! is trait listed in locus and pedigree files?\n"
			 "\n"
			 "18                       ! batch item number\n"
			 "0                        ! number of quantitative variables in pedigree file\n");
  fclose(BATCH2);

  INFO("Finished successfully!");

}


