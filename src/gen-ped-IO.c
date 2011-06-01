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
#include "gen-ped-IO.h"
#include "util.h"
#include "log.h"
#include <string.h>
#include <stdlib.h>

/**
 * Read a file with the following structure:
 * n_loci
 * n_indiv
 * gen_1
 * gen_2
 * ...
 * gen_n
 * id_father id_mother id_children
 * ...
 * id_father id_mother id_children
 *
 * It ignores any row starting with a "#".
 **/


#define LEN_BUFFER 10001

pindiv indiv_create(const int n_loci, const int n_individuals) {
  my_assert(n_loci>0);
  my_assert(n_individuals>0);
  pindiv pi= PALLOC(struct _indiv);
  pi->id= -1;
  pi->fi= pi->mi= -1;
  pi->f= pi->m= NULL;
  pi->g= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->g, 0, n_loci);
  pi->p= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->p, 0, n_loci);
  pi->df= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->df, 0, n_loci);
  pi->dm= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->dm, 0, n_loci);
  pi->rf= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->rf, 0, n_loci);
  pi->rm= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->rm, 0, n_loci);
  pi->h_f= pi->h_m= -1;

  pi->hv_f= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->hv_f, -1, n_loci);
  pi->hv_m= NPALLOC(int, n_loci);
  _SET_VECTOR(pi->hv_m, -1, n_loci);

  pi->n_children= 0;
  pi->children= NPALLOC(pindiv, n_individuals);
  _SET_VECTOR(pi->children, NULL, n_individuals);

  return pi;
}


void indiv_destroy(pindiv pi) {
  my_assert(pi!=NULL);
  pfree(pi->g);
  pfree(pi->p);
  pfree(pi->df);
  pfree(pi->dm);
  pfree(pi->rf);
  pfree(pi->rm);
  pfree(pi->hv_f);
  pfree(pi->hv_m);
  pfree(pi->children);
  pfree(pi);
}


pgenped gp_create(const int n_loci, const int n_individuals) {
  my_assert(n_loci>0);
  my_assert(n_individuals>0);
  pgenped gp= PALLOC(struct _genped);
  gp->n_loci= n_loci;
  gp->n_indiv= n_individuals;
  gp->individuals= NPALLOC(pindiv, n_individuals);
  for (int i= 0; i<n_individuals; ++i) {
	 gp->individuals[i]=indiv_create(n_loci, n_individuals);
	 gp->individuals[i]->id= i;
  }
  return gp;
}

static pindiv
indiv_copy(pindiv src, const unsigned int n_loci, const unsigned int n_individuals) {
  pindiv ris= indiv_create(n_loci, n_individuals);
  ris->id= src->id;
  ris->f= src->f;
  ris->m= src->m;
  ris->fi= src->fi;
  ris->mi= src->mi;
  ris->h_f= src->h_f;
  ris->h_m= src->h_m;
  ris->n_children= src->n_children;
  for (unsigned int i= 0; i<n_loci; ++i) {
	 ris->g[i]= src->g[i];
	 ris->p[i]= src->p[i];
	 ris->df[i]= src->df[i];
	 ris->dm[i]= src->dm[i];
	 ris->rf[i]= src->rf[i];
	 ris->rm[i]= src->rm[i];
	 ris->hv_f[i]= src->hv_f[i];
	 ris->hv_m[i]= src->hv_m[i];
  }
  for (unsigned int i= 0; i<ris->n_children; ++i) {
	 ris->children[i]= src->children[i];
  }
  return ris;
}

pgenped gp_copy(pgenped src) {
  pgenped ris= PALLOC(struct _genped);
  ris->n_loci= src->n_loci;
  ris->n_indiv= src->n_indiv;
  ris->individuals= NPALLOC(pindiv, src->n_indiv);
  for (unsigned int i= 0; i<src->n_indiv; ++i) {
	 ris->individuals[i]=
		indiv_copy(src->individuals[i],
					  src->n_loci,
					  src->n_indiv);
  }
  for (unsigned int i= 0; i<src->n_indiv; ++i) {
	 pindiv const c= ris->individuals[i];
	 if (c->fi>=0)
		c->f= ris->individuals[c->fi];
	 if (c->mi>=0)
		c->m= ris->individuals[c->mi];
	 for (unsigned int j= 0; j<c->n_children; ++j) {
		c->children[j]= ris->individuals[c->children[j]->id];
	 }
  }
  return ris;
}

void gp_destroy(pgenped gp) {
  my_assert(gp!=NULL);
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 indiv_destroy(gp->individuals[i]);
  }
  pfree(gp->individuals);
  pfree(gp);
}

void gp_add_trio(pgenped gp,
					  const unsigned int fi,
					  const unsigned int mi,
					  const unsigned int ci) {
  my_assert(gp!=NULL);
  assert_ulim(fi, gp->n_indiv);
  assert_ulim(mi, gp->n_indiv);
  assert_ulim(ci, gp->n_indiv);
  pindiv f= gp->individuals[fi];
  pindiv m= gp->individuals[mi];
  pindiv c= gp->individuals[ci];
  c->fi= fi;
  c->f= f;
  c->mi= mi;
  c->m= m;
  f->children[f->n_children]= c;
  f->n_children++;
  m->children[m->n_children]= c;
  m->n_children++;
}

static int
char2gen(const char g) {
  switch (g) {
	 case '0': return 0;
	 case '1': return 1;
	 case '2': return 2;
	 default:
		ERROR("Read an invalid genotype symbol >%c<. Terminating.", g);
		fail();
  }
}


pgenped gp_read_from_file(FILE* fin) {
  my_assert(fin!=NULL);

  int row= 0;
  int n_loci;
  int n_indiv;
  size_t size= LEN_BUFFER;
  char * BUFF= c_palloc(size);
  int len;

  bool read_n_loci= false;
  bool read_n_indiv= false;
  while (!read_n_loci) {
	 len= my_getline(&BUFF, &size, fin);
	 ++row;
	 if (len>0 && BUFF[0]!='#') {
		read_n_loci= sscanf(BUFF, "%d", &n_loci) > 0 && n_loci>0;
	 }
  }
  while (!read_n_indiv) {
	 len= my_getline(&BUFF, &size, fin);
	 ++row;
	 if (len>0 && BUFF[0]!='#') {
		read_n_indiv= sscanf(BUFF, "%d", &n_indiv) > 0 && n_indiv>0;
	 }
  }

  DEBUG("The genotyped pedigree has %d loci and %d individuals.", n_loci, n_indiv);

  pgenped gp= gp_create(n_loci, n_indiv);

// Reading the genotypes
  int i= 0;
  while (i<n_indiv) {
	 len= my_getline(&BUFF, &size, fin);
	 ++row;
	 if (len>0 && BUFF[0]!='#') {
		if (len!=n_loci) {
		  ERROR("Read a genotype at line %d with a number of characters different to the number of loci. Read >%s<.", row, BUFF);
		  ERROR("Terminating.");
		  fail();
		} else {
		  for (int j= 0; j<n_loci; ++j) {
			 gp->individuals[i]->g[j]= char2gen(BUFF[j]);
		  }
		  ++i;
		}
	 }
  }

  while (!feof(fin)) {
	 len= my_getline(&BUFF, &size, fin);
	 ++row;
	 if (len>0 && BUFF[0]!='#') {
		int fi, mi, ci;
		int el= sscanf(BUFF, "%d %d %d", &fi, &mi, &ci);
		if (el<3) {
		  ERROR("The trio read at line %d is invalid. Read >%s<. Terminating.", row, BUFF);
		  fail();
		} else {
		  gp_add_trio(gp, fi, mi, ci);
		  TRACE("Added trio (F, M, C)= (%4d, %4d, %4d).", fi, mi, ci);
		}
	 }
  }

  pfree(BUFF);

  return gp;
}





#undef LEN_BUFFER

