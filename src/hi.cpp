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
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "hi-options.h"
#include "log.h"
#include "util.h"
#include "gen-ped-IO.h"
#include "data.hpp"
#include "m4ri/m4ri.h"
#include "log-build-info.h"

#define LEN_BUFFER 10001

static void
read_events_from_file(FILE* fev, ped::e_vars_t& evs,
							 const size_t n_indiv, const size_t n_loci) {
  my_assert(fev!=NULL);
  size_t row= 0;
  int len= 0;
  size_t size= LEN_BUFFER;
  char * BUFF= c_palloc(size);
  char * skind= c_palloc(LEN_BUFFER);
  while (!feof(fev)) {
	 len= my_getline(&BUFF, &size, fev);
	 ++row;
	 if (len>0 && BUFF[0]=='=') {
		int read= sscanf(BUFF+1, "%s", skind);
// Look for the name
		bool found= false;
		unsigned int kind= 0;
		unsigned int i= 0;
		while (read>0 && i<ped::e_variable_t::N_KINDS && !found) {
		  if (strcmp(skind, ped::kind_names[i])==0) {
			 kind= i;
			 found= true;
		  }
		  ++i;
		}
		if (found) {
		  DEBUG("New event read \"%s\".", BUFF);
		  unsigned int indiv, locus, parent;
		  read= sscanf(BUFF+1, "%s %u %u %u", skind, &indiv, &locus, &parent);
		  if (read<4) {
			 ERROR("Malformed row " ST_FMT " \"%s\". Row ignored.", row, BUFF);
		  } else {
			 if (indiv>=n_indiv) {
				ERROR("Event \"%s\" refers to an invalid individual. Ignored.", BUFF);
			 } else if (locus>=n_loci) {
				ERROR("Event \"%s\" refers to an invalid locus. Ignored.", BUFF);
			 } else if (parent!=0 && parent!=1) {
				ERROR("Event \"%s\" refers to an invalid parent. Ignored.", BUFF);
			 } else {
				ped::e_variable_t e(kind, indiv, locus, parent);
				evs.insert(e);
			 }
		  }
		} else {
		  DEBUG("Row ignored \"%s\".", BUFF);
		}
	 }
  }
  INFO("Read " ST_FMT " events.", evs.size());

  pfree(BUFF);
  pfree(skind);
}

#undef LEN_BUFFER


static void
integrate_events_in_pedigree(pgenped gp, ped::e_vars_t evs) {
// Initialize structs.
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 for (unsigned int l= 0; l<gp->n_loci; ++l) {
		gp->individuals[i]->df[l]= 0;
		gp->individuals[i]->dm[l]= 0;
		gp->individuals[i]->rm[l]= 0;
		gp->individuals[i]->rf[l]= 0;
	 }
  }
// Fill structs
  for (ped::e_vars_t::const_iterator it= evs.begin();
		 it!=evs.end();
		 ++it) {
	 const ped::e_variable_t& e= *it;
	 if (e.kind()==ped::e_variable_t::MUT) {
		if (e.p==0) {
		  gp->individuals[e.i]->df[e.l]= 1;
		} else {
		  gp->individuals[e.i]->dm[e.l]= 1;
		}
	 } else if (e.kind()==ped::e_variable_t::REC) {
		int * r= NULL;
		if (e.p==0) {
		  r= gp->individuals[e.i]->rf;
		} else {
		  r= gp->individuals[e.i]->rm;
		}
		for (unsigned int l= e.l; l<gp->n_loci; ++l) {
		  r[l]= (r[l]+1)%2;
		}
	 }
  }
}

static int
get_from_map(map<int*, int> varsmap, int* pt) {
  my_assert(varsmap.find(pt)!=varsmap.end());
  return varsmap[pt];
}

static void
solve_pedigree(pgenped gp) {
  my_assert(gp!=NULL);
  pindiv* const is= gp->individuals;
// prepare known phases
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 for (unsigned int l= 0; l<gp->n_loci; ++l) {
		int phase= -1;
		int np;
		if (is[i]->g[l]!=GEN2) {
// i is homozygous
		  phase= is[i]->g[l];
		} else if (is[i]->f!=NULL && is[i]->f->g[l]!=GEN2) {
// father homozygous
		  DEBUG("Father of %3d homozygous at locus %3d with phase %d.", i, l, is[i]->f->g[l]);
		  np= (is[i]->f->g[l] + is[i]->df[l]) % 2;
		  my_assert(phase == -1 || phase == np);
		  phase= np;
		} else if (is[i]->m!=NULL && is[i]->m->g[l]!=GEN2) {
// mother homozygous
		  DEBUG("Mother of %3d homozygous at locus %3d with phase %d.", i, l, is[i]->m->g[l]);
		  np= (is[i]->m->g[l] + is[i]->dm[l] + 1) % 2;
		  my_assert(phase == -1 || phase == np);
		  phase= np;
		}
		is[i]->p[l]= phase;
		if (phase!=-1) {
		  DEBUG("Individual %3d predetermined at locus %3d at phase %d.", i, l, phase);
		}
	 }
  }

// create the variable set and count the equations
  size_t tmpnv= 0;
  size_t tmpne= 0;
  map<int*, int> varsmap;
  vector<int*> vars;
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 TRACE("Trio F %3d  M %3d -> C %3d", is[i]->fi, is[i]->mi, i);
	 for (unsigned int l= 0; l<gp->n_loci; ++l) {
// i is heterozygous
		if (is[i]->p[l]==-1) {
		  varsmap[&is[i]->p[l]]= tmpnv;
		  vars.push_back(&(is[i]->p[l]));
		  DEBUG("Variable " ST_FMTL(5) "  p_%u_%u", tmpnv, i, l);
		  ++tmpnv;
		}
	 }
	 if (is[i]->fi!=-1) {
// i is not a founder
		varsmap[&is[i]->h_f]= tmpnv;
		vars.push_back(&(is[i]->h_f));
		DEBUG("Variable " ST_FMTL(5) "  hf_%u.", tmpnv, i);
		++tmpnv;
		varsmap[&is[i]->h_m]= tmpnv;
		vars.push_back(&(is[i]->h_m));
		DEBUG("Variable " ST_FMTL(5) "  hm_%u.", tmpnv, i);
		++tmpnv;
		for (unsigned int l= 0; l<gp->n_loci; ++l) {
		  if (is[i]->f->g[l]==GEN2)
			 ++tmpne;
		  if (is[i]->m->g[l]==GEN2)
			 ++tmpne;
		}
	 }
  }
  const size_t nvar= tmpnv;
  const size_t neq= tmpne;
  INFO("The linear system has " ST_FMTL(5) "  variables and " ST_FMTL(5) "  equations.", nvar, neq);
// create the binary matrix A and the vector B
  mzd_t *A= mzd_init(neq, nvar);
  mzd_t *B= mzd_init(neq, 1);
  for (size_t i= 0; i<neq; ++i) {
	 mzd_row_clear_offset(A, i, 0);
	 mzd_write_bit(B, i, 0, 0);
  }
  tmpne= 0;
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 if (is[i]->fi!=-1) {
// i is not a founder
		const int hf= get_from_map(varsmap, &is[i]->h_f);
		const int hm= get_from_map(varsmap, &is[i]->h_m);
		for (unsigned int l= 0; l<gp->n_loci; ++l) {
		  DEBUG("Individual %4u,  locus %4u", i, l);
		  if (is[i]->f->g[l]==GEN2) {
			 DEBUG_NOTN("eq " ST_FMTL(5) " ) ", tmpne);
			 int kt= 0;
			 if (is[i]->f->p[l]!=-1) {
				kt+= is[i]->f->p[l];
			 } else {
				mzd_write_bit(A, tmpne, get_from_map(varsmap, &is[i]->f->p[l]), 1);
				DEBUG_MSG("p_%d_%d + ", is[i]->fi, l);
			 }
			 mzd_write_bit(A, tmpne, hf, 1);
			 DEBUG_MSG("hf_%d + ", i);
			 if (is[i]->p[l]==-1) {
				mzd_write_bit(A, tmpne, get_from_map(varsmap, &is[i]->p[l]), 1);
				DEBUG_MSG("p_%d_%d + ", i, l);
			 } else {
				kt+= is[i]->p[l];
			 }
			 kt= (kt + is[i]->df[l] + is[i]->rf[l]) % 2;
			 mzd_write_bit(B, tmpne, 0, kt);
			 DEBUG_MSG("%d = 0\n", kt);
			 ++tmpne;
		  }
		  if (is[i]->m->g[l]==GEN2) {
			 DEBUG_NOTN("eq " ST_FMTL(5) " ) ", tmpne);
			 int kt= (is[i]->g[l]==GEN2?1:0);
			 if (is[i]->m->p[l]!=-1) {
				kt= kt + is[i]->m->p[l];
			 } else {
				mzd_write_bit(A, tmpne, get_from_map(varsmap, &is[i]->m->p[l]), 1);
				DEBUG_MSG("p_%d_%d + ", is[i]->mi, l);
			 }
			 mzd_write_bit(A, tmpne, hm, 1);
			 DEBUG_MSG("hm_%d + ", i);
			 if (is[i]->p[l]==-1) {
				mzd_write_bit(A, tmpne, get_from_map(varsmap, &is[i]->p[l]), 1);
				DEBUG_MSG("p_%d_%d + ", i, l);
			 } else {
				kt= kt + is[i]->p[l];
			 }
			 kt= (kt + is[i]->dm[l] + is[i]->rm[l]) % 2;
			 mzd_write_bit(B, tmpne, 0, kt);
			 DEBUG_MSG("%d = 0\n", kt);
			 ++tmpne;
		  }
		}
	 }
  }
  mzd_t* Aext= mzd_init(neq, nvar+1);
  for (size_t i= 0; i<neq; ++i) {
	 for (size_t j= 0; j<nvar; ++j)
		mzd_write_bit(Aext, i, j, mzd_read_bit(A, i, j));
	 mzd_write_bit(Aext, i, nvar, mzd_read_bit(B, i, 0));
  }
  const size_t rank= mzd_echelonize_m4ri(A, 0, 0);
  const size_t rank2= mzd_echelonize_m4ri(Aext, 0, 0);

  INFO("Rank " ST_FMT "", rank);
  if (rank != rank2) {
	 FATAL("Rank of coefficient matrix and complete matrix differ.");
	 fail();
  }
  for (size_t i= 0; i<neq; ++i) {
	 mzd_write_bit(B, i, 0, mzd_read_bit(Aext, i, nvar));
  }
  mzd_free(Aext);
  for (int i= rank-1; i>=0; --i) {
	 int first_one= i;
	 while (mzd_read_bit(A, i, first_one)==0)
		++first_one;
	 int kt= mzd_read_bit(B, i, 0);
	 for (int j= 0; j<=i; ++j) {
		if (mzd_read_bit(A, j, first_one)==1) {
		  mzd_write_bit(B, j, 0, (mzd_read_bit(B, j, 0)+kt)%2);
		  mzd_write_bit(A, j, first_one, 0);
		}
	 }
	 *(vars[first_one])= kt;
  }
  for (size_t i= 0; i<neq; ++i) {
	 for (size_t j= i; j<nvar; ++j) {
		if (mzd_read_bit(A, i, j)==1) {
		  DEBUG("Redundant variable " ST_FMTL(5) "  arbitrarly set to 0.", j);
		  *(vars[j])= 0;
		}
	 }
  }
  mzd_free(A);
  mzd_free(B);
}

static void
compute_haplotypes(pgenped gp) {
  pindiv* const is= gp->individuals;
// prepare known phases
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 for (unsigned int l= 0; l<gp->n_loci; ++l) {
		is[i]->hv_f[l]= is[i]->p[l];
		is[i]->hv_m[l]= (is[i]->p[l]+(is[i]->g[l]==GEN2?1:0))%2;
		bool c1= (is[i]->g[l]!=GEN2 || is[i]->hv_f[l]!=is[i]->hv_m[l]);
		bool c2= (is[i]->g[l]!=GEN0 ||
					 ((is[i]->hv_f[l]==is[i]->hv_m[l]) &&
					  (is[i]->hv_f[l]==0)));
		bool c3= (is[i]->g[l]!=GEN1 ||
					 ((is[i]->hv_f[l]==is[i]->hv_m[l]) &&
					  (is[i]->hv_f[l]==1)));
		if (!(c1 && c2 && c3)) {
		  FATAL("Genotype-haplotype consistency violated in individual "
				  "%5u at locus %5u.", i, l);
		  fail();
		}
	 }
  }
}

static void
print_haplotypes(pgenped gp) {
  pindiv* const is= gp->individuals;
// prepare known phases
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 if (is[i]->fi == -1) {
		printf("%5d FOUNDER\n", i);
	 }
	 printf("%5u G ", i);
	 for (unsigned int l= 0; l<gp->n_loci; ++l) {
		printf("%d", is[i]->g[l]);
	 }
	 printf("\n");
	 printf("%5u P ", i);
	 for (unsigned int l= 0; l<gp->n_loci; ++l) {
		printf("%d", is[i]->hv_f[l]);
	 }
	 printf("\n");
	 printf("%5u M ", i);
	 for (unsigned int l= 0; l<gp->n_loci; ++l) {
		printf("%d", is[i]->hv_m[l]);
	 }
	 printf("\n");
  }
}


static bool
is_present(const int* h, const int* p1, const int* p2, const size_t n_loci,
			  const int* m, const int* r) {
  bool f1= true;
  const int* t[2]= {p1, p2};
  for (size_t l= 0; l<n_loci; ++l) {
	 f1= f1 && h[l]==((t[r[l]][l]+m[l])%2);
  }
  bool f2= true;
  for (size_t l= 0; l<n_loci; ++l) {
	 f2= f2 && h[l]==((t[1-r[l]][l]+m[l])%2);
  }
  return f1||f2;
}

static void
check_inheritance(pgenped gp) {
  pindiv* const is= gp->individuals;
  for (unsigned int i= 0; i<gp->n_indiv; ++i) {
	 if (is[i]->fi!=-1) {
		DEBUG("Analysis individual %5u.", i);
		bool pat= is_present(is[i]->hv_f, is[i]->f->hv_f, is[i]->f->hv_m, gp->n_loci,
									is[i]->df, is[i]->rf);
		bool mat= is_present(is[i]->hv_m, is[i]->m->hv_f, is[i]->m->hv_m, gp->n_loci,
									is[i]->dm, is[i]->rm);
		if (!pat) {
		  FATAL("Paternal haplotype of individual %5d is not inherited.", i);
		  fail();
		}
		if (!mat) {
		  FATAL("Maternal haplotype of individual %5d is not inherited.", i);
		  fail();
		}
	 }
  }
}


int main(int argc, char ** argv) {
  INFO("TO-HAP started.");
  PRINT_SYSTEM_INFORMATIONS;
  gengetopt_args_info args_info;

  if (cmdline_parser(argc, argv, &args_info) != 0) {
	 FATAL("Parameter error. "
			 "Try option --help for additional intormation about parameters.");
	 fail();
  }

  INFO("Pedigree file: \"%s\"", args_info.pedigree_arg);
  INFO("Events file:   \"%s\"", args_info.events_arg);
  FILE * fped;
  bool stin= false;
  if (strcmp(args_info.pedigree_arg, "-")==0) {
	 stin= true;
	 fped= stdin;
  } else {
	 fped= fopen(args_info.pedigree_arg, "r");
	 if (fped==NULL) {
		FATAL("Opening file \"%s\" failed.", args_info.pedigree_arg);
		fail();
	 }
  }

  FILE * fev;
  if (strcmp(args_info.events_arg, "-")==0) {
	 if (stin) {
		FATAL("At most one of the files can be \"-\".");
		fail();
	 }
	 fev= stdin;
  } else {
	 fev= fopen(args_info.events_arg, "r");
	 if (fev==NULL) {
		FATAL("Opening file \"%s\" failed.", args_info.events_arg);
		if (!stin)
		  fclose(fped);
		fail();
	 }
  }

  pgenped gp= gp_read_from_file(fped);

  ped::e_vars_t evs;
  read_events_from_file(fev, evs, gp->n_indiv, gp->n_loci);

  if (!stin)
	 fclose(fped);
  if (fev!=stdin)
	 fclose(fev);

  integrate_events_in_pedigree(gp, evs);
  solve_pedigree(gp);
  compute_haplotypes(gp);
  print_haplotypes(gp);
  check_inheritance(gp);

  gp_destroy(gp);
  cmdline_parser_free(&args_info);

  INFO("TO-HAP terminated.");
  return 0;
}
