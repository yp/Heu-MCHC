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
#include "locus-graph.hpp"

#include "log.h"
#include "my_time.h"
#include "util.h"
#include "gen-ped-IO.h"
#include "solve.hpp"
#include "log-build-info.h"

#define DEFAULT_GAMMA 0.20

#ifndef MAX_TENT
#define MAX_TENT 10
#endif

#ifndef MAX_EVENTS
#define MAX_EVENTS 2000
#endif

using namespace std;

pmytime pt_gauss;
#if (defined NO_MUTATIONS) || (defined NO_RECOMBINATIONS)
// if mutations or recombinations are not admitted, do no prefer
// anything
static const BP_T p1min[2]= {0.001, 0.001};
static const BP_T p1max[2]= {0.125, 0.125};
#else
// if mutations and recombinations are admitted, prefer recombinations
static const BP_T p1min[2]= {0.00075, 0.001};
static const BP_T p1max[2]= {0.10000, 0.125};
#endif

int main(int argc, char * argv[]) {
  double gamma= DEFAULT_GAMMA;
  char* foutname= "out.txt";
  if (argc>1) {
	 foutname= argv[1];
  }
  if (argc>2) {
	 gamma= atoi(argv[2])/100.0;
  }
  INFO("PED-HI started.");
  PRINT_SYSTEM_INFORMATIONS;
  pmytime pt_gen= MYTIME_create_with_name("general timer");
  pt_gauss= MYTIME_create_with_name("gauss timer");
  MYTIME_start(pt_gen);
  pgenped gp= gp_read_from_file(stdin);

  BP_T* p1= NPALLOC(BP_T, ped::e_variable_t::N_KINDS);
  ped::e_vars_t* mut_vars= NULL;
  size_t max_mut= MAX_EVENTS;
  for (size_t tent= 0; tent<MAX_TENT; ++tent) {
	 srand(tent*1981);
	 for (size_t i= 0; i<ped::e_variable_t::N_KINDS; ++i) {
#if MAX_TENT > 1
		  p1[i]= p1min[i]+((p1max[i]-p1min[i])*(MAX_TENT-tent-1)/(MAX_TENT-1));
#else
		  p1[i]= p1max[i];
#endif
	 }
	 try {
		INFO("**** Starting new trial.");
		ped::e_vars_t* ris= calculate_minimum_solution(gp, max_mut, gamma, p1);
		INFO("Trial " ST_FMTL(3) " . No. of events= " ST_FMT "", tent+1, ris->size());
		STATS("# trial " ST_FMT "", tent+1);
		STATS("# predicted events " ST_FMT "", ris->size());
		if (mut_vars==NULL || ris->size()<mut_vars->size()) {
		  if (mut_vars!=NULL)
			 delete mut_vars;
		  mut_vars= ris;
		  max_mut= mut_vars->size();
		  FILE* fout= fopen(foutname, "w");
		  fprintf(fout, "= RESULTING EVENTS\n");
		  for(ped::e_vars_t::const_iterator it= mut_vars->begin();
				it!=mut_vars->end();
				++it) {
			 fprintf(fout, "= %s %4d %4d %2d\n", ped::kind_names[(*it).kind()], (*it).i, (*it).l, (*it).p);
		  }
		  fclose(fout);
		} else {
		  delete ris;
		}
	 } catch (ped_hi_exception& e) {
		ERROR("Trial " ST_FMTL(3) "  failed. \"%s\"", tent+1, e.get_message().c_str());
	 }
  }
  if (mut_vars!=NULL) {
	 FILE* fout= fopen(foutname, "w");
	 fprintf(fout, "= RESULTING EVENTS\n");
	 for(ped::e_vars_t::const_iterator it= mut_vars->begin();
		  it!=mut_vars->end();
		  ++it) {
		fprintf(fout, "= %s %4d %4d %2d\n", ped::kind_names[(*it).kind()], (*it).i, (*it).l, (*it).p);
	 }
	 fclose(fout);
	 delete mut_vars;
  } else {
	 FATAL("No trials ended with a valid mutation set.");
  }

  pfree(p1);
  gp_destroy(gp);
  MYTIME_stop(pt_gen);
  MYTIME_LOG(pt_gauss);
  MYTIME_LOG(pt_gen);
  INFO("PED-HI stopped.");
  MYTIME_destroy(pt_gen);
  MYTIME_destroy(pt_gauss);
  resource_usage_log();
  return 0;
}

