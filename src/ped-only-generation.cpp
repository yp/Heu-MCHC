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
#include <vector>
#include <list>
#include <set>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <ctime>
#include <cstring>

#include "log.h"
#include "util.h"
#include "ped-only-generation-opt.h"
#include "log-build-info.h"

using namespace std;

size_t n_indiv;

double p_children;
double p_multimating;

char* file_prefix;
int seed;

class p_hap_t {
public:
  p_hap_t(int _id)
		:id(_id), fid(-1), mid(-1), gender(0)
  {}

  const static int MALE= 1;
  const static int FEMALE= 2;

  int id;
  int fid;
  int mid;
  int gender;   // 1 male 2 female
};


int main(int argc, char** argv) {
  INFO("PEDIGREE-ONLY GENERATION started.");
  PRINT_SYSTEM_INFORMATIONS;
  time_t rawtime;
  struct tm * timeinfo;
  int base_time= ((int)time (&rawtime));

  gengetopt_args_info args_info;

  if (cmdline_parser(argc, argv, &args_info) != 0) {
	 FATAL("Parameter error. "
			 "Try option --help for additional intormation about parameters.");
	 fail();
  }
// Check parameters
  if (args_info.n_individuals_arg<=0) {
	 FATAL("Parameter n-individuals has an invalid value (%d).",
			 args_info.n_individuals_arg);
	 fail();
  } else {
	 n_indiv= args_info.n_individuals_arg;
	 INFO("n-individuals:    " ST_FMTL(4) " ", n_indiv);
  }

  if ((args_info.children_probability_arg<0.0) ||
		(args_info.children_probability_arg>1.0)) {
	 FATAL("Parameter children-probability has an invalid value (%e).",
			 args_info.children_probability_arg);
	 fail();
  } else {
	 p_children= args_info.children_probability_arg;
	 INFO("children-probability: %+6f", p_children);
  }
  if ((args_info.multimating_probability_arg<0.0) ||
		(args_info.multimating_probability_arg>1.0)) {
	 FATAL("Parameter multimating-probability has an invalid value (%e).",
			 args_info.multimating_probability_arg);
	 fail();
  } else {
	 p_multimating= args_info.multimating_probability_arg;
	 INFO("multimating-probability: %+6f", p_multimating);
  }

  file_prefix= args_info.file_prefix_arg;
  INFO("file-prefix: %s", file_prefix);


  if (args_info.seed_arg<0) {
	 FATAL("Parameter seed has an invalid value (%d).",
			 args_info.seed_arg);
	 fail();
  } else {
	 if (args_info.seed_arg==0){
		seed= base_time;
	 } else {
		seed= args_info.seed_arg;
	 }
	 INFO("seed: %d", seed);
  }
  srand(seed);

  size_t n= 0; // no of individuals

  vector<p_hap_t> indivs;

  vector<int> indivs_as_founder;
  list<int> tmp_children;
// Topology creation
// first individual
  indivs.push_back(p_hap_t(0));
  n= 1;
  indivs_as_founder.push_back(0);

  while (n<n_indiv) {
// create a family given a child
		int f= indivs_as_founder[0];
		indivs_as_founder.erase(indivs_as_founder.begin());
		bool rip= false;
		do {
		  int m;
		  if (indivs_as_founder.size()>1 && rnd_01()<0.25) {
			 int pos= rnd_in_range(0, indivs_as_founder.size());
			 m= indivs_as_founder[pos];
			 INFO("Re-using pre-existing individual %4d.", m);
			 indivs_as_founder.erase(indivs_as_founder.begin()+pos);
		  } else {
			 m= n;
			 n= n+1;
			 indivs.push_back(p_hap_t(m));
		  }
		  if ((!rip) && (rnd_01()<0.5)) {
			 MY_SWAP(int, f, m);
		  }
		  double pc= 1.0;
		  while ((rnd_01()<pc) && ((n<n_indiv) || (pc==1.0))) {
			 int c= n;
			 n= n+1;
			 indivs.push_back(p_hap_t(c));
			 my_assert((indivs[f].gender == 0) ||
						  (indivs[f].gender == p_hap_t::MALE));
			 my_assert((indivs[m].gender == 0) ||
						  (indivs[m].gender == p_hap_t::FEMALE));
			 indivs[c].fid= f;
			 if (indivs[f].gender==0)
				indivs[f].gender= p_hap_t::MALE;
			 indivs[c].mid= m;
			 if (indivs[m].gender==0)
				indivs[m].gender= p_hap_t::FEMALE;
			 tmp_children.push_back(c);
			 pc= pc*p_children;
		  }
		  rip= rnd_01() < p_multimating;
		} while (rip && (n<n_indiv));
// Accoda i nuovi figli
		while (!tmp_children.empty()) {
		  indivs_as_founder.push_back(tmp_children.front());
		  tmp_children.pop_front();
		}
  }

  if (args_info.print_pedigree_given) {
	 FILE * pout;
	 if (strcmp(args_info.print_pedigree_file_arg, "-")==0) {
		pout= stdout;
	 } else {
		pout= fopen(args_info.print_pedigree_file_arg, "w");
	 }
	 if (pout==NULL) {
		ERROR("Opening file \"%s\" for writing failed!",
				args_info.print_pedigree_file_arg);
	 } else {
		set<int> mmnode;
		fprintf(pout, "digraph G {\n");
		for (size_t i= 0; i<indivs.size(); ++i) {
		  if (indivs[i].fid!=-1) {
			 if (mmnode.find(indivs[i].fid*n+indivs[i].mid)==mmnode.end()) {
				fprintf(pout, "  f%03dm%03d[label=\"\", shape=\"point\"];\n",
						  indivs[i].fid, indivs[i].mid);
				fprintf(pout, "  i%03d[shape=\"box\"];\n",
						  indivs[i].fid);
				fprintf(pout, "  i%03d[shape=\"ellipse\"];\n",
						  indivs[i].mid);
				fprintf(pout, "  f%03dm%03d[label=\"\", shape=\"point\"];\n",
						  indivs[i].fid, indivs[i].mid);
				fprintf(pout, "  i%03d -> f%03dm%03d;\n",
						  indivs[i].fid, indivs[i].fid, indivs[i].mid);
				fprintf(pout, "  i%03d -> f%03dm%03d;\n",
						  indivs[i].mid, indivs[i].fid, indivs[i].mid);
				mmnode.insert(indivs[i].fid*n+indivs[i].mid);
			 }
			 fprintf(pout, "  f%03dm%03d -> i%03d;\n",
					  indivs[i].fid, indivs[i].mid, indivs[i].id);
		  }
		}
		fprintf(pout, "}\n");
		INFO("Pedigree DOT description wrote on \"%s\".",
			  args_info.print_pedigree_file_arg);
		if (pout!=stdout)
		  fclose(pout);
	 }
  }

// Print the genotyped pedigree
  ostringstream sfilename;
  sfilename << file_prefix << ".txt";
  string filename= sfilename.str();
  DEBUG("Saving the pedigree on file \"%s\".", filename.c_str());
  ostringstream oss;
  oss << base_time;
  timeinfo = localtime ( &rawtime );
  ofstream fout(filename.c_str());
  fout << "# Pedigree generated " << asctime (timeinfo);
  fout << "# n-individuals:        " << n << endl;
  fout << "# children-probability: " << p_children << endl;
  fout << "# multimating-prob:     " << p_multimating << endl;
  fout << "# file-prefix:          " << file_prefix << endl;
  fout << "# seed:                 " << seed << endl;
  fout << "# format: id father mother gender(1=male,2=female)" << endl;
  fout << "input pedigree size " << n << endl;
  fout << "input pedigree record names 3 integers 1" << endl;
  fout << "input pedigree record father mother" << endl;
  fout << "input pedigree record gender present" << endl;
  fout << "input pedigree record observed absent" << endl;
  fout << "********" << endl;
  for (size_t i=0; i<n; ++i) {
	 if (indivs[i].gender==0)
		indivs[i].gender= p_hap_t::MALE;
	 fout << i+1 << " " << indivs[i].fid+1 << " "
			<< indivs[i].mid+1 << " " << indivs[i].gender << endl;
  }
  fout.close();
  INFO("PEDIGREE-ONLY GENERATION terminated.");
}



