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
#include "gt-ped-gen-opt.h"
#include "log-build-info.h"

using namespace std;

size_t n_indiv;
size_t len_hap;

double p_children;
double p_multimating;

double p_mutation;
double p_recombination;
size_t n_confs;

char* file_prefix;
int seed;

typedef vector<bool> hap_t;

typedef vector<int> gen_t;

class p_hap_t {
public:
  p_hap_t(int _id, size_t hap_len)
		:id(_id), fid(-1), mid(-1),
		 h1(hap_len, false), h2(hap_len, false),
		 mf(hap_len, false), mm(hap_len, false),
		 rf(hap_len, false), rm(hap_len, false)
  {}

  int id;
  int fid;
  int mid;
  hap_t h1;
  hap_t h2;
  hap_t mf;
  hap_t mm;
  hap_t rf;
  hap_t rm;
};

ostream&
operator<<(ostream& out, p_hap_t& ph) {
  for (size_t i= 0; i<len_hap; ++i) {
	 if (ph.h1[i]==ph.h2[i]) {
		out << (ph.h1[i]?1:0);
	 } else {
		out << 2;
	 }
  }
  return out;
}

ostream&
operator<<(ostream& out, hap_t& ph) {
  for (size_t i= 0; i<len_hap; ++i) {
	 if (ph[i]) {
		out << '1';
	 } else {
		out << '0';
	 }
  }
  return out;
}

static void
randomize_mutations(hap_t& m, hap_t& mut, bool infinite_site) {
  for (size_t i= 0; i<m.size(); ++i) {
	 if ((!infinite_site || !mut[i]) && rnd_01()<p_mutation) {
		m[i]= true;
		mut[i]= true;
	 } else {
		m[i]= false;
	 }
  }
}

static void
randomize_haplotype(hap_t& h) {
  for (size_t i= 0; i<h.size(); ++i) {
	 if (rnd_01()<0.5) {
		h[i]= true;
	 } else {
		h[i]= false;
	 }
  }
}

static void
randomize_recombinations(hap_t& r, const hap_t& h1, const hap_t& h2) {
  bool state= false;
  for (size_t i= 1; i<r.size(); ++i) {
	 if ((h1[i]!=h2[i]) && (rnd_01()<p_recombination)) {
		state= !state;
	 }
	 r[i]= state;
  }
}

static void
calculate_haplotype(hap_t& h,
						  const hap_t& m, const hap_t& r,
						  hap_t& h1p, hap_t& h2p) {
  hap_t* h1= &h1p;
  hap_t* h2= &h2p;
  if (rnd_01()<0.5) {
	 h1= &h2p;
	 h2= &h1p;
  }
  for (size_t i= 0; i<r.size(); ++i) {
	 bool state= (*h1)[i];
	 if (r[i])
		state= (*h2)[i];
	 if (m[i])
		state= !state;
	 h[i]= state;
  }
}

int main(int argc, char** argv) {
  INFO("GENERAL-PED-GENERATION started.");
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
  if (args_info.haplotype_length_arg<=0) {
	 FATAL("Parameter haplotype-length has an invalid value (%d).",
			 args_info.haplotype_length_arg);
	 fail();
  } else {
	 len_hap= args_info.haplotype_length_arg;
	 INFO("haplotype-length: " ST_FMTL(4) " ", len_hap);
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

  if ((args_info.mutation_probability_arg<0.0) ||
		(args_info.mutation_probability_arg>1.0)) {
	 FATAL("Parameter mutation-probability has an invalid value (%e).",
			 args_info.mutation_probability_arg);
	 fail();
  } else {
	 p_mutation= args_info.mutation_probability_arg;
	 INFO("mutation-probability: %+6f", p_mutation);
  }
  if ((args_info.recombination_probability_arg<0.0) ||
		(args_info.recombination_probability_arg>1.0)) {
	 FATAL("Parameter recombination-probability has an invalid value (%e).",
			 args_info.recombination_probability_arg);
	 fail();
  } else {
	 p_recombination= args_info.recombination_probability_arg;
	 INFO("recombination-probability: %+6f", p_recombination);
  }
  if (args_info.infinite_site_flag==0) {
	 INFO("Infinite-site assumption: NO");
  } else {
	 INFO("Infinite-site assumption: YES");
  }
  if (args_info.n_configurations_arg<=0) {
	 FATAL("Parameter n-configurations has an invalid value (%d).",
			 args_info.n_configurations_arg);
	 fail();
  } else {
	 n_confs= args_info.n_configurations_arg;
	 INFO("n-configuration: " ST_FMT, n_confs);
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
  indivs.push_back(p_hap_t(n, len_hap));
  indivs.push_back(p_hap_t(n, len_hap));
  indivs.push_back(p_hap_t(n, len_hap));
  n= 3;
  indivs_as_founder.push_back(0);
  indivs_as_founder.push_back(1);
  indivs_as_founder.push_back(2);

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
			 indivs.push_back(p_hap_t(m, len_hap));
		  }
		  if ((!rip) && (rnd_01()<0.5)) {
			 MY_SWAP(int, f, m);
		  }
		  double pc= 1.0;
		  while ((rnd_01()<pc) && ((n<n_indiv) || (pc==1.0))) {
			 int c= n;
			 n= n+1;
			 indivs.push_back(p_hap_t(c, len_hap));
			 indivs[c].fid= f;
			 indivs[c].mid= m;
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

  for (size_t j= 0; j<n_confs; ++j) {
// Calculate a random haplotype configuration
	 INFO("Calculate random haplotype configuration " ST_FMTL(3) " /" ST_FMTL(3) " .", j, n_confs);
	 vector<bool> iniz(n, false);
	 hap_t tmut(len_hap, false);
	 list<int> queue;
	 for (size_t i= 0; i<n; ++i) {
		if (!iniz[i]) {
		  queue.push_front(i);
		  while (!queue.empty()) {
			 int p= queue.front();
			 queue.pop_front();
			 if (iniz[p])
				continue;
			 if (indivs[p].fid!=-1) {
// not a founder
				int f= indivs[p].fid;
				int m= indivs[p].mid;
				if (!iniz[f] || !iniz[m]) {
// parents not initialized
				  queue.push_front(p);
				  if (!iniz[m])
					 queue.push_front(m);
				  if (!iniz[f])
					 queue.push_front(f);
				} else {
// parents initialized
// initialize events
				  randomize_mutations(indivs[p].mf, tmut, args_info.infinite_site_flag!=0);
				  randomize_mutations(indivs[p].mm, tmut, args_info.infinite_site_flag!=0);
				  randomize_recombinations(indivs[p].rf, indivs[f].h1, indivs[f].h2);
				  randomize_recombinations(indivs[p].rm, indivs[m].h1, indivs[m].h2);
// calculate haplotypes
				  calculate_haplotype(indivs[p].h1,
											 indivs[p].mf, indivs[p].rf,
											 indivs[f].h1, indivs[f].h2);
				  calculate_haplotype(indivs[p].h2,
											 indivs[p].mm, indivs[p].rm,
											 indivs[m].h1, indivs[m].h2);
				  iniz[p]= true;
				}
			 } else {
// founder
				randomize_haplotype(indivs[p].h1);
				randomize_haplotype(indivs[p].h2);
				iniz[p]= true;
			 }
		  }
		}
	 }
// Print the genotyped pedigree
	 ostringstream sfilename;
	 sfilename << file_prefix << "conf-" << j << ".txt";
	 string filename= sfilename.str();
	 DEBUG("Saving the genotyped pedigree on file \"%s\".", filename.c_str());
	 ostringstream oss;
	 oss << base_time;
	 timeinfo = localtime ( &rawtime );
	 ofstream fout(filename.c_str());
	 fout << "# General tree pedigree generated " << asctime (timeinfo);
	 fout << "# configuration number: " << j << endl;
	 fout << "# mutation-probability: " << p_mutation << endl;
	 fout << "# recomb-probability:   " << p_recombination << endl;
	 fout << "# n-individuals:        " << n << endl;
	 fout << "# haplotype-length:     " << len_hap << endl;
	 fout << "# children-probability: " << p_children << endl;
	 fout << "# multimating-prob:     " << p_multimating << endl;
	 fout << "# n-configuration:      " << n_confs << endl;
	 fout << "# file-prefix:          " << file_prefix << endl;
	 fout << "# seed:                 " << seed << endl;
	 fout << len_hap << endl;
	 fout << n << endl;
	 fout << "# genotypes" << endl;
	 for (size_t i=0; i<n; ++i) {
		fout << "#| " << i << " G " << indivs[i] << endl;
		fout << "#| " << i << " P " << indivs[i].h1 << endl;
		fout << "#| " << i << " M " << indivs[i].h2 << endl;
		fout << indivs[i] << endl;
	 }
	 fout << "# events " << endl;
	 for (size_t i=0; i<n; ++i) {
		for (size_t j= 0; j<len_hap; ++j) {
		  if (indivs[i].mf[j]) {
			 fout << "#* MUT "<< i << " " << j << " 0" << endl;
		  }
		  if (indivs[i].mm[j]) {
			 fout << "#* MUT "<< i << " " << j << " 1" << endl;
		  }
		}
	 }
	 for (size_t i=0; i<n; ++i) {
		for (size_t j= 1; j<len_hap; ++j) {
		  if (indivs[i].rf[j] != indivs[i].rf[j-1]) {
			 fout << "#* REC "<< i << " " << j << " 0" << endl;
		  }
		  if (indivs[i].rm[j] != indivs[i].rm[j-1]) {
			 fout << "#* REC "<< i << " " << j << " 1" << endl;
		  }
		}
	 }
	 fout << "# pedigree Father Mother Child" << endl;
	 for (size_t i=0; i<n; ++i) {
		if (indivs[i].fid!=-1) {
		  fout << indivs[i].fid << " " << indivs[i].mid << " "
				 << indivs[i].id << endl;
		}
	 }
	 fout << "# end" << endl;
	 fout.close();
  }
  INFO("GENERAL-PED-GENERATION terminated.");
}



