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
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <ctime>

#include "util.h"

#include "log.h"
#include "log-build-info.h"

using namespace std;

const double p_mutation= 0.008;
const double p_recombination= 0.002;
const size_t len_hap= 60;
const double probability_factor= 0.94;
const double MAX_LEVEL= 5.0;
int SEED;

typedef vector<bool> hap_t;

typedef vector<int> gen_t;

class p_hap_t {
public:
  ~p_hap_t() {
	 delete h1;
	 delete h2;
  }
  hap_t* h1;
  hap_t* h2;
  hap_t* mutations;
  hap_t* recombinations;
  bool swapped;
};

ostream&
operator<<(ostream& out, p_hap_t& ph) {
  for (size_t i= 0; i<len_hap; ++i) {
	 if ((*ph.h1)[i]==(*ph.h2)[i]) {
		out << ((*ph.h1)[i]?1:0);
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


class trio {
public:
  int f;
  int m;
  int c;

  trio(int f_, int m_, int c_)
		:f(f_), m(m_), c(c_)
  {}
};

ostream&
operator<<(ostream& out, const trio& t) {
  out << t.f << ' ' << t.m << ' ' << t.c;
  return out;
}

static bool
rnd_bool(void) {
  static int DIV= (RAND_MAX>>1);
  return (rand()>=DIV);
}

static string&
from_hap_to_str(const hap_t& h) {
  static string s;
  s= string(len_hap, '0');
  for (size_t i= 0; i<len_hap; ++i) {
	 if(h[i])
		s[i]= '1';
	 else
		s[i]= '0';
  }
  return s;
}

static hap_t*
generate_rnd_hap(void) {
  hap_t* ph= new hap_t(len_hap, false);
  hap_t& h= *ph;
  for (size_t i= 0; i<len_hap; ++i) {
	 h[i]= rnd_bool();
  }
  DEBUG("Generated random haplotype %s", from_hap_to_str(h).c_str());
  return ph;
}

static p_hap_t*
generate_rnd_p_hap(void) {
  p_hap_t* p= new p_hap_t();
  p->h1= generate_rnd_hap();
  p->h2= generate_rnd_hap();
  p->mutations= new hap_t(len_hap, false);
  p->recombinations= new hap_t(len_hap, false);
  return p;
}

static void
mutate_haplotype1(size_t id, size_t role, p_hap_t* p,
						ostringstream& osm, size_t& n_mutations,
						double p_mutation) {
// Pointwise mutations
  for (size_t i= 0; i<len_hap; ++i) {
	 if (rnd_01()<p_mutation) {
		(*(p->h1))[i]= !((*(p->h1))[i]);
		(*(p->mutations))[i]= true;
		++n_mutations;
		osm << "# (MUT " << id << ", " << i << ", " << role <<
		  ")" << endl;
		  // "  (h1,h2) -> (h1', h2) = (" <<
		  // ((*(p->h1))[i]?'0':'1') << ", " <<
		  // ((*(p->h2))[i]?'1':'0') << ") -> (" <<
		  // ((*(p->h1))[i]?'1':'0') << ", " <<
		  // ((*(p->h2))[i]?'1':'0') << ")" << endl;
	 }
  }
}

static void
recombinate_haplotype1(size_t id, size_t role, p_hap_t* p,
							  ostringstream& osm, size_t& n_recombinations,
							  double p_recombination) {
// Recombinations
  bool swap= false;
  for (size_t i= 1; i<len_hap; ++i) {
	 if ((*(p->h1))[i]!=(*(p->h2))[i]) {
		if (rnd_01()<p_recombination) {
		  swap= !swap;
		  ++n_recombinations;
		  osm << "# (REC " << id << ", " << i << ", " <<
			 role <<")" << endl;
		}
	 }
	 (*(p->recombinations))[i]= swap;
	 if (swap) {
		bool tmp= (*(p->h1))[i];
		(*(p->h1))[i]= (*(p->h2))[i];
		(*(p->h2))[i]= tmp;
	 }
  }
}

static void
random_swapping(p_hap_t* p) {
  // SWAP its haplotypes with probability 0.5
  if (rnd_01()<0.5) {
	 hap_t* temp= p->h1;
	 p->h1= p->h2;
	 p->h2= temp;
	 p->swapped= !p->swapped;
  }
}


static void
generate_rnd_trio(p_hap_t*& c, p_hap_t*& f, p_hap_t*& m,
						size_t id, ostringstream& osm,
						size_t& n_mutations, double p_mutation,
						size_t& n_recombinations, double p_recombination) {
  f= new p_hap_t();
  f->swapped= false;
  m= new p_hap_t();
  m->swapped= false;
  f->h1= new hap_t(*(c->h1));
  m->h1= new hap_t(*(c->h2));
  f->h2= generate_rnd_hap();
  m->h2= generate_rnd_hap();
  f->mutations= new hap_t(len_hap, false);
  m->mutations= new hap_t(len_hap, false);
  f->recombinations= new hap_t(len_hap, false);
  m->recombinations= new hap_t(len_hap, false);
  mutate_haplotype1(id, 0, f, osm, n_mutations, p_mutation);
  recombinate_haplotype1(id, 0, f, osm, n_recombinations, p_recombination);
  random_swapping(f);
  mutate_haplotype1(id, 1, m, osm, n_mutations, p_mutation);
  recombinate_haplotype1(id, 1, m, osm, n_recombinations, p_recombination);
  random_swapping(m);
}

static void
prepare_founder(p_hap_t*& p) {
  bool swap= p->swapped;
  for (size_t i= 0; i<len_hap; ++i) {
	 swap= p->swapped ^ (*(p->recombinations))[i];
	 if ((*(p->mutations))[i]) {
		if (swap) {
		  (*(p->h1))[i]= (*(p->h2))[i];
		} else {
		  (*(p->h2))[i]= (*(p->h1))[i];
		}
	 }
  }
}

int main() {
#ifdef _MOD_RND
  INFO("Generation of a random pedigree");
#else // _MOD_FULL
  INFO("Generation of a full binary pedigree of depth %2.0f.", MAX_LEVEL);
#endif
  PRINT_SYSTEM_INFORMATIONS;

  time_t rawtime;
  struct tm * timeinfo;
  int base_time= ((int)time (&rawtime));

#ifdef EXT_SEED
  SEED= EXT_SEED;
#else
  SEED= base_time;
#endif
  srand(SEED);

  INFO("Seed=%16d \t Time=%16d", SEED, base_time);

//Mutation data
  ostringstream osm;
  size_t n_mutations= 0;
  size_t n_recombinations= 0;

  list<p_hap_t*> happ;
  list<trio> trios;

  p_hap_t* p;
  p= generate_rnd_p_hap();
  happ.push_back(p);

  p_hap_t* f;
  p_hap_t* m;
  generate_rnd_trio(p, f, m, 0, osm,
						  n_mutations, p_mutation,
						  n_recombinations, p_recombination);
  happ.push_back(f);
  happ.push_back(m);
  trios.push_back(trio(1, 2, 0));

  list<p_hap_t*> q_i;
  list<double> probs;
  list<int> ids;
  q_i.push_back(f);
  q_i.push_back(m);
  probs.push_back(1.0);
  probs.push_back(1.0);
  ids.push_back(1);
  ids.push_back(2);

  int curr_id= 3;

  while (!q_i.empty()) {
	 p= q_i.front();
	 q_i.pop_front();
	 int id= ids.front();
	 ids.pop_front();
	 double prob= probs.front();
	 probs.pop_front();
#ifdef _MOD_RND
	 DEBUG("Extracted individual %4d with probability %.5e.", id, prob);
	 if (rnd_01()<prob) {
#else // _MOD_FULL
	 DEBUG("Extracted individual %4d of level %2.0f.", id, prob);
	 if (prob<MAX_LEVEL) {
#endif
		DEBUG("Individual %4d selected to be a children.", id);
		generate_rnd_trio(p, f, m, id, osm,
								n_mutations, p_mutation,
								n_recombinations, p_recombination);
		happ.push_back(f);
		happ.push_back(m);
		trios.push_back(trio(curr_id, curr_id+1, id));
		q_i.push_back(f);
		ids.push_back(curr_id);
		q_i.push_back(m);
		ids.push_back(curr_id+1);
#ifdef _MOD_RND
		probs.push_back(prob*probability_factor);
		probs.push_back(prob*probability_factor);
#else // _MOD_FULL
		probs.push_back(prob+1);
		probs.push_back(prob+1);
#endif
		curr_id= curr_id+2;
	 } else {
		DEBUG("Individual %d selected to be a founder.", id);
		DEBUG("Since individual %d has been selected to be a founder we set it as all-homozygote.", id);
		prepare_founder(p);
	 }
  }
// Salva i risultati
  ostringstream oss;
  oss << base_time;
  timeinfo = localtime ( &rawtime );

  ofstream fout(string(string("gp-mutated-")+oss.str()+string(".txt")).c_str());
#ifdef _MOD_RND
  fout << "# Random pedigree generated " << asctime (timeinfo);
#else // _MOD_FULL
  fout << "# Full pedigree generated " << asctime (timeinfo);
#endif
  fout << "# MUTATED with probability " << p_mutation << endl;
  fout << "# RECOMBINATED with probability " << p_recombination << endl;
  fout << "# N. individuals= " << curr_id << "   N. loci= " << len_hap << endl;
  fout << "# Seed= " << SEED << endl;
  fout << len_hap << endl;
  fout << curr_id << endl;
  fout << "# genotypes" << endl;
  size_t j= 0;
  for (list<p_hap_t*>::iterator it= happ.begin();
		 it!= happ.end();
		 ++it) {
	 fout << "# " << j << "a) "<< *((*it)->h1) << endl;
	 fout << "# " << j << "b) "<< *((*it)->h2) << endl;
	 fout << **it << endl;
	 ++j;
  }
  fout << "# mutations " << endl;
  fout << osm.str();
  fout << "# mutations end. total number "<< n_mutations<< endl;
  fout << "# recombinations end. total number "<< n_recombinations << endl;
  fout << "# pedigree" << endl;
  for (list<trio>::const_iterator it= trios.begin();
		 it!= trios.end();
		 ++it) {
	 fout << *it << endl;
  }
  fout << "# end" << endl;
  fout.close();

  INFO("Generation complete. Seed=%16d \t Time=%16d", SEED, base_time);
  INFO("No. of mutations=     " ST_FMTL(4), n_mutations);
  INFO("No. of recombinations=" ST_FMTL(4), n_recombinations);

}
