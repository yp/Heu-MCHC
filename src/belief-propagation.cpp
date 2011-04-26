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
#include "belief-propagation.hpp"

#include "log.h"
#include "util.h"

#define MIN_DEV (1E-10)

size_t row_h::cont= 0;

size_t col_h::cont= 0;

struct int_cmp {
  bool operator()(const int* a, const int *b) const {
	 return *a<*b;
  }
};


#define TANH_ONE (1.9061547465e+01)

#if defined(USE_MKL) && defined(USE_FRAMEWAVE)
#error Invalid settings of preprocessor flags. You must specify one of USE_MKL and USE_FRAMEWAVE
#endif

#ifdef USE_MKL

#include <mkl_vml_functions.h>
#include <mkl_cblas.h>

// **************************************************************
// *
// *  MKL VERSION
// *
// **************************************************************/
void
belief_propagation(system_desc& sd, pbp_config bp) {
// Count the ones
  const size_t n_rows= sd.n_rows();
  const size_t n_cols= sd.n_cols();
  const size_t nc= n_cols-n_rows;
  const double GAMMA= bp->gamma;
  const double ONE_MINUS_GAMMA= 1.0-bp->gamma;
  const double DOUBLE_ONE_MINUS_GAMMA= 2.0 * ONE_MINUS_GAMMA;
  size_t* rowslen= NPALLOC(size_t, n_rows+1);
  size_t abs_posr= 0;
  size_t r= 0;
  for (system_desc::rows_node* rn= sd.rows_begin();
		 rn != sd.rows_end();
		 rn= rn->next()) {
	 rowslen[r]= abs_posr;
	 abs_posr+= (**rn).size()-1;
	 (**rn).header()->pos= r;
	 ++r;
  }
  rowslen[r]= abs_posr;

  size_t* colslen= NPALLOC(size_t, n_cols+1);
  size_t abs_posc= 0;
  size_t c= 0;
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end() && !(**cn).header()->phantom;
		 cn= cn->next()){
	 colslen[c]= abs_posc;
	 abs_posc+= (**cn).size();
	 (**cn).header()->pos= c;
	 ++c;
  }
  colslen[c]= abs_posc;
  my_assert(abs_posc == abs_posr);

  BP_T* const rm= NPALLOC(BP_T, abs_posc);
  BP_T* const rmt= NPALLOC(BP_T, abs_posc);
  BP_T* const qm= NPALLOC(BP_T, abs_posc);
  BP_T* const qmt= NPALLOC(BP_T, abs_posc);
  BP_T* const p1= NPALLOC(BP_T, nc);
  BP_T* const q1= NPALLOC(BP_T, nc);
  BP_T* const llv= NPALLOC(BP_T, nc);
  int* const synd= NPALLOC(int, n_rows);
  size_t** const asdr= NPALLOC(size_t*, n_rows);
  size_t** const asdc= NPALLOC(size_t*, nc);
  r= 0;
  for (system_desc::rows_node* rs= sd.rows_begin();
		 rs != sd.rows_end();
		 rs= rs->next(), ++r) {
	 asdr[r]= NPALLOC(size_t, rowslen[r+1]-rowslen[r]);
	 size_t v= 0;
	 for (system_desc::row_node* rn= (**rs).begin();
			rn != (**rs).end() && !(**rn).col().header()->phantom;
			rn= rn->next(), ++v) {
		asdr[r][v]= (**rn).col().header()->pos;
	 }
	 synd[r]= (**(**rs).end()->prev()).col().header()->assigned?-1:1;
  }
  c=0;
  for (system_desc::cols_node* cs= sd.cols_begin();
		 cs != sd.cols_end() && !(**cs).header()->phantom;
		 cs= cs->next(), ++c) {
	 asdc[c]= NPALLOC(size_t, colslen[c+1]-colslen[c]);
	 size_t v= 0;
	 for (system_desc::col_node* cn= (**cs).begin();
			cn != (**cs).end(); cn= cn->next(), ++v) {
		asdc[c][v]= (**cn).row().header()->pos;
	 }
	 p1[c]= (**cs).header()->p1;
  }

  size_t* cols= NPALLOC(size_t, n_rows);
  size_t* rows= NPALLOC(size_t, nc);
//Inizializzazione
  for (size_t i= 0; i<n_rows; cols[i]=0, ++i) {};
  for (size_t var= 0; var<nc; ++var) {
	 q1[var]= (BP1-p1[var])/p1[var];
  }
  vdLn(nc, q1, q1);
  cblas_dcopy(nc, q1, 1, llv, 1);
  cblas_dscal(nc, BP0_5, llv, 1);
  vdTanh(nc, llv, llv);
  for (size_t var= 0; var<nc; ++var) {
	 const BP_T ll= q1[var];
	 const BP_T llt= llv[var];
	 llv[var]= ll;
	 for (size_t cons= 0; cons<colslen[var+1]-colslen[var]; ++cons) {
		const size_t base_pos= asdc[var][cons];
		qm[rowslen[base_pos]+cols[base_pos]]= ll;
		qmt[rowslen[base_pos]+cols[base_pos]]= llt;
		++cols[base_pos];
	 }
  }
//Ciclo
  size_t iter= 0;
  BP_T tot_dev= BP0;
  do {
// Calcolo rm
	 for (size_t i= 0; i<nc; rows[i]=0, ++i) {};
	 for (size_t cons= 0; cons<n_rows; ++cons) {
		for (size_t var= 0; var<rowslen[cons+1]-rowslen[cons]; ++var) {
		  const size_t pos= colslen[asdr[cons][var]]+rows[asdr[cons][var]];
		  BP_T tmp= BP1*synd[cons];
		  const size_t v1start= rowslen[cons];
		  const size_t v1end1= rowslen[cons]+var;
		  const size_t v1end2= rowslen[cons+1];
		  for (size_t v1= v1start; v1<v1end1; ++v1) {
			 tmp*= qmt[v1];
		  }
		  for (size_t v1= v1end1+1; v1<v1end2; ++v1) {
			 tmp*= qmt[v1];
		  }
		  rmt[pos]= tmp;
		  ++rows[asdr[cons][var]];
//		  fprintf(stderr, "rm %5d %5d = %+e\n", cons, asdr[cons][var], new_rm);
		}
	 }
	 if (iter>0) {
		vdAtanh(abs_posc, rmt, rmt);
		for (size_t i= 0; i<abs_posc; ++i) {
		  if (!isfinite(rmt[i])) {
			 const BP_T sgn= (signbit(rmt[i])!=0)?-BP1:BP1;
			 rmt[i]= sgn*TANH_ONE;
		  }
		}
		cblas_dscal(abs_posc, GAMMA, rm, 1);
		cblas_daxpy(abs_posc, DOUBLE_ONE_MINUS_GAMMA, rmt, 1, rm, 1);
	 } else {
		vdAtanh(abs_posc, rmt, rm);
		for (size_t i= 0; i<abs_posc; ++i) {
		  if (!isfinite(rm[i])) {
			 const BP_T sgn= (signbit(rm[i])!=0)?-BP1:BP1;
			 rm[i]= sgn*TANH_ONE;
		  }
		}
		cblas_dscal(abs_posc, 2.0, rm, 1);
	 }
// Calcolo qm
	 for (size_t i= 0; i<n_rows; cols[i]=0, ++i) {};
	 for (size_t var= 0; var<nc; ++var) {
		const BP_T init_tmp= log((BP1-p1[var])/p1[var]);
		for (size_t cons= 0; cons<colslen[var+1]-colslen[var]; ++cons) {
		  BP_T tmp= init_tmp;
		  const size_t c1start= colslen[var];
		  const size_t c1end= colslen[var]+cons;
		  const size_t c1end2= colslen[var+1];
		  for (size_t c1= c1start; c1<c1end; ++c1) {
			 tmp+= rm[c1];
		  }
		  for (size_t c1= c1end+1; c1<c1end2; ++c1) {
			 tmp+= rm[c1];
		  }
		  const size_t pos= rowslen[asdc[var][cons]]+cols[asdc[var][cons]];
		  const BP_T new_qm= (iter>0)?
			 (ONE_MINUS_GAMMA * tmp) + (GAMMA * qm[pos])
			 :tmp;
		  qm[pos]= new_qm;
		  ++cols[asdc[var][cons]];
//		  fprintf(stderr, "qm %5d %5d = %+e\n", asdc[var][cons], var, new_qm);
		  fail_if(!isfinite(new_qm));
		}
	 }
	 cblas_dcopy(abs_posc, qm, 1, qmt, 1);
	 cblas_dscal(abs_posc, BP0_5, qmt, 1);
	 vdTanh(abs_posc, qmt, qmt);

// Calcolo posteriori
	 tot_dev= BP0;
	 for (size_t var= 0; var<nc; ++var) {
		BP_T tmp= llv[var];
		const size_t c1start= colslen[var];
		const size_t c1end= colslen[var+1];
		for (size_t c1= c1start; c1<c1end; ++c1) {
		  tmp+= rm[c1];
		}
		const BP_T new_q1= (tmp*ONE_MINUS_GAMMA)+(q1[var]*GAMMA);
		tot_dev+= BPABS(q1[var]-new_q1);
		q1[var]= new_q1;
		fail_if(!isfinite(tmp));
//		fprintf(stderr, "%+.6e (%+.2e)  ", tmp, p1[var]);
	 }
//	 fprintf(stderr, "\n");
	 ++iter;
  } while (iter< bp->max_iter &&
			  (tot_dev>MIN_DEV || iter==1 || !isfinite(tot_dev)));
  c= 0;
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end() && !(**cn).header()->phantom;
		 cn= cn->next()){
	 (**cn).header()->q1= q1[c];
	 ++c;
  }
  for (size_t cons= 0; cons<n_rows; ++cons) {
	 pfree(asdr[cons]);
  }
  for (size_t var= 0; var<nc; ++var) {
	 pfree(asdc[var]);
  }
  pfree(rows);
  pfree(cols);
  pfree(rm);
  pfree(rmt);
  pfree(qm);
  pfree(qmt);
  pfree(p1);
  pfree(q1);
  pfree(llv);
  pfree(asdr);
  pfree(synd);
  pfree(asdc);
  pfree(colslen);
  pfree(rowslen);
}


#elif defined(USE_FRAMEWAVE)


// **************************************************************
// *
// *  FRAMEWAVE VERSION
// *
// **************************************************************/

#include <fwBase.h>
#include <fwSignal.h>

void
belief_propagation(system_desc& sd, pbp_config bp) {
// Count the ones
  const size_t n_rows= sd.n_rows();
  const size_t n_cols= sd.n_cols();
  const size_t nc= n_cols-n_rows;
  const double GAMMA= bp->gamma;
  const double ONE_MINUS_GAMMA= 1.0-bp->gamma;
  const double DOUBLE_ONE_MINUS_GAMMA= 2.0 * ONE_MINUS_GAMMA;
  size_t* rowslen= NPALLOC(size_t, n_rows+1);
  size_t abs_posr= 0;
  size_t r= 0;
  for (system_desc::rows_node* rn= sd.rows_begin();
		 rn != sd.rows_end();
		 rn= rn->next()) {
	 rowslen[r]= abs_posr;
	 abs_posr+= (**rn).size()-1;
	 (**rn).header()->pos= r;
	 ++r;
  }
  rowslen[r]= abs_posr;

  size_t* colslen= NPALLOC(size_t, n_cols+1);
  size_t abs_posc= 0;
  size_t c= 0;
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end() && !(**cn).header()->phantom;
		 cn= cn->next()){
	 colslen[c]= abs_posc;
	 abs_posc+= (**cn).size();
	 (**cn).header()->pos= c;
	 ++c;
  }
  colslen[c]= abs_posc;
  my_assert(abs_posc == abs_posr);

  Fw64f* const rm= fwsMalloc_64f(abs_posc);
  Fw64f* const rmt= fwsMalloc_64f(abs_posc);
  Fw64f* const qm= fwsMalloc_64f(abs_posc);
  Fw64f* const qmt= fwsMalloc_64f(abs_posc);
  Fw64f* const p1= fwsMalloc_64f(nc);
  Fw64f* const q1= fwsMalloc_64f(nc);
  Fw64f* const llv= fwsMalloc_64f(nc);
  int* const synd= NPALLOC(int, n_rows);
  size_t** const asdr= NPALLOC(size_t*, n_rows);
  size_t** const asdc= NPALLOC(size_t*, nc);
  r= 0;
  for (system_desc::rows_node* rs= sd.rows_begin();
		 rs != sd.rows_end();
		 rs= rs->next(), ++r) {
	 asdr[r]= NPALLOC(size_t, rowslen[r+1]-rowslen[r]);
	 size_t v= 0;
	 for (system_desc::row_node* rn= (**rs).begin();
			rn != (**rs).end() && !(**rn).col().header()->phantom;
			rn= rn->next(), ++v) {
		asdr[r][v]= (**rn).col().header()->pos;
	 }
	 synd[r]= (**(**rs).end()->prev()).col().header()->assigned?-1:1;
  }
  c=0;
  for (system_desc::cols_node* cs= sd.cols_begin();
		 cs != sd.cols_end() && !(**cs).header()->phantom;
		 cs= cs->next(), ++c) {
	 asdc[c]= NPALLOC(size_t, colslen[c+1]-colslen[c]);
	 size_t v= 0;
	 for (system_desc::col_node* cn= (**cs).begin();
			cn != (**cs).end(); cn= cn->next(), ++v) {
		asdc[c][v]= (**cn).row().header()->pos;
	 }
	 p1[c]= (**cs).header()->p1;
  }

  size_t* cols= NPALLOC(size_t, n_rows);
  size_t* rows= NPALLOC(size_t, nc);
//Inizializzazione
  for (size_t i= 0; i<n_rows; cols[i]=0, ++i) {};
  for (size_t var= 0; var<nc; ++var) {
	 q1[var]= (BP1-p1[var])/p1[var];
  }
  fwsLn_64f_A50(q1, q1, nc);
  fwsMulC_64f(q1, BP0_5, llv, nc);
  fwsTanh_64f_A50(llv, llv, nc);
  for (size_t var= 0; var<nc; ++var) {
	 const Fw64f ll= q1[var];
	 const Fw64f llt= llv[var];
	 llv[var]= ll;
	 for (size_t cons= 0; cons<colslen[var+1]-colslen[var]; ++cons) {
		const size_t base_pos= asdc[var][cons];
		qm[rowslen[base_pos]+cols[base_pos]]= ll;
		qmt[rowslen[base_pos]+cols[base_pos]]= llt;
		++cols[base_pos];
	 }
  }
//Ciclo
  size_t iter= 0;
  Fw64f tot_dev= BP0;
  do {
// Calcolo rm
	 for (size_t i= 0; i<nc; rows[i]=0, ++i) {};
	 for (size_t cons= 0; cons<n_rows; ++cons) {
		for (size_t var= 0; var<rowslen[cons+1]-rowslen[cons]; ++var) {
		  const size_t pos= colslen[asdr[cons][var]]+rows[asdr[cons][var]];
		  Fw64f tmp= BP1*synd[cons];
		  const size_t v1start= rowslen[cons];
		  const size_t v1end1= rowslen[cons]+var;
		  const size_t v1end2= rowslen[cons+1];
		  for (size_t v1= v1start; v1<v1end1; ++v1) {
			 tmp*= qmt[v1];
		  }
		  for (size_t v1= v1end1+1; v1<v1end2; ++v1) {
			 tmp*= qmt[v1];
		  }
		  rmt[pos]= tmp;
		  ++rows[asdr[cons][var]];
//		  fprintf(stderr, "rm %5d %5d = %+e\n", cons, asdr[cons][var], new_rm);
		}
	 }
	 if (iter>0) {
		fwsAtanh_64f_A50(rmt, rmt, abs_posc);
		for (size_t i= 0; i<abs_posc; ++i) {
		  if (!isfinite(rmt[i])) {
			 const Fw64f sgn= (signbit(rmt[i])!=0)?-BP1:BP1;
			 rmt[i]= sgn*TANH_ONE;
		  }
		}
		fwsMulC_64f_I(GAMMA, rm, abs_posc);
		fwsMulC_64f_I(DOUBLE_ONE_MINUS_GAMMA, rmt, abs_posc);
		fwsAdd_64f_I(rmt, rm, abs_posc);
	 } else {
		fwsAtanh_64f_A50(rmt, rm, abs_posc);
		for (size_t i= 0; i<abs_posc; ++i) {
		  if (!isfinite(rm[i])) {
			 const Fw64f sgn= (signbit(rm[i])!=0)?-BP1:BP1;
			 rm[i]= sgn*TANH_ONE;
		  }
		}
		fwsMulC_64f_I(2.0, rm, abs_posc);
	 }
// Calcolo qm
	 for (size_t i= 0; i<n_rows; cols[i]=0, ++i) {};
	 for (size_t var= 0; var<nc; ++var) {
		const Fw64f init_tmp= log((BP1-p1[var])/p1[var]);
		for (size_t cons= 0; cons<colslen[var+1]-colslen[var]; ++cons) {
		  Fw64f tmp= init_tmp;
		  const size_t c1start= colslen[var];
		  const size_t c1end= colslen[var]+cons;
		  const size_t c1end2= colslen[var+1];
		  for (size_t c1= c1start; c1<c1end; ++c1) {
			 tmp+= rm[c1];
		  }
		  for (size_t c1= c1end+1; c1<c1end2; ++c1) {
			 tmp+= rm[c1];
		  }
		  const size_t pos= rowslen[asdc[var][cons]]+cols[asdc[var][cons]];
		  const Fw64f new_qm= (iter>0)?
			 (ONE_MINUS_GAMMA * tmp) + (GAMMA * qm[pos])
			 :tmp;
		  qm[pos]= new_qm;
		  ++cols[asdc[var][cons]];
//		  fprintf(stderr, "qm %5d %5d = %+e\n", asdc[var][cons], var, new_qm);
		  fail_if(!isfinite(new_qm));
		}
	 }
	 fwsMulC_64f(qm, BP0_5, qmt, abs_posc);
	 fwsTanh_64f_A50(qmt, qmt, abs_posc);

// Calcolo posteriori
	 tot_dev= BP0;
	 for (size_t var= 0; var<nc; ++var) {
		Fw64f tmp= llv[var];
		const size_t c1start= colslen[var];
		const size_t c1end= colslen[var+1];
		for (size_t c1= c1start; c1<c1end; ++c1) {
		  tmp+= rm[c1];
		}
		const Fw64f new_q1= (tmp*ONE_MINUS_GAMMA)+(q1[var]*GAMMA);
		tot_dev+= BPABS(q1[var]-new_q1);
		q1[var]= new_q1;
		fail_if(!isfinite(tmp));
//		fprintf(stderr, "%+.6e (%+.2e)  ", tmp, p1[var]);
	 }
//	 fprintf(stderr, "\n");
	 ++iter;
  } while (iter< bp->max_iter &&
			  (tot_dev>MIN_DEV || iter==1 || !isfinite(tot_dev)));
  c= 0;
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end() && !(**cn).header()->phantom;
		 cn= cn->next()){
	 (**cn).header()->q1= q1[c];
	 ++c;
  }
  for (size_t cons= 0; cons<n_rows; ++cons) {
	 pfree(asdr[cons]);
  }
  for (size_t var= 0; var<nc; ++var) {
	 pfree(asdc[var]);
  }
  pfree(rows);
  pfree(cols);
  pfree(rm);
  pfree(rmt);
  pfree(qm);
  pfree(qmt);
  pfree(p1);
  pfree(q1);
  pfree(llv);
  pfree(asdr);
  pfree(synd);
  pfree(asdc);
  pfree(colslen);
  pfree(rowslen);
}

#else // not defined USE_MKL and not defined USE_FRAMEWAVE

// **************************************************************
// *
// *  PLAIN VERSION
// *
// **************************************************************/

void
belief_propagation(system_desc& sd, pbp_config bp) {
// Count the ones
  const size_t n_rows= sd.n_rows();
  const size_t n_cols= sd.n_cols();
  const size_t nc= n_cols-n_rows;
  const double GAMMA= bp->gamma;
  size_t* rowslen= NPALLOC(size_t, n_rows+1);
  size_t abs_posr= 0;
  size_t r= 0;
  for (system_desc::rows_node* rn= sd.rows_begin();
		 rn != sd.rows_end();
		 rn= rn->next()) {
	 rowslen[r]= abs_posr;
	 abs_posr+= (**rn).size()-1;
	 (**rn).header()->pos= r;
	 ++r;
  }
  rowslen[r]= abs_posr;

  size_t* colslen= NPALLOC(size_t, n_cols+1);
  size_t abs_posc= 0;
  size_t c= 0;
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end() && !(**cn).header()->phantom;
		 cn= cn->next()){
	 colslen[c]= abs_posc;
	 abs_posc+= (**cn).size();
	 (**cn).header()->pos= c;
	 ++c;
  }
  colslen[c]= abs_posc;
  my_assert(abs_posc == abs_posr);

  BP_T* const rm= NPALLOC(BP_T, abs_posc);
  BP_T* const qm= NPALLOC(BP_T, abs_posc);
  BP_T* const qmt= NPALLOC(BP_T, abs_posc);
  BP_T* const p1= NPALLOC(BP_T, nc);
  BP_T* const q1= NPALLOC(BP_T, nc);
  BP_T* const llv= NPALLOC(BP_T, nc);
  int* const synd= NPALLOC(int, n_rows);
  size_t** const asdr= NPALLOC(size_t*, n_rows);
  size_t** const asdc= NPALLOC(size_t*, nc);
  r= 0;
  for (system_desc::rows_node* rs= sd.rows_begin();
		 rs != sd.rows_end();
		 rs= rs->next(), ++r) {
	 asdr[r]= NPALLOC(size_t, rowslen[r+1]-rowslen[r]);
	 size_t v= 0;
	 for (system_desc::row_node* rn= (**rs).begin();
			rn != (**rs).end() && !(**rn).col().header()->phantom;
			rn= rn->next(), ++v) {
		asdr[r][v]= (**rn).col().header()->pos;
	 }
	 synd[r]= (**(**rs).end()->prev()).col().header()->assigned?-1:1;
  }
  c=0;
  for (system_desc::cols_node* cs= sd.cols_begin();
		 cs != sd.cols_end() && !(**cs).header()->phantom;
		 cs= cs->next(), ++c) {
	 asdc[c]= NPALLOC(size_t, colslen[c+1]-colslen[c]);
	 size_t v= 0;
	 for (system_desc::col_node* cn= (**cs).begin();
			cn != (**cs).end(); cn= cn->next(), ++v) {
		asdc[c][v]= (**cn).row().header()->pos;
	 }
	 p1[c]= (**cs).header()->p1;
  }

  size_t* cols= NPALLOC(size_t, n_rows);
  size_t* rows= NPALLOC(size_t, nc);
//Inizializzazione
  for (size_t i= 0; i<n_rows; cols[i]=0, ++i) {};
  for (size_t var= 0; var<nc; ++var) {
	 const BP_T ll= log((BP1-p1[var])/p1[var]);
	 const BP_T llt= tanh(BP0_5*ll);
	 q1[var]= ll;
	 llv[var]= ll;
	 for (size_t cons= 0; cons<colslen[var+1]-colslen[var]; ++cons) {
		qm[rowslen[asdc[var][cons]]+cols[asdc[var][cons]]]= ll;
		qmt[rowslen[asdc[var][cons]]+cols[asdc[var][cons]]]= llt;
		++cols[asdc[var][cons]];
	 }
  }
//Ciclo
  size_t iter= 0;
  BP_T tot_dev= BP0;
  do {
// Calcolo rm
	 for (size_t i= 0; i<nc; rows[i]=0, ++i) {};
	 for (size_t cons= 0; cons<n_rows; ++cons) {
		for (size_t var= 0; var<rowslen[cons+1]-rowslen[cons]; ++var) {
		  BP_T tmp= BP1*synd[cons];
		  const size_t v1start= rowslen[cons];
		  const size_t v1end1= rowslen[cons]+var;
		  const size_t v1end2= rowslen[cons+1];
		  for (size_t v1= v1start; v1<v1end1; ++v1) {
			 tmp*= qmt[v1];
		  }
		  for (size_t v1= v1end1+1; v1<v1end2; ++v1) {
			 tmp*= qmt[v1];
		  }
		  if (tmp<=-1.0)
			 tmp= 2.0*-TANH_ONE;
		  else if (tmp>=1.0)
			 tmp= 2.0*TANH_ONE;
		  else
			 tmp= 2.0*atanh(tmp);
		  const BP_T new_rm= (iter>0)?
			 ((1.0-GAMMA)*tmp)+
			 (GAMMA*rm[colslen[asdr[cons][var]]+rows[asdr[cons][var]]])
			 :tmp;
		  rm[colslen[asdr[cons][var]]+rows[asdr[cons][var]]]= new_rm;
		  ++rows[asdr[cons][var]];
//		  fprintf(stderr, "rm %5d %5d = %+e\n", cons, asdr[cons][var], new_rm);
		  fail_if(!isfinite(new_rm));
		}
	 }
// Calcolo qm
	 for (size_t i= 0; i<n_rows; cols[i]=0, ++i) {};
	 for (size_t var= 0; var<nc; ++var) {
		const BP_T init_tmp= log((BP1-p1[var])/p1[var]);
		for (size_t cons= 0; cons<colslen[var+1]-colslen[var]; ++cons) {
		  BP_T tmp= init_tmp;
		  const size_t c1start= colslen[var];
		  const size_t c1end= colslen[var]+cons;
		  const size_t c1end2= colslen[var+1];
		  for (size_t c1= c1start; c1<c1end; ++c1) {
			 tmp+= rm[c1];
		  }
		  for (size_t c1= c1end+1; c1<c1end2; ++c1) {
			 tmp+= rm[c1];
		  }
		  const BP_T new_qm= (iter>0)?
			 ((1.0-GAMMA) * tmp)
			 + (GAMMA * qm[rowslen[asdc[var][cons]]+cols[asdc[var][cons]]])
			 :tmp;
		  qm[rowslen[asdc[var][cons]]+cols[asdc[var][cons]]]= new_qm;
		  ++cols[asdc[var][cons]];
//		  fprintf(stderr, "qm %5d %5d = %+e\n", asdc[var][cons], var, new_qm);
		  fail_if(!isfinite(new_qm));
		}
	 }
	 for (size_t i= 0; i<abs_posc; ++i) qmt[i]= BP0_5*qm[i];
	 for (size_t i= 0; i<abs_posc; ++i) qmt[i]= tanh(qmt[i]);

// Calcolo posteriori
	 tot_dev= BP0;
	 for (size_t var= 0; var<nc; ++var) {
		BP_T tmp= llv[var];
		const size_t c1start= colslen[var];
		const size_t c1end= colslen[var+1];
		for (size_t c1= c1start; c1<c1end; ++c1) {
		  tmp+= rm[c1];
		}
		const BP_T new_q1= (tmp*(1.0-GAMMA))+(q1[var]*GAMMA);
		tot_dev+= BPABS(q1[var]-new_q1);
		q1[var]= new_q1;
		fail_if(!isfinite(tmp));
//		fprintf(stderr, "%+.6e (%+.2e)  ", tmp, p1[var]);
	 }
//	 fprintf(stderr, "\n");
	 ++iter;
  } while (iter< bp->max_iter &&
			  (tot_dev>MIN_DEV || iter==1 || !isfinite(tot_dev)));
  c= 0;
  for (system_desc::cols_node* cn= sd.cols_begin();
		 cn != sd.cols_end() && !(**cn).header()->phantom;
		 cn= cn->next()){
	 (**cn).header()->q1= q1[c];
	 ++c;
  }
  for (size_t cons= 0; cons<n_rows; ++cons) {
	 pfree(asdr[cons]);
  }
  for (size_t var= 0; var<nc; ++var) {
	 pfree(asdc[var]);
  }
  pfree(rows);
  pfree(cols);
  pfree(rm);
  pfree(qm);
  pfree(qmt);
  pfree(p1);
  pfree(q1);
  pfree(llv);
  pfree(asdr);
  pfree(synd);
  pfree(asdc);
  pfree(colslen);
  pfree(rowslen);
}

#endif

