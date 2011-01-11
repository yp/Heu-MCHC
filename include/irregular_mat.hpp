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
#ifndef _IRREGULAR_MAT_HPP_
#define _IRREGULAR_MAT_HPP_

#include "util.h"

template <typename T>
class matrix {
private:
  const size_t R;
  const size_t C;
  T* mat;
public:
  matrix(const size_t _R, const size_t _C, const T& initial_value):
		R(_R), C(_C), mat(NULL)
  {
	 my_assert(0<R);
	 my_assert(0<C);
	 mat= NPALLOC(T, R*C);
	 
	 for (size_t i= 0; i < R*C; ++i)
		mat[i]= initial_value;
  }

  ~matrix() {
	 pfree(mat);
  }

  T& operator()(const size_t row, const size_t col) {
	 my_assert(row<R);
	 my_assert(col<C);
	 return mat[C*row+col];
  }

};

template <typename T>
class vect {
private:
  const size_t N;
  T* mat;
public:
  vect(const size_t _N, const T& initial_value):
		N(_N), mat(NULL)
  {
	 my_assert(0<N);
	 mat= new T[N];
	 for (size_t i= 0; i < N; ++i)
		mat[i]= initial_value;
  }

  vect(const vect<T>& orig):
		N(orig.N), mat(NULL)
  {
	 mat= new T[N];
	 for (size_t i= 0; i < N; ++i)
		mat[i]= orig.mat[i];
  }

  ~vect() {
	 delete [] mat;
  }

  T& operator()(const size_t i) {
	 my_assert(i<N);
	 return mat[i];
  }

  const T& operator()(const size_t i) const {
	 my_assert(i<N);
	 return mat[i];
  }

  size_t size(void) const {
	 return N;
  }

  void set_all(const T& value) {
	 for(size_t i= 0; i<N; ++i)
		mat[i]= value;
  }
};


typedef vect<size_t> uint_vec;
typedef matrix<size_t> uint_mat;


template <typename T>
class irregular_mat {
private:
  uint_vec lens;
  T** mat;

public:
  irregular_mat(const uint_vec& len_rows, const T& initial_value)
		:lens(len_rows), mat(NULL)
  {
	 mat= new T*[len_rows.size()];
	 for (size_t i= 0; i<len_rows.size(); ++i) {
		mat[i]= new T[len_rows(i)];
		for (size_t j= 0; j<len_rows(i); ++j)
		  mat[i][j]= initial_value;
	 }
  }

  ~irregular_mat() {
	 for (size_t i= 0; i<lens.size(); ++i) {
		delete [] mat[i];
	 }
	 delete [] mat;
  }

  T& operator()(const size_t row, const size_t col) {
	 my_assert(row<lens.size());
	 my_assert(col<lens(row));
	 return (mat[row])[col];
  }

  const T get(const size_t row, const size_t col) const {
	 my_assert(row<lens.size());
	 my_assert(col<lens(row));
	 return (mat[row])[col];
  }

  void set_all(const T& value) {
	 for (size_t i= 0; i<lens.size(); ++i) {
		for (size_t j= 0; j<lens(i); ++j)
		  mat[i][j]= value;
	 }
  }
};

template <class header_t, class data_t>
class list_t {
public:
  class node {
  private:
	 list_t* const l;
	 data_t* data_;
	 node* next_;
	 node* prev_;
  public:
	 node(list_t* const _l, data_t* _data)
		  :l(_l), data_(_data), next_(this), prev_(this)
	 {}

	 data_t& data(void) {
		return *data_;
	 }

	 data_t* pdata(void) const {
		return data_;
	 }

	 const data_t& operator*(void) {
		return *data_;
	 }

	 node* next(void) const {
		return next_;
	 }

	 node* prev(void) const {
		return prev_;
	 }

	 node* insert_after(data_t* new_data) {
		node* nn= new node(l, new_data);
		nn->next_= next_;
		nn->prev_= this;
		next_->prev_= nn;
		next_= nn;
		++l->size_;
		return nn;
	 }

	 node* insert_before(data_t* new_data) {
		return prev_->insert_after(new_data);
	 }

	 node* remove(void) {
		my_assert(next_ != this);
		node* n= next_;
		node* p= prev_;
		--l->size_;
		delete this;
		p->next_= n;
		n->prev_= p;
		return n;
	 }

	 template <typename data_free>
	 node* free_and_remove(const data_free& fr) {
		fr(data_);
		return remove();
	 }
  };

private:
  size_t size_;
  header_t header_;
  node* base;
public:
  list_t(const header_t& _header)
		:size_(0), header_(_header), base(new node(this, NULL))
  {
  }

  ~list_t(void) {
	 clear();
	 delete base;
  }

  const header_t& header(void) const {
	 return header_;
  }

  node* begin(void) const {
	 return base->next();
  }

  node* end(void) const {
	 return base;
  }

  size_t size(void) const {
	 return size_;
  }

  template <typename data_free>
  void free(const data_free& fr) {
	 while (begin()!=end())
		begin()->free_and_remove(fr);
  }

  void clear(void) {
	 while (begin()!=end())
		begin()->remove();
  }
};

template <class row_header_t, class col_header_t, class T>
class irr_mat {

public:
  class entry;

  typedef list_t<row_header_t, entry > row_t;
  typedef list_t<col_header_t, entry > col_t;

  typedef typename list_t<size_t, row_t>::node rows_node;
  typedef typename list_t<size_t, col_t>::node cols_node;

  typedef typename row_t::node row_node;
  typedef typename col_t::node col_node;

  class entry {
  private:
	 row_t& r_;
	 col_t& c_;
	 row_node* rn_;
	 col_node* cn_;
	 T* const data_;
  public:
	 entry(row_t& _r, col_t& _c, T* const _data)
		  :r_(_r), c_(_c), rn_(NULL), cn_(NULL), data_(_data)
	 {}

	 void set_position(row_node* _rn, col_node* _cn) {
		my_assert(rn_==NULL);
		my_assert(cn_==NULL);
		rn_= _rn;
		cn_= _cn;
	 }

	 T& data(void) const {
		return *data_;
	 }

	 const T& operator*(void) const {
		return *data_;
	 }

	 T* pdata(void) const {
		return data_;
	 }

	 row_t& row(void) const {
		return r_;
	 }

	 col_t& col(void) const {
		return c_;
	 }

	 row_node* row_n(void) const {
		return rn_;
	 }

	 col_node* col_n(void) const {
		return cn_;
	 }

  };


private:
  list_t<size_t, col_t> cols_;
  list_t<size_t, row_t> rows_;

public:
  irr_mat(void)
		:cols_(0), rows_(0)
  {}

  ~irr_mat(void)
  {
	 clear();
  }

  void insert_at(row_t& r, col_t& c, T* data) {
	 entry* e= new entry(r, c, data);
	 row_node* rn= r.end()->insert_before(e);
	 col_node* cn= c.end()->insert_before(e);
	 e->set_position(rn, cn);
  }

  rows_node* rows_begin(void) const {
	 return rows_.begin();
  }

  rows_node* rows_end(void) const {
	 return rows_.end();
  }

  rows_node* append_row(const row_header_t& rh) {
	 row_t* r= new row_t(rh);
	 return rows_.end()->insert_before(r);
  }

  cols_node* cols_begin(void) const {
	 return cols_.begin();
  }

  cols_node* cols_end(void) const {
	 return cols_.end();
  }

  cols_node* append_col(const col_header_t& ch) {
	 col_t* c= new col_t(ch);
	 return cols_.end()->insert_before(c);
  }

  cols_node* insert_col_before(const col_header_t& ch, cols_node* const cn) {
	 col_t* c= new col_t(ch);
	 return cn->insert_before(c);
  }

  size_t n_rows(void) const {
	 return rows_.size();
  }

  size_t n_cols(void) const {
	 return cols_.size();
  }


  void clear(void) {
	 rows_.clear();
	 cols_.clear();
  }

  template <typename data_free>
  void free(const data_free& fr) {
	 rows_node* r= rows_.begin();
	 while (r != rows_.end()) {
		row_node* n= (**r).begin();
		while (n != (**r).end()) {
		  fr((**n).pdata());
		  delete (*n).pdata();
		  n= n->next();
		}
		fr((**r).header());
		delete &**r;
		r= r->next();
	 }
	 rows_.clear();
	 cols_node* c= cols_.begin();
	 while (c != cols_.end()) {
		fr((**c).header());
		delete &**c;
		c= c->next();
	 }
	 cols_.clear();
  }


};

// #include <functional>

// struct intdel
//   :public unary_function<int*, void> {
//   void operator()(int* pt) const {
// 	 delete pt;
//   }
// };

// #define M 200
// int main() {
//   typedef irr_mat<int, int, int> imat;
//   int* h= new int[M];
//   for (int ci= 0; ci<M; ++ci)
// 	 h[ci]= ci;
//   imat m;
//   int ri;
//   for (int ci= 0; ci<M; ++ci)
// 	 m.append_col(h[ci]);
//   for (int ri= 0; ri<M; ++ri)
// 	 m.append_row(h[ri]);
//   int i= 0;
//   int ones= 0;
//   imat::cols_node* c;
//   imat::rows_node* r= m.rows_begin();
//   while (r != m.rows_end()) {
// 	 c= m.cols_begin();
// 	 while (c != m.cols_end()) {
// 		if (rand()>(RAND_MAX/2)) {
// 		  m.insert_at(r->data(), c->data(), new int(i));
// 		  printf("%3d ", i);
// 		  ++ones;
// 		} else {
// 		  printf("    ");
// 		}
// 		++i;
// 		c= c->next();
// 	 }
// 	 r= r->next();
// 	 printf("\n");
//   }
//   printf("\n");
//   r= m.rows_begin();
//   while (r != m.rows_end()) {
// 	 imat::row_node* n= (**r).begin();
// 	 while (n != (**r).end()) {
// 		printf("%3d ", (**n).data());
// 		n= n->next();
// 	 }
// 	 printf("\n");
// 	 r= r->next();
//   }
//   printf("\n");
//   c= m.cols_begin();
//   while (c != m.cols_end()) {
// 	 imat::col_node* n= (**c).begin();
// 	 while (n != (**c).end()) {
// 		printf("%3d ", (**n).data());
// 		n= n->next();
// 	 }
// 	 printf("\n");
// 	 c= c->next();
//   }
//   m.free(intdel());
//   delete [] h;
//   fprintf(stderr, "%d %e\n", ones, ((double)ones)/(M*M));
// }




#endif // _IRREGULAR_MAT_HPP_
