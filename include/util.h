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
/**
 *
 * @file util.h
 *
 * Funzioni di utilita' generale.
 *
 **/

#ifdef __cplusplus
using namespace std;
extern "C" {
#endif

#ifndef _UTIL_H_
#define _UTIL_H_


#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <stdbool.h>

#ifdef i386
#define ST_FMTL( l ) "%" #l "u"
#define ST_FMT "%u"
#else
#define ST_FMTL( l ) "%" #l "lu"
#define ST_FMT "%lu"
#endif

double rnd_01(void);

unsigned int
rnd_in_range(unsigned int min, unsigned int max);

#define NPALLOC( type, dim ) (type*)palloc((dim)*sizeof(type))
#define PALLOC( type ) (type*)palloc(sizeof(type))

#define _SET_VECTOR( V, val, lim )					\
  do {														\
	 for (int ___i= 0; ___i<lim; V[___i]=val, ++___i);	\
  } while (0)

#define MY_SWAP( type, el1, el2 ) \
  {										 \
	 type tmp= (el1);					 \
	 (el1)= el2;						 \
	 (el2)= tmp;						 \
  }

#define MIN( x, y ) (((x)<=(y))?(x):(y))
#define MAX( x, y ) (((x)>=(y))?(x):(y))

void* palloc(size_t dim) __attribute__ ((malloc));

void pfree(void* p)
#ifndef __ICC
  __attribute__ ((nonnull))
#endif
  ;

char* c_palloc(size_t dim) __attribute__ ((malloc));

char* alloc_and_copy(char* source)  __attribute__ ((malloc));

void noop_free(void* el __attribute__ ((unused)) ) __attribute__ ((const));

char* substring(const int, const char* ) __attribute__ ((pure));

// Chiama la getline e rimuove i caratteri non stampabili
// finali (ad es. \n)
ssize_t my_getline(char **lineptr, size_t *n, FILE *stream)
#ifndef __ICC
  __attribute__ ((nonnull))
#endif
  ;

void print_repetitions(FILE* f, const char c, int rep);

void resource_usage_log(void);

#ifndef MYNDEBUG

#define fail() \
  do {			\
	 fprintf(stderr, "Failing at %s:%d.\n", __FILE__, __LINE__);		\
	 char* pc= NULL;																		\
	 *pc= '\65';																			\
	 exit(1);																				\
  } while (0)

#define my_assert( cond )																\
  do {																						\
	 if (!(cond)) {																		\
		fprintf(stderr, "Assertion " #cond " failed at %s:%d.\n", __FILE__, __LINE__); \
		fail();																				\
	 }																							\
  } while (0)

#else

#define fail() \
  do {			\
	 fprintf(stderr, "Failing at %s:%d.\n", __FILE__, __LINE__);		\
	 char* pc= NULL;																		\
	 *pc= '\65';																			\
	 exit(1);																				\
  } while (0)
/* #define fail() \ */
/*   do {			\ */
/* 	 fprintf(stderr, "Failing at %s:%d.\n", __FILE__, __LINE__);		\ */
/* 	 exit(1);																				\ */
/*   } while (0) */

#define my_assert( cond ) do { } while (0)

#endif

#define assert_lims(lbound, var, ubound) { my_assert(lbound<=(var));	\
	 my_assert((var)<ubound);														\
  } while (0);

#define assert_ulim(var, ubound) {					\
	 my_assert((var)<ubound);							\
  } while (0);



#endif

#define fail_if( cond ) if (cond) {					\
  fail();													\
}

#define fail_if_msg( cond , msg ) if (cond) {	\
	 fprintf(stderr, "%s", msg);						\
	 fail();													\
}

#ifdef __cplusplus
}
#endif
