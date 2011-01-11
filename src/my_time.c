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
#include "my_time.h"
#include "util.h"
#include <stdio.h>


#define DTYPE unsigned long

static const char* default_timer_name= "generic timer";


struct _mytime {
  const char* timer_name;
  DTYPE interval;
  struct timeval start;
  struct timeval stop;
  bool active;
};



static DTYPE
diff_usec(struct timeval start, struct timeval stop) {
  DTYPE start_usec= ((DTYPE)start.tv_sec*(DTYPE)1000000)+(DTYPE)start.tv_usec;
  DTYPE stop_usec= ((DTYPE)stop.tv_sec*(DTYPE)1000000)+(DTYPE)stop.tv_usec;
  return stop_usec-start_usec;
}

static void
compute_difference(pmytime pt) {
  pt->interval+= diff_usec(pt->start, pt->stop);
}

void
MYTIME_print_interval(FILE* file, pmytime pt) {
  my_assert(pt!=NULL);
  my_assert(file!=NULL);
  DTYPE diff= pt->interval;
  if (diff>=(DTYPE)1000) {
// in secs
	 fprintf(file, "@Timer %s. Time elapsed: ", pt->timer_name);
	 DTYPE min= diff/60000000;
	 diff= diff%60000000;
	 if (min>0)
		fprintf(file, "%lum ", min);
	 fprintf(file, "%.3fs\n", ((double)diff/1000000.0));
  } else {
	 fprintf(file, "@Timer %s. Time elapsed: %lumicrosec", pt->timer_name, diff);
  }
}

void
MYTIME_print_current(FILE* file, pmytime pt) {
  MYTIME_stop(pt);
  MYTIME_print_interval(file, pt);
  MYTIME_start(pt);
}



pmytime
MYTIME_create(void) {
  return MYTIME_create_with_name(default_timer_name);
}

pmytime
MYTIME_create_with_name(const char* timer_name) {
  my_assert(timer_name!=NULL);
  pmytime pt= PALLOC(struct _mytime);
  pt->active= false;
  pt->interval= 0;
  pt->timer_name= timer_name;
  return pt;
}

void
MYTIME_destroy(pmytime pt) 
{
  my_assert(pt!=NULL);
  pfree(pt);
}

void
MYTIME_start(pmytime pt) {
  my_assert(pt!=NULL);
  pt->active= true;
  gettimeofday(&pt->start, NULL);
}

void
MYTIME_reset(pmytime pt)
{
  my_assert(pt!=NULL);
  pt->active= false;
  pt->interval= 0;
}

void
MYTIME_stop(pmytime pt)
{
  my_assert(pt!=NULL);
  gettimeofday(&(pt->stop), NULL);
  pt->active= false;
  compute_difference(pt);
}

const char*
MYTIME_getname(pmytime pt)
{
  my_assert(pt!=NULL);
  return pt->timer_name;
}

DTYPE
MYTIME_getinterval(pmytime pt)
{
  my_assert(pt!=NULL);
  return pt->interval;
}



#undef TOUT
#undef DTYPE
