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
#ifndef _MYTIME_H_
#define _MYTIME_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>

#include "log.h"

#define MYTIME_LOG( timer ) \
  do {														\
	 STATS("timer %-22s %15lu microsec",			\
			 MYTIME_getname(timer),						\
			 MYTIME_getinterval(timer));				\
  } while (0)


typedef struct _mytime* pmytime;

void
MYTIME_print_interval(FILE*, pmytime);

void
MYTIME_print_current(FILE*, pmytime);

pmytime
MYTIME_create(void);

pmytime
MYTIME_create_with_name(const char*);

void
MYTIME_destroy(pmytime pt);

void
MYTIME_start(pmytime pt);

void
MYTIME_reset(pmytime pt);

void
MYTIME_stop(pmytime pt);

const char*
MYTIME_getname(pmytime pt);

unsigned long
MYTIME_getinterval(pmytime pt);

#ifdef __cplusplus
}
#endif

#endif
