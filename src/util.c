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
#include "util.h"
#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/resource.h>
#include <unistd.h>

#ifndef NDEBUG
static unsigned long int memsize= 0;
static unsigned long int prev_memsize= 0;
#endif

void* palloc(size_t size)
{
#ifndef NDEBUG
  memsize+= size;
#endif

  void* p= malloc(size);
  if (p==NULL) {
	 fprintf(stderr, "Allocation memory error. Trying to allocate %zu bytes.\n", size);
	 fail();
  }
#ifndef NDEBUG
  if ((memsize-prev_memsize)>>20>=256) {
	 WARN("Allocated %lu MB", memsize>>20);
	 prev_memsize= memsize;
  }
  if ((memsize>>20)>4096) {
	 WARN("Allocated %lu MB", memsize>>20);
	 FATAL("Currently allocated more than 2GB. Abort.");
	 fail();
  }
#endif
  return p;
}

char* c_palloc(size_t size)
{
  return (char*)palloc(size*sizeof(char));
}

void pfree(void* p)
{
  if (p==NULL) {
	 fprintf(stderr, "Cannot free a NULL pointer.\n");
	 fail();
  }
  free(p);
}

char* alloc_and_copy(char* source)
{
  size_t len= strlen(source);
  char* ris= c_palloc(len+1);
  strncpy(ris, source, len+1);
  return ris;
}

void noop_free(void* el)
{}

char* substring(const int index, const char* string){
  my_assert(index>=0);
  my_assert(string != NULL);
  size_t slen= strlen(string);
  char* ris = c_palloc(slen-index+1);
  ris = (char*)memcpy(ris, (void*)(string+index), slen-index);
  ris[slen-index] = '\0';
  return ris;
}//end subString

ssize_t my_getline(char **lineptr, size_t *n, FILE *stream) {
  ssize_t ris= getline(lineptr, n, stream);
  while (ris>0 && (*lineptr)[ris-1]<' ') {
	 --ris;
	 (*lineptr)[ris]= '\0';
  }
  return ris;
}

void print_repetitions(FILE* f, const char c, int rep) {
  char* s= c_palloc(rep+1);
  for (int i= 0; i<rep; ++i)
	 s[i]= c;
  s[rep]='\0';
  fprintf(f, "%s", s);
  pfree(s);
}

void resource_usage_log(void) {
  struct rusage *rusage= PALLOC(struct rusage);
  if (getrusage(RUSAGE_SELF, rusage)==0) {
	 STATS("user time   %10ld s %7ld millisec.", rusage->ru_utime.tv_sec, rusage->ru_utime.tv_usec/1000);
	 STATS("system time %10ld s %7ld millisec.", rusage->ru_stime.tv_sec, rusage->ru_stime.tv_usec/1000);
	 char buf[1000];
	 snprintf(buf, 1000, "/proc/%u/statm", (unsigned)getpid());
	 FILE* pf = fopen(buf, "r");
	 if (pf) {
		unsigned int size; //       total program size
		fscanf(pf, "%u", &size);
		STATS("used memory %10u KB", size);
	 }
	 fclose(pf);
  }
  pfree(rusage);
}

double rnd_01(void) {
  return ((double)rand())/RAND_MAX;
}

unsigned int
rnd_in_range(unsigned int min, unsigned int max) {
  return min+(unsigned int)(rnd_01()*(max-min));
}
