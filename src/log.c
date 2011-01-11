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
#include "log.h"

#include "util.h"
#include "string.h"

void log_format_str(const char* funz,
						  const char* srcf,
						  char* dst,
						  int max_len1,
						  int max_len2) {
  const int lenf= strlen(funz);
  const int lens= strlen(srcf);
  const int tot= max_len1+max_len2+1;
  int max_lenf;
  if (max_len1>max_len2-lens)
	 max_lenf= max_len1;
  else
	 max_lenf= max_len2-lens;
  strncpy(dst, funz, max_lenf);
  if (lenf>max_lenf)
	 dst[max_lenf-2]= dst[max_lenf-1]= '.';
  int pos= MIN(lenf, max_lenf);
  dst[pos]= '@';
  ++pos;
  if (pos+lens>tot) {
	 strcpy(dst+pos, srcf+lens+pos-tot);
	 dst[pos]= dst[pos+1]= '.';
  } else {
	 strncpy(dst+pos, srcf, max_len1+max_len2+1-pos);
	 pos= MIN(max_len1+max_len2+1, pos+lens);
	 for (int i= pos; i<max_len1+max_len2+1; ++i) {
		dst[i]=' ';
	 }
  }
  dst[tot]= '\0';
}
