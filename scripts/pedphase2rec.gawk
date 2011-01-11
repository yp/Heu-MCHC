#!/usr/bin/gawk -f
####
#
#
#                               Heu-MCHC
#
# A fast and accurate heuristic algorithm for the haplotype inference
# problem on pedigree data with recombinations and mutations
#
# Copyright (C) 2009,2010,2011  Yuri PIROLA
#
# Distributed under the terms of the GNU General Public License (GPL)
#
#
# This file is part of Heu-MCHC.
#
# Heu-MCHC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Heu-MCHC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
#
####

BEGIN {
	 print "= PEDPHASE EVENTS";
	 lastid=-1;
}


{
	 if (lastid<($2-1)) {
		  lastid=$2-1;
		  pos=0;
		  s=$(NF-1);
		  last=substr(s,1,1);
		  while (s!="") {
				if (last!=substr(s,1,1))
					 printf "= REC %4d %4d  0\n", $2-1, pos;
				last=substr(s,1,1);
				s=substr(s,2);
				pos++;
		  }
		  pos=0;
		  s=$(NF);
		  last=substr(s,1,1);
		  while (s!="" && (substr(s,1,1)=="0" || substr(s,1,1)=="1")) {
				if (last!=substr(s,1,1))
					 printf "= REC %4d %4d  1\n", $2-1, pos;
				last=substr(s,1,1);
				s=substr(s,2);
				pos++;
		  }
	 }
}
