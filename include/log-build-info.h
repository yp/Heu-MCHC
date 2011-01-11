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
/*
** log-build-info.h
**
** Made by Yuri Pirola
** Login   <yuri@yuri>
**
** Started on  Mon Aug 24 00:38:28 2009 Yuri Pirola
** Last update Mon Aug 24 00:38:28 2009 Yuri Pirola
*/

#ifndef __LOG_BUILD_INFO_H__
#define __LOG_BUILD_INFO_H__

#include "log.h"

#ifndef __SRC_DESC
#define __SRC_DESC "not available"
#endif

#ifndef __BUILD_DESC
#define __BUILD_DESC "not available"
#endif

#ifndef __BUILD_DATETIME
#define __BUILD_DATETIME "unknown"
#endif

#ifndef __BUILD_HOST
#define __BUILD_HOST "unknown"
#endif

#ifndef __COMPILER_VER
#ifdef __VERSION__
#define QUOTEME( x ) #x
#define __COMPILER_VER QUOTEME(__VERSION__)
#else
#define __COMPILER_VER "unknown"
#endif
#endif

#define PRINT_SYSTEM_INFORMATIONS													\
  do {																						\
	 INFO("Program compiled %s @ %s.", __BUILD_DATETIME,						\
			__BUILD_HOST);																	\
	 INFO("Source version >%s<.", __SRC_DESC);									\
	 INFO("Configuration  %s.", __BUILD_DESC);									\
	 INFO("Compiler %s", __COMPILER_VER);											\
  } while (0)


#endif /* !LOG-BUILD-INFO_H_ */
