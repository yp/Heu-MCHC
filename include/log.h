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
 * @file log.h
 *
 * Funzioni per la stampa di messaggi di log.
 *
 **/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _LOG_H_
#define _LOG_H_


#define LOG_LEVEL_FATAL (0)
#define LOG_LEVEL_ERROR (1)
#define LOG_LEVEL_WARN (2)
#define LOG_LEVEL_INFO (3)
#define LOG_LEVEL_STATS (3)
#define LOG_LEVEL_DEBUG (4)
#define LOG_LEVEL_TRACE (5)
#define LOG_LEVEL_FINETRACE (6)

#ifndef LOG_THRESHOLD
#define LOG_THRESHOLD LOG_LEVEL_INFO
#endif

#ifndef LOG_PREFIX
#define LOG_PREFIX "* "
#endif

#include <stddef.h>

  void log_format_str(const char* funz, const char* srcf, char* dst, int max_len1, int max_len2);


#endif // _LOG_H_


#ifdef LOG
#undef __INTERNAL_LOG
#undef LOG
#undef LOG_NOTN
#undef FATAL
#undef ERROR
#undef WARN
#undef INFO
#undef STATS
#undef DEBUG
#undef TRACE
#undef FINETRACE
#undef FATAL_NOTN
#undef ERROR_NOTN
#undef WARN_NOTN
#undef INFO_NOTN
#undef STATS_NOTN
#undef DEBUG_NOTN
#undef TRACE_NOTN
#undef FINETRACE_NOTN
#undef FATAL_IF_NOTN
#undef ERROR_IF_NOTN
#undef WARN_IF_NOTN
#undef INFO_IF_NOTN
#undef STATS_IF_NOTN
#undef DEBUG_IF_NOTN
#undef TRACE_IF_NOTN
#undef FINETRACE_IF_NOTN
#undef LOG_FATAL_ENABLED
#undef LOG_ERROR_ENABLED
#undef LOG_WARN_ENABLED
#undef LOG_INFO_ENABLED
#undef LOG_STATS_ENABLED
#undef LOG_DEBUG_ENABLED
#undef LOG_TRACE_ENABLED
#undef LOG_FINETRACE_ENABLED
#endif

#ifdef LOG_MSG

#include <stdio.h>
#include <string.h>

#define MAX_LEN_FUNC_NAME 16
#define MAX_LEN_FILE_NAME 12


#define BASIC_LOG(fmt, ...) fprintf(stderr, fmt, ## __VA_ARGS__)

#define __BASIC_STATS(format, ...) do {											\
	 BASIC_LOG("@ STATS "																\
				  format "%s",																\
				  ## __VA_ARGS__);														\
  } while (0)

#define LOG_NOTN(level, prefix, ...) __INTERNAL_LOG(level, prefix, ## __VA_ARGS__, "")
#define LOG(level, prefix, ...) __INTERNAL_LOG(level, prefix, ## __VA_ARGS__, "\n")

#define __INTERNAL_LOG(level, prefix, format, ...) do {		\
	 if (level<=LOG_THRESHOLD) {														\
		char __my_internal_buff__[MAX_LEN_FUNC_NAME+MAX_LEN_FILE_NAME+2];	\
		log_format_str(__func__, __FILE__, __my_internal_buff__,				\
							MAX_LEN_FUNC_NAME, MAX_LEN_FILE_NAME);					\
		BASIC_LOG(LOG_PREFIX prefix "(%s:%4d) "							\
					format "%s",															\
					__my_internal_buff__, __LINE__,									\
					## __VA_ARGS__);														\
	 }																							\
  } while (0)

#if (LOG_LEVEL_FATAL<=LOG_THRESHOLD)
#define LOG_FATAL_ENABLED
#endif

#if (LOG_LEVEL_ERROR<=LOG_THRESHOLD)
#define LOG_ERROR_ENABLED
#endif

#if (LOG_LEVEL_WARN<=LOG_THRESHOLD)
#define LOG_WARN_ENABLED
#endif

#if (LOG_LEVEL_INFO<=LOG_THRESHOLD)
#define LOG_INFO_ENABLED
#endif

#if (LOG_LEVEL_STATS<=LOG_THRESHOLD)
#define LOG_STATS_ENABLED
#endif

#if (LOG_LEVEL_DEBUG<=LOG_THRESHOLD)
#define LOG_DEBUG_ENABLED
#endif

#if (LOG_LEVEL_TRACE<=LOG_THRESHOLD)
#define LOG_TRACE_ENABLED
#endif

#if (LOG_LEVEL_FINETRACE<=LOG_THRESHOLD)
#define LOG_FINETRACE_ENABLED
#endif

#else

#define LOG(level, prefix, ...) do { } while (0)
#define LOG_NOTN(level, prefix, ...) do { } while (0)

#endif

#ifdef LOG_FATAL_ENABLED
#define FATAL(...) LOG(LOG_LEVEL_FATAL, "FATAL", ## __VA_ARGS__)
#define FATAL_MSG(fmt, ...) BASIC_LOG(fmt, ## __VA_ARGS__)
#define FATAL_NOTN(...) LOG_NOTN(LOG_LEVEL_FATAL, "FATAL", ## __VA_ARGS__)
#define FATAL_IF(cond, ...) if ((cond)) {			\
  FATAL(##__VA_ARGS__);									\
}
#define FATAL_IF_NOTN(cond, ...) if ((cond)) {			\
  FATAL_NOTN(##__VA_ARGS__);									\
}
#else
#define FATAL(...) do {} while (0)
#define FATAL_MSG(...) do {} while (0)
#define FATAL_IF(...) do {} while (0)
#define FATAL_NOTN(...) do {} while (0)
#define FATAL_IF_NOTN(...) do {} while (0)
#endif

#ifdef LOG_ERROR_ENABLED
#define ERROR(...) LOG(LOG_LEVEL_ERROR, "ERROR", ## __VA_ARGS__)
#define ERROR_MSG(fmt, ...) BASIC_LOG(fmt,## __VA_ARGS__)
#define ERROR_NOTN(...) LOG_NOTN(LOG_LEVEL_ERROR, "ERROR", ## __VA_ARGS__)
#define ERROR_IF(cond, ...) if ((cond)) {			\
	 ERROR( __VA_ARGS__);								\
}
#define ERROR_IF_NOTN(cond, ...) if ((cond)) {			\
	 ERROR_NOTN( __VA_ARGS__);								\
}
#else
#define ERROR(...) do {} while (0)
#define ERROR_MSG(...) do {} while (0)
#define ERROR_IF(...) do {} while (0)
#define ERROR_NOTN(...) do {} while (0)
#define ERROR_IF_NOTN(...) do {} while (0)
#endif

#ifdef LOG_WARN_ENABLED
#define WARN(...) LOG(LOG_LEVEL_WARN, "WARN ", ## __VA_ARGS__)
#define WARN_MSG(fmt, ...) BASIC_LOG(fmt,## __VA_ARGS__)
#define WARN_NOTN(...) LOG_NOTN(LOG_LEVEL_WARN, "WARN ", ## __VA_ARGS__)
#define WARN_IF(cond, ...) if ((cond)) {			\
	 WARN( __VA_ARGS__);									\
}
#define WARN_IF_NOTN(cond, ...) if ((cond)) {			\
	 WARN_NOTN( __VA_ARGS__);									\
}
#else
#define WARN(...) do {} while (0)
#define WARN_MSG(...) do {} while (0)
#define WARN_IF(...) do {} while (0)
#define WARN_NOTN(...) do {} while (0)
#define WARN_IF_NOTN(...) do {} while (0)
#endif

#ifdef LOG_INFO_ENABLED
#define INFO(...) LOG(LOG_LEVEL_INFO, "INFO ", ## __VA_ARGS__)
#define INFO_MSG(fmt, ...) BASIC_LOG(fmt,## __VA_ARGS__)
#define INFO_NOTN(...) LOG_NOTN(LOG_LEVEL_INFO, "INFO ", ## __VA_ARGS__)
#define INFO_IF(cond, ...) if ((cond)) {			\
	 INFO( __VA_ARGS__);									\
}
#define INFO_IF_NOTN(cond, ...) if ((cond)) {			\
	 INFO_NOTN( __VA_ARGS__);									\
}
#else
#define INFO(...) do {} while (0)
#define INFO_MSG(...) do {} while (0)
#define INFO_IF(...) do {} while (0)
#define INFO_NOTN(...) do {} while (0)
#define INFO_IF_NOTN(...) do {} while (0)
#endif

#ifdef LOG_STATS_ENABLED
#define STATS(fmt, ...) __BASIC_STATS(fmt,##__VA_ARGS__, "\n")
#define STATS_MSG(fmt, ...) BASIC_LOG(fmt,##__VA_ARGS__)
#define STATS_NOTN(...) __BASIC_STATS(__VA_ARGS__, "")
#define STATS_IF(cond, ...) if ((cond)) {				\
	 STATS( __VA_ARGS__);									\
  }
#define STATS_IF_NOTN(cond, ...) if ((cond)) {			\
	 STATS_NOTN( __VA_ARGS__);									\
  }
#else
#define STATS(...) do {} while (0)
#define STATS_MSG(...) do {} while (0)
#define STATS_IF(...) do {} while (0)
#define STATS_NOTN(...) do {} while (0)
#define STATS_IF_NOTN(...) do {} while (0)
#endif

#ifdef LOG_DEBUG_ENABLED
#define DEBUG(...) LOG(LOG_LEVEL_DEBUG, "DEBUG", ## __VA_ARGS__)
#define DEBUG_MSG(fmt, ...) BASIC_LOG(fmt,## __VA_ARGS__)
#define DEBUG_NOTN(...) LOG_NOTN(LOG_LEVEL_DEBUG, "DEBUG", ## __VA_ARGS__)
#define DEBUG_IF(cond, ...) if ((cond)) {	\
	 DEBUG( __VA_ARGS__);						\
}
#define DEBUG_IF_NOTN(cond, ...) if ((cond)) {	\
	 DEBUG_NOTN( __VA_ARGS__);						\
}
#else
#define DEBUG(...) do {} while (0)
#define DEBUG_MSG(...) do {} while (0)
#define DEBUG_IF(...) do {} while (0)
#define DEBUG_NOTN(...) do {} while (0)
#define DEBUG_IF_NOTN(...) do {} while (0)
#endif

#ifdef LOG_TRACE_ENABLED
#define TRACE(...) LOG(LOG_LEVEL_TRACE, "TRACE", ## __VA_ARGS__)
#define TRACE_MSG(fmt, ...) BASIC_LOG(fmt,## __VA_ARGS__)
#define TRACE_NOTN(...) LOG_NOTN(LOG_LEVEL_TRACE, "TRACE", ## __VA_ARGS__)
#define TRACE_IF(cond, ...) if ((cond)) {			\
	 TRACE( __VA_ARGS__);								\
}
#define TRACE_IF_NOTN(cond, ...) if ((cond)) {			\
	 TRACE_NOTN( __VA_ARGS__);								\
}
#else
#define TRACE(...) do {} while (0)
#define TRACE_MSG(...) do {} while (0)
#define TRACE_IF(...) do {} while (0)
#define TRACE_NOTN(...) do {} while (0)
#define TRACE_IF_NOTN(...) do {} while (0)
#endif

#ifdef LOG_FINETRACE_ENABLED
#define FINETRACE(...) LOG(LOG_LEVEL_FINETRACE, "TRAC+", ## __VA_ARGS__)
#define FINETRACE_MSG(fmt, ...) BASIC_LOG(fmt,## __VA_ARGS__)
#define FINETRACE_NOTN(...) LOG_NOTN(LOG_LEVEL_FINETRACE, "TRAC+", ## __VA_ARGS__)
#define FINETRACE_IF(cond, ...) if ((cond)) {			\
	 FINETRACE( __VA_ARGS__);								\
}
#define FINETRACE_IF_NOTN(cond, ...) if ((cond)) {			\
	 FINETRACE_NOTN( __VA_ARGS__);								\
}
#else
#define FINETRACE(...) do {} while (0)
#define FINETRACE_MSG(...) do {} while (0)
#define FINETRACE_IF(...) do {} while (0)
#define FINETRACE_NOTN(...) do {} while (0)
#define FINETRACE_IF_NOTN(...) do {} while (0)
#endif


#ifdef __cplusplus
}
#endif
