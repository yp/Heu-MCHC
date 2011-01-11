.PHONY: all reall clean build ped-hi reped-hi ped-gen reped-gen gt-ped-gen regt-ped-gen to-hap reto-hap conv reconv gen-ped-gen regen-ped-gen
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

DEFAULT_COMPILER=gcc
DEFAULT_MATH_LIB=none

ifndef STATUS
STATUS=debug
endif

ifndef LINKING
LINKING=dynamic
endif

ifneq ($(LINKING), dynamic)
ifneq ($(LINKING), static)
LINKING=dynamic
endif
endif


ifneq ($(COMPILER), icc)
ifneq ($(COMPILER), gcc)

ifndef DEFAULT_COMPILER
$(error No compiler selected and no default compiler configured.)
endif

$(warning No compiler selected. Using the default compiler ${DEFAULT_COMPILER}.)
COMPILER=$(DEFAULT_COMPILER)

endif
endif

$(info Including definitions for the compiler >${COMPILER}<)
include Makefile.$(COMPILER)

ifneq ($(MATH_LIB), mkl)
ifneq ($(MATH_LIB), framewave)
ifneq ($(MATH_LIB), none)

ifndef DEFAULT_MATH_LIB
$(error No default mathematical library defined. Possible values: mkl, framewave, none)
endif
$(warning No mathematical library selected. Using default >${DEFAULT_MATH_LIB}<)
override MATH_LIB:=$(DEFAULT_MATH_LIB)

endif
endif
endif


$(info Including definitions for the mathematical library >${MATH_LIB}<)
include Makefile.math_lib_$(MATH_LIB)


# Types of events: recombinations and/or mutations
ifneq ($(RECOMBINATIONS), yes)
ifneq ($(RECOMBINATIONS), no)
RECOMBINATIONS=yes
endif
endif

ifneq ($(MUTATIONS), yes)
ifneq ($(MUTATIONS), no)
MUTATIONS=yes
endif
endif

DEVENTS=
ifeq ($(MUTATIONS), no)
ifeq ($(RECOMBINATIONS), no)
$(error "Both recombinations and mutations are excluded.")
endif
DEVENTS=-DNO_MUTATIONS
endif

ifeq ($(RECOMBINATIONS), no)
DEVENTS=-DNO_RECOMBINATIONS
endif

# Simplification of cycle constraints via Gauss
ifneq ($(SIMPLIFICATION), yes)
ifneq ($(SIMPLIFICATION), no)
$(warning Assuming simplification of cycle-constraints.)
SIMPLIFICATION=yes
endif
endif

ifeq ($(SIMPLIFICATION), yes)
$(info Simplification of cycle-constraints enabled.)
DSIMPL=-D_SIMPLIFY_CYCLE_CONSTRAINTS
else
$(info Simplification of cycle-constraints disabled.)
DSIMPL=-U_SIMPLIFY_CYCLE_CONSTRAINTS
endif



COMPILING_DESC=$(COMPILER)-$(CC)-$(CXX)-math_lib_$(MATH_LIB_DESC)-$(LINKING)-rec_$(RECOMBINATIONS)-mut_$(MUTATIONS)-simpl_$(SIMPLIFICATION)-$(STATUS)
BASE_BIN_DIR= bin
BIN_DIR= $(BASE_BIN_DIR)/_tmp/$(COMPILING_DESC)
SRC_DIR= src
BASE_OBJ_DIR= obj
OBJ_DIR= $(BASE_OBJ_DIR)/$(COMPILING_DESC)
TEST_DIR= test
OBJT_DIR= test-obj/$(COMPILING_DESC)
BINT_DIR= test-bin
INCLUDE_DIR= include
ALL_DIR= $(BIN_DIR) $(SRC_DIR) $(OBJ_DIR)

#####################
# Compiler options for production code
#
ifeq ($(STATUS), production)
PART_DFLAGS=-DLOG_MSG -DLOG_THRESHOLD=LOG_LEVEL_INFO -DNDEBUG -DMYNDEBUG
endif
#####################

#####################
# Compiler options for verbose code
#
ifeq ($(STATUS), verbose)
PART_DFLAGS=-DLOG_MSG -DLOG_THRESHOLD=LOG_LEVEL_DEBUG -DNDEBUG
endif
#####################

#####################
# Compiler options for debugging code
#
ifeq ($(STATUS), debug)
PART_DFLAGS=-DLOG_MSG -DLOG_THRESHOLD=LOG_LEVEL_DEBUG
endif
#####################

COMPFLAGS=$(OPTP) $(OPTD) $(OPTPROF)

INCLUDE=-I. -I$(INCLUDE_DIR)/ $(INCLUDE_MATH) -Ipers-lib/local/include/

ifeq ($(LINKING), static)
LIBS_M4RI=pers-lib/local/lib/libm4ri.a
#LFLAGS=-static
else
LIBS_M4RI=-lm4ri
LFLAGS=
endif

LIBS=-Lpers-lib/local/lib $(LIBS_M4RI) $(LIBS_MATH) -lm $(LFLAGS)

#####################
# System information
__DATETIME=`LANG=C date`
__HOST=`hostname`
__SRC_DESC=`git branch -v | grep '^\*'`
__COMPILER_VER=`$(CC) --version | head -n 1`
DSYSINFO=-D__BUILD_DATETIME="\"$(__DATETIME)\"" -D__BUILD_HOST="\"$(__HOST)\""
DSYSINFO+=-D__SRC_DESC="\"$(__SRC_DESC)\""
DSYSINFO+=-D__BUILD_DESC="\"compiler=$(COMPILER) math_lib=$(MATH_LIB) status=$(STATUS) recombinations=$(RECOMBINATIONS) mutations=$(MUTATIONS) cycle_constraints_simplification=$(CCSIMPL) linking=$(LINKING) compflags='$(COMPFLAGS)' cflags='$(OPTBCC)' cxxflags='$(OPTBCXX)' \""
DSYSINFO+=-D__COMPILER_VER="\"$(__COMPILER_VER)\""
#####################

#####################
# Simboli del preprocessore
#
# NDEBUG Disabilita vari controlli
# LOG_GRAPHS Abilita l'output dei locus graph e del constrint graph
# LOG_MSG Abilita il log dei messaggi
# LOG_THRESHOLD Livello di visualizzazione dei messaggi di log.
# Puo' assumere il valore LOG_LEVEL_XXX con XXX uguale a 
# FATAL, ERROR, WARN, INFO, DEBUG, TRACE
# Assume il valore di LOG_LEVEL_INFO se non definito.
# NO_RECOMBINATIONS does not considers recombinations as permitted events
# NO_MUTATIONS does not considers mutations as permitted events
#
DFLAGS=-D_GNU_SOURCE $(DEVENTS) $(DSIMPL) $(PART_DFLAGS) $(DFLAGSMATH) $(DSYSINFO) $(EXT_DEFINE)
#####################

CFLAGS=$(OPTBCC) $(COMPFLAGS) $(DFLAGS) $(INCLUDE)
CXXFLAGS=$(OPTBCXX) $(COMPFLAGS) $(DFLAGS) $(INCLUDE)

include_FILES= \
        $(INCLUDE_DIR)/belief-propagation.hpp \
        $(INCLUDE_DIR)/bp_double.h \
        $(INCLUDE_DIR)/bp_gmp.h \
        $(INCLUDE_DIR)/gen-ped-IO.h \
        $(INCLUDE_DIR)/irregular_mat.hpp \
        $(INCLUDE_DIR)/locus-graph.hpp \
        $(INCLUDE_DIR)/log.h \
        $(INCLUDE_DIR)/my_time.h \
        $(INCLUDE_DIR)/util.h \
        $(INCLUDE_DIR)/solve.hpp \
        $(INCLUDE_DIR)/data.hpp \


base_SOURCE= \
	$(SRC_DIR)/my_time.c \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/gen-ped-IO.c \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/data.cpp \

base_OBJ= \
	$(OBJ_DIR)/my_time.o \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/gen-ped-IO.o \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/data.o \

ped_hi_SOURCE= \
	$(SRC_DIR)/locus-graph.cpp \
	$(SRC_DIR)/belief-propagation.cpp \
	$(SRC_DIR)/solve.cpp \
	$(SRC_DIR)/ped-hi-main.cpp

ped_hi_OBJ= \
	$(OBJ_DIR)/locus-graph.o \
	$(OBJ_DIR)/belief-propagation.o \
	$(OBJ_DIR)/solve.o \
	$(OBJ_DIR)/ped-hi-main.o

ped_hi_PROG=$(BIN_DIR)/ped-hi

ped_gen_SOURCE= \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/ped-generation.cpp

ped_gen_OBJ= \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/ped-generation.o

ped_gen_PROG=$(BIN_DIR)/gen-ped

gt_ped_gen_SOURCE= \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/gt-ped-generation.cpp

gt_ped_gen_OBJ= \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/gt-ped-gen-opt.o \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/gt-ped-generation.o

gt_ped_gen_PROG=$(BIN_DIR)/gen-gt-ped


to_hap_SOURCE= \
	$(SRC_DIR)/hi-options.c \
	$(SRC_DIR)/hi.cpp

to_hap_OBJ= \
	$(OBJ_DIR)/hi-options.o \
	$(OBJ_DIR)/hi.o

to_hap_PROG=$(BIN_DIR)/to-hap


min_sol_SOURCE= \
	$(SRC_DIR)/hi-options.c \
	$(SRC_DIR)/minimal-solution.cpp

min_sol_OBJ= \
	$(OBJ_DIR)/hi-options.o \
	$(OBJ_DIR)/minimal-solution.o

min_sol_PROG=$(BIN_DIR)/min-sol


conv_SOURCE= \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/gen-ped-IO.c \
	$(SRC_DIR)/conversion-2pedphase.c

conv_OBJ= \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/gen-ped-IO.o \
	$(OBJ_DIR)/conversion-2pedphase.o

conv_PROG=$(BIN_DIR)/ped2pedphase


simconv_SOURCE= \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/gen-ped-IO.c \
	$(SRC_DIR)/conversion-2simwalk.c

simconv_OBJ= \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/gen-ped-IO.o \
	$(OBJ_DIR)/conversion-2simwalk.o

simconv_PROG=$(BIN_DIR)/ped2simwalk


bobconv_SOURCE= \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/gen-ped-IO.c \
	$(SRC_DIR)/conversion-2bob.c

bobconv_OBJ= \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/gen-ped-IO.o \
	$(OBJ_DIR)/conversion-2bob.o

bobconv_PROG=$(BIN_DIR)/ped2bob


gen_ped_gen_SOURCE= \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/gen-ped-generation.cpp

gen_ped_gen_OBJ= \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/gt-ped-gen-opt.o \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/gen-ped-generation.o

gen_ped_gen_PROG=$(BIN_DIR)/gen-gen-ped


ped_only_gen_SOURCE= \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/ped-only-generation.cpp

ped_only_gen_OBJ= \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/ped-only-generation-opt.o \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/ped-only-generation.o

ped_only_gen_PROG=$(BIN_DIR)/ped-gen


all	: build
	@echo '     All compiled!'

.make   :  Makefile Makefile.*
	@rm -rf $(BASE_OBJ_DIR)/* $(BASE_BIN_DIR)/*; \
	echo " !   Makefile modified. Cleaning enforced."; \
        touch .make

$(SRC_DIR)/hi-options.c: $(SRC_DIR)/hi-options.ggo
	@echo " -   Generating the configuration parser..."; \
	gengetopt --output-dir=`dirname $@` \
	--file-name=`basename $@ .c` < $<; \
	rm -f $(INCLUDE_DIR)/`basename $@ .c`.h; \
	mv $(SRC_DIR)/`basename $@ .c`.h $(INCLUDE_DIR)/

$(SRC_DIR)/gt-ped-gen-opt.c: $(SRC_DIR)/gt-ped-gen-opt.ggo
	@echo " -   Generating the configuration parser..."; \
	gengetopt --output-dir=`dirname $@` \
	--file-name=`basename $@ .c` < $<; \
	rm -f $(INCLUDE_DIR)/`basename $@ .c`.h; \
	mv $(SRC_DIR)/`basename $@ .c`.h $(INCLUDE_DIR)/

$(SRC_DIR)/ped-only-generation-opt.c: $(SRC_DIR)/ped-only-generation-opt.ggo
	@echo " -   Generating the configuration parser..."; \
	gengetopt --output-dir=`dirname $@` \
	--file-name=`basename $@ .c` < $<; \
	rm -f $(INCLUDE_DIR)/`basename $@ .c`.h; \
	mv $(SRC_DIR)/`basename $@ .c`.h $(INCLUDE_DIR)/

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(include_FILES)
	@echo ' _   Compiling   $<'; \
	mkdir -pv $(dir $@) > /dev/null; \
	$(CC) $(CFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(include_FILES)
	@echo ' _   Compiling   $<'; \
	mkdir -pv $(dir $@) > /dev/null; \
	$(CXX) $(CXXFLAGS) -o $@ -c $<


build	: .make ped-hi ped-gen to-hap min-sol gt-ped-gen conv simconv bobconv gen-ped-gen ped-only-gen
	@echo "     Config. C:    ${CFLAGS}"; \
	echo "     Config. CXX:  ${CXXFLAGS}"; \
	echo "     Compiler CC:  ${CC}"; \
	echo "     Compiler CXX: ${CXX}"; \
	echo "     Symbols:      ${DFLAGS}"; \
	echo "     Status:       ${STATUS}"; \
	echo "     Math lib:     ${MATH_LIB_DESC}"; \
	echo '     All built!'

ped-hi	: $(ped_hi_PROG)
	@ln -f $(ped_hi_PROG) $(BASE_BIN_DIR)

reped-hi	: clean ped-hi
	@echo 'Cleaned and rebuilt!'

$(ped_hi_PROG)	: $(base_OBJ) $(ped_hi_OBJ)
	@echo ' *   Linking    ' `basename $(ped_hi_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	find $(BIN_DIR) -name $(notdir $@) -delete ; \
	$(CXX) $^ $(LIBS) -o $@ $(CXXFLAGS) ; \


ped-gen	: $(ped_gen_PROG)
	@ln -f $(ped_gen_PROG) $(BASE_BIN_DIR)

reped-gen	: clean ped-gen
	@echo '     Cleaned and rebuilt!'

$(ped_gen_PROG)	: $(ped_gen_OBJ)
	@echo ' *   Linking    ' `basename $(ped_gen_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ -o $@ $(CXXFLAGS) ; \


gt-ped-gen	: $(gt_ped_gen_PROG)
	@ln -f $(gt_ped_gen_PROG) $(BASE_BIN_DIR)

regt-ped-gen	: clean gt-ped-gen
	@echo '      Cleaned and rebuilt!'

$(gt_ped_gen_PROG)	: $(gt_ped_gen_OBJ)
	@echo ' *   Linking    ' `basename $(gt_ped_gen_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ -o $@ $(CXXFLAGS) ; \


to-hap	: $(to_hap_PROG)
	@ln -f $(to_hap_PROG) $(BASE_BIN_DIR)

reto-hap	: clean to-hap
	@echo '      Cleaned and rebuilt!'

$(to_hap_PROG)	: $(base_OBJ) $(to_hap_OBJ)
	@echo ' *   Linking    ' `basename $(to_hap_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ $(LIBS) -o $@ $(CXXFLAGS) ; \


min-sol	: $(min_sol_PROG)
	@ln -f $(min_sol_PROG) $(BASE_BIN_DIR)

remin-sol	: clean min-sol
	@echo '      Cleaned and rebuilt!'

$(min_sol_PROG)	: $(base_OBJ) $(min_sol_OBJ)
	@echo ' *   Linking    ' `basename $(min_sol_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ $(LIBS) -o $@ $(CXXFLAGS) ; \


conv	: $(conv_PROG)
	@ln -f $(conv_PROG) $(BASE_BIN_DIR)

reconv	: clean conv
	@echo '      Cleaned and rebuilt!'

$(conv_PROG)	: $(conv_OBJ)
	@echo ' *   Linking    ' `basename $(conv_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ -o $@ $(CXXFLAGS) ; \


simconv	: $(simconv_PROG)
	@ln -f $(simconv_PROG) $(BASE_BIN_DIR)

resimconv	: clean simconv
	@echo '      Cleaned and rebuilt!'

$(simconv_PROG)	: $(simconv_OBJ)
	@echo ' *   Linking    ' `basename $(simconv_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ -o $@ $(CXXFLAGS) ; \


bobconv	: $(bobconv_PROG)
	@ln -f $(bobconv_PROG) $(BASE_BIN_DIR)

rebobconv	: clean bobconv
	@echo '      Cleaned and rebuilt!'

$(bobconv_PROG)	: $(bobconv_OBJ)
	@echo ' *   Linking    ' `basename $(bobconv_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ -o $@ $(CXXFLAGS) ; \


gen-ped-gen	: $(gen_ped_gen_PROG)
	@ln -f $(gen_ped_gen_PROG) $(BASE_BIN_DIR)

regen-ped-gen	: clean gen-ped-gen
	@echo '     Cleaned and rebuilt!'

$(gen_ped_gen_PROG)	: $(gen_ped_gen_OBJ)
	@echo ' *   Linking    ' `basename $(gen_ped_gen_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ -o $@ $(CXXFLAGS) ; \

ped-only-gen	: $(ped_only_gen_PROG)
	@ln -f $(ped_only_gen_PROG) $(BASE_BIN_DIR)

reped-only-gen	: clean ped-only-gen
	@echo '     Cleaned and rebuilt!'

$(ped_only_gen_PROG)	: $(ped_only_gen_OBJ)
	@echo ' *   Linking    ' `basename $(ped_only_gen_PROG)`; \
	mkdir -pv $(BIN_DIR) > /dev/null; \
	$(CXX) $^ -o $@ $(CXXFLAGS) ; \


reall	: clean all
	@echo '     Cleaned and rebuilt all!'

clean 	:
	@echo '     Cleaning objects and programs' ; \
	rm -rf $(BASE_OBJ_DIR)/* $(BASE_BIN_DIR)/*; \

