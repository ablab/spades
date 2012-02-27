############################################################################
# Title:          common.mak
# Author:         Mark Chaisson, Glenn Tesler
# Created:        2007
# Last modified:  11/28/2009
# 
# Copyright (c) 2007-2009 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################


############################################################################
# Debugging options
############################################################################

### Which compiler
CPP = g++
#CPP = g++44


### Compiler warnings
CPPWARN =
# For GNU compiler:
#CPPWARN = -Wall -Wextra -Wconversion -Wold-style-cast
# Additional for GNU compiler on OS X:
#   CPPWARN += -Wshorten-64-to-32

### Preserve copies of files that are rewritten by renaming with suffix
### .v1, .v2, ...  E.g., fixed/reads.fasta.spect is written 3-5 timews
NOCLOBBER =
# NOCLOBBER = -DNOCLOBBER

#ASSERT = -DNDEBUG
ASSERT =

############################################################################
# Common stuff
############################################################################

SHELL=/bin/sh

# CPPOPTS is the options used on $(CPP) when compiling a .cpp source file.
# CPPOPTS_lib is the options used on $(CPP) when linking together
# a bunch of .o files and -lxxx files, rather than compiling source.
# On most platforms, these should be the same, but sparc needs them
# to be different.
# If CPPOPTS_lib == "_default_" at the end of this file, it will
# be replaced by CPPOPTS.
# But if one of the platforms (sparc) sets  CPPOPTS_lib to something else,
# it will stay different.

CPPOPTS_lib = _default_

ifneq (, $(filter $(MAKECMDGOALS), debug))
	CPPOPTS := -g # $(ASSERT) #@CPPOPTS Flag to replace
	ifeq (, $(filter %_debug, $(MACHTYPE)))
		MACHTYPE := $(join ${MACHTYPE},_debug)
	endif
else
	CPPOPTS := -O3 $(ASSERT)
endif

ifneq (, $(filter $(MAKECMDGOALS), parallel))
	CPPOPTS += -g -fopenmp
endif



ifneq (, $(filter $(MAKECMDGOALS), debugclean))
	ifeq (, $(filter %_debug, $(MACHTYPE)))
		MACHTYPE := $(join ${MACHTYPE},_debug)
	endif
endif


SRC_DIR=$(EUSRC)
BIN_DIR=$(SRC_DIR)/$(MACHTYPE)/bin

#MYSQL_INC = -I${MYSQL_INC_DIR}
#MYSQL_LIB = -L${MYSQL_LIB_DIR} -lmysqlclient -lz

#where to look
INCLUDEDIR = -I${SRC_DIR}/lib -I${SRC_DIR}/lib/${MACHTYPE}

ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(OBJS))
#ARCH_OBJS  = $(OBJS)


LIBBASE    = $(SRC_DIR)/lib
LIBPATH    = $(SRC_DIR)/lib/$(MACHTYPE)
LIBDIR     = -L$(SRC_DIR)/lib/$(MACHTYPE)


define MAKE_SUBDIR
$(1):
	cd $(1) && ${MAKE} ${MAKECMDGOALS}
endef

define MAKE_DIR_CLEAN
$(join clean, $(1)):
	cd $(1) && $${MAKE} $${MAKECMDGOALS} clean
CLEAN_TARGETS += $(join clean, $(1))
endef


############################################################################
# Architecture specific options.
############################################################################

ifeq ($(MACHTYPE), LINUX)
	CPPOPTS += -D_3_3_ 
	ARCH_LIBS   = $(LIBS) 
	LIB_EXT = so
	LIB_CREATE = -shared
endif

ifeq ($(MACHTYPE), i686)
	LIB_EXT = a
	LIB_CREATE = rcs
	LIB_TOOL = ar 
endif

ifeq ($(MACHTYPE), i386)
	ARCH_LIBS = $(LIBS)
	LIB_EXT = so
	LIB_CREATE = -shared
endif

ifeq ($(MACHTYPE), DARWIN)
	ARCH_LIBS = $(LIBS)
	LIB_EXT = dylib
	LIB_CREATE = -dynamic
endif

ifeq ($(MACHTYPE), powerpc)
	CPPOPTS += -D_3_3_ -D_POWERPC_
	ARCH_LIBS = $(LIBS)
	LIB_EXT = a
	LIB_CREATE = -static -o
	LIB_TOOL = libtool
endif

ifeq ($(MACHTYPE), x86_64)
  CPPOPTS += -m64
	LIB_EXT = a
	LIB_CREATE = rcs
	LIB_TOOL = ar 
endif

ifeq ($(MACHTYPE), x86_64_debug)
  CPPOPTS += -m64
	LIB_EXT = a
	LIB_CREATE = rcs
	LIB_TOOL = ar 
endif


ifeq ($(MACHTYPE), ALPHA)
	CPPOPTS += -O3 $(ASSERT)
endif

# For debugging
CPPOPTS += $(CPPWARN) $(NOCLOBBER)

ifeq ($(MACHTYPE), sparc)
	CPPOPTS += -m64
	CPPOPTS_lib := $(CPPOPTS)

# Solaris needs help to find the GNU assembler "as" and linker "ld"
# instead of the Solaris assembler and linker.
# Use "g++ -print-prog-name=as" and "g++ -print-prog-name=ld"
# to see which versions are being used.
# /usr/ccs/bin/as and /usr/ccs/bin/ld are the Solaris assembler / linker,
# not the GNU versions.  If it indicates these are in use,
# setting "CPPOPTS += -B/pathname" may be sufficient to fix it.
# If not, see the following sites for additional information:
#	http://gcc.gnu.org/faq.html#gas
#	http://gcc.gnu.org/install/configure.html#with-gnu-as
#		(also on that page, see --with-as  --with-glu-ld   --with-ld)
#	http://gcc.gnu.org/ml/gcc/2004-12/msg00346.html

	CPPOPTS += -B/usr/local/bin

	LIB_EXT = a
	LIB_CREATE = rcs
	LIB_TOOL = ar
endif

ifeq ($(CPPOPTS_lib), _default_)
	CPPOPTS_lib = $(CPPOPTS)
endif

MKDEP       = ${EUSRC}/mkdep.pl
MKFILES     = ${EUSRC}/mkfiles.pl

.PHONY: clean


MACHEXECS = $(addprefix $(MACHTYPE)/, $(NOARCH_EXECS))


############################################################################
# define some library locations
############################################################################

LIBCOMMON = $(SRC_DIR)/lib/$(MACHTYPE)/libcommon.${LIB_EXT} 
LIBLAV    = $(SRC_DIR)/lib/$(MACHTYPE)/liblav.${LIB_EXT}
LIBAXT    = $(SRC_DIR)/lib/$(MACHTYPE)/libaxt.${LIB_EXT}
LIBNET    = $(SRC_DIR)/lib/$(MACHTYPE)/libnet.${LIB_EXT}
LIBALIGN =  $(SRC_DIR)/lib/$(MACHTYPE)/libalign.${LIB_EXT} 
LIBBLOCKDB = $(SRC_DIR)/lib/$(MACHTYPE)/libblockdb.${LIB_EXT}
LIBBLOCKS  = $(SRC_DIR)/lib/$(MACHTYPE)/libblocks.${LIB_EXT}
LIBEMBOSS = $(SRC_DIR)/lib/$(MACHTYPE)/libemboss.${LIB_EXT}
LIBHASH    = $(SRC_DIR)/lib/$(MACHTYPE)/libhash.${LIB_EXT}
LIBPARSEBLAST = $(SRC_DIR)/lib/$(MACHTYPE)/libparseblast.${LIB_EXT}
LIBBBBWT = $(SRC_DIR)/lib/$(MACHTYPE)/libbbbwt.${LIB_EXT}
LIBCHAIN = $(SRC_DIR)/lib/$(MACHTYPE)/libchain.${LIB_EXT}
LIBNET = $(SRC_DIR)/lib/$(MACHTYPE)/libnet.${LIB_EXT}
LIBTREE = $(SRC_DIR)/lib/$(MACHTYPE)/libtree.${LIB_EXT}
LIBASSEMBLE = $(SRC_DIR)/lib/$(MACHTYPE)/libassemble.${LIB_EXT}
LIBREGEXP = $(SRC_DIR)/lib/$(MACHTYPE)/libboost_regex.${LIB_EXT}
