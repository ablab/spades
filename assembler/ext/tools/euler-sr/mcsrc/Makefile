############################################################################
# Title:          Makefile
# Author:         Mark Chaisson, Glenn Tesler
# Created:        2007
# Last modified:  11/26/2009
# 
# Copyright (c) 2007-2009 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################

include common.mak

ifndef EUSRC
 $(error You must set EUSRC to the root directory where euler was unpacked (probably the current directory))
endif


VALID_ARCH = "x86_64, i686, powerpc, x86_64_debug, i686_debug, i386, i386_debug, sparc, sparc_debug"

ifndef MACHTYPE 
 $(error You must set MACHTYPE to describe your computer, one of $(VALID_ARCH). If you have a PC or a Mac running Intel, this should be x86_64)
endif

ifeq (, $(findstring $(MACHTYPE), $(VALID_ARCH)))
  $(error This may only compile for the following machine types: $(VALID_ARCH), yours is set to "$(MACHTYPE)".  To fix this you may simply need to redefine the MACHTYPE environment variable to one on this list that most closely matches your system.)
endif

DIR_LIST = $(sort $(dir $(wildcard */Makefile)))
DIR_BASE = $(subst /,, $(DIR_LIST))
SUBDIRS  = $(filter-out $(MACHTYPE), $(DIR_BASE))
SUBDIR_LIBS = $(addprefix lib, $(SUBDIRS))
SUBDIR_CLEAN = $(addprefix clean, $(SUBDIRS))


.PHONY: $(SUBDIRS)
.PHONY: debug debugclean parallel compatibility

all: compatibility $(SUBDIRS)
debug: all
parallel: all


$(foreach subdir,$(SUBDIRS),$(eval $(call MAKE_SUBDIR,$(subdir))))
$(foreach subdir,$(SUBDIRS),$(eval $(call MAKE_DIR_CLEAN,$(subdir))))


clean: $(CLEAN_TARGETS)
debugclean: $(CLEAN_TARGETS)

compatibility: lib/${MACHTYPE}/compatibility.h

lib/${MACHTYPE}/compatibility.h: lib/CreateCompatibility.cpp lib/compatibility_template.h
	@ if [ ! -d $(CURDIR)/lib/$(MACHTYPE) ]; then  mkdir $(CURDIR)/lib/$(MACHTYPE); fi
	@ touch $(CURDIR)/lib/$(MACHTYPE)/.exists
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) lib/CreateCompatibility.cpp -o lib/$(MACHTYPE)/CreateCompatibility
	lib/$(MACHTYPE)/CreateCompatibility < lib/compatibility_template.h > lib/$(MACHTYPE)/compatibility.h
