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

include ../common.mak

#CPP_FILES  = $(sort $(wildcard *.cpp))
CPP_FILES  = $(filter-out CreateCompatibility.cpp, $(sort $(wildcard *.cpp)))
LIBOBJS    = $(patsubst %.cpp, %.o, $(CPP_FILES))
ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))


# The following commands configure make to recursively build subdirectories
# so that libraries with names  libsubdir.so are created.


DIR_LIST = $(sort $(dir $(wildcard */Makefile)))
DIR_BASE = $(subst /,, $(DIR_LIST))
SUBDIRS  = $(filter-out $(MACHTYPE), $(DIR_BASE))
SUBDIR_LIBS = $(addprefix lib, $(SUBDIRS))
SUBDIR_CLEAN = $(addprefix clean, $(SUBDIRS))

define MAKE_DIRECTORY
$(1):
	cd $(1) && $$(MAKE) $(MAKECMDGOALS)	
endef

define MAKE_LIBRARY
$(join lib, $(1)):
	cd $(1) && $$(MAKE) $(MAKECMDGOALS)
endef

BUILDLIBS = common $(SUBDIR_LIBS) 


all: compatibility $(BUILDLIBS)
	@ echo compatibility $(BUILDLIBS)

# rules for running make 'libNNN'
$(foreach libdir,$(SUBDIRS),$(eval $(call MAKE_LIBRARY,$(libdir))))
# rules for running make 'NNN'
$(foreach libdir,$(SUBDIRS),$(eval $(call MAKE_DIRECTORY,$(libdir))))
# rules for cleaning
$(foreach libdir,$(SUBDIRS),$(eval $(call MAKE_DIR_CLEAN,$(libdir))))

.PHONY: $(BUILDLIBS) $(SUBDIRS) compatibility


include ../make.rules
include $(CURDIR)/make.dep


common: arch_dir $(MACHTYPE)/libcommon.$(LIB_EXT)


$(MACHTYPE)/libcommon.$(LIB_EXT): $(ARCH_OBJS)
	@ echo mt: $(MACHTYPE)
	@ echo lt $(LIB_TOOL) 
	@ echo lc $(LIB_CREATE)
	$(LIB_TOOL) $(LIB_CREATE) $@ $(ARCH_OBJS)

clean: $(SUBDIR_CLEAN)
	@ echo $(SUBDIR_CLEAN)
	-rm -f $(MACHTYPE)/*.o $(MACHTYPE)/*.$(LIB_EXT)

compatibility: $(MACHTYPE)/compatibility.h

$(MACHTYPE)/compatibility.h: CreateCompatibility.cpp compatibility_template.h
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) CreateCompatibility.cpp -o $(MACHTYPE)/CreateCompatibility
	$(MACHTYPE)/CreateCompatibility < compatibility_template.h > $(MACHTYPE)/compatibility.h
