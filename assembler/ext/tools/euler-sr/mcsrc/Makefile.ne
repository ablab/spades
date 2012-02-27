include common.mak

ifeq '$(CC)' 'cc'
    $(error You must set CC to your g++ compiler, it is now $(CC))
endif

ifndef CC
 $(error You must set CC to your g++ compiler)
endif

ifndef EUSRC
 $(error You must set EUSRC to the root directory where euler was unpacked (probably the current directory))
endif

ifndef PERL5LIB
 $(error You must include the path to "FwgLib.pm" in your PERL5LIB environment variable.  Although this is not required for the build, it is for running helper scripts.  This is $(EUSRC)/fwg_scripts/.)
endif

ifeq (, $(findstring fwg_scripts, $(PERL5LIB)))
 $(error You must include the path to FwgLib.pm in your PERL5LIB environment variable.  Although this is not required for the build, it is for running helper scripts. The path currently contains $(PERL5LIB).)
endif

VALID_ARCH = "x86_64, i686, powerpc"

ifndef MACHTYPE 
 $(error You must set MACHTYPE to describe your architecture, one of $(VALID_ARCH))
endif

ifeq (, $(findstring $(MACHTYPE), $(VALID_ARCH)))
  $(error This may only compile for the following machine types: $(VALID_ARCH), yours is set to "$(MACHTYPE)".  To fix this you may simply need to redefine the MACHTYPE environment variable to that on this list that most closely matches your system.)
endif

DIR_LIST = $(sort $(dir $(wildcard */Makefile)))
DIR_BASE = $(subst /,, $(DIR_LIST))
SUBDIRS  = $(filter-out $(MACHTYPE), $(DIR_BASE))
SUBDIR_LIBS = $(addprefix lib, $(SUBDIRS))
SUBDIR_CLEAN = $(addprefix clean, $(SUBDIRS))


#SUBDIRS = blocks lav

#define MAKE_SUBDIR
#$(1):
#	cd $(1) && $$(MAKE) all
#endef

.PHONY: $(SUBDIRS)

#$(call MAKE_SUBDIR, lav)

all: $(SUBDIRS)


$(foreach subdir, $(SUBDIRS), $(eval $(call MAKE_SUBDIR, $(subdir))))
$(foreach subdir, $(SUBDIRS), $(eval $(call MAKE_DIR_CLEAN, $(subdir))))


clean: $(CLEAN_TARGETS)


