include ../common.mak
#CPP += -p 
LIBOBJS = ${MACHTYPE}/DeBruijnGraph.o \
	 ${MACHTYPE}/ReadIntervals.o \
	 ${MACHTYPE}/IntervalGraph.o \
	 ${MACHTYPE}/GraphReader.o \
	 ${MACHTYPE}/Spectrum.o \
	 ${MACHTYPE}/ReadPaths.o \
	 ${MACHTYPE}/ReadMap.o \
	 ${MACHTYPE}/ReadPos.o \
	 ${MACHTYPE}/SortedTupleList.o \
	 ${MACHTYPE}/MapContigs.o \
	 ${MACHTYPE}/ContigMap.o \
	 ${MACHTYPE}/Tuple.o \
	 ${MACHTYPE}/ListSpectrum.o \
	 ${MACHTYPE}/StringTuple.o \
	 ${MACHTYPE}/HashedSpectrum.o \
	 ${MACHTYPE}/BitSpectrum.o \
	 ${MACHTYPE}/NumericTuple.o \
   ${MACHTYPE}/IntegralTuple.o \
	 ${MACHTYPE}/PrintGraph.o \
	 ${MACHTYPE}/BVertex.o \
	 ${MACHTYPE}/RuleList.o \
	 ${MACHTYPE}/Voting.o \
	 ${MACHTYPE}/Scaffold.o \
   ${MACHTYPE}/AlternativeEdge.o 



PATH_LIBOBJS = \
	 ${MACHTYPE}/PathLib.o \
	 ${MACHTYPE}/PathBranch.o \
	 ${MACHTYPE}/MateLibrary.o 


DIR_LIST = $(sort $(dir $(wildcard */Makefile)))
DIR_BASE = $(subst /,, $(DIR_LIST))
SUBDIRS  = $(filter-out $(MACHTYPE), $(DIR_BASE))

MODE     = $(filter debug, $(MAKECMDGOALS))

# uncomment to profile
#CPPOPTS += -p

.PHONY: arch_dir $(EXECS) $(SUBDIRS) Makefile.execs


include Makefile.targets

all: $(CURDIR)/make.dep libpath libassemble Makefile.targets ${EXECS} $(SUBDIRS) 
	@echo "building execs: ${EXECS}"

# This builds both Makefile.targets and Makefile.execs
# This way only a subset of all the possible programs are built.

$(LIBCOMMON):	
	cd $(SRC_DIR)/lib && $(MAKE) $(MODE) 

$(LIBALIGN):
	cd $(SRC_DIR)/lib/align && $(MAKE) $(MODE) 

$(LIBHASH):
	cd $(SRC_DIR)/lib/hash && $(MAKE) $(MODE) 

$(LIBPARSEBLAST):
	cd $(SRC_DIR)/lib/parseblast && $(MAKE) $(MODE)


Makefile.targets: *.cpp ExecList.txt 
								@echo "building exec list" $(MAKECMDGOALS)
								../CreateExecBuildCommands.pl ExecList.txt
							  $(MAKE) $(MAKECMDGOALS) libassemble
								$(MAKE) $(MAKECMDGOALS) all 

define MAKE_DIRECTORY
$(1):
	cd $(1) && $${MAKE} $${MAKECMDGOALS}
endef

$(foreach subdir, $(SUBDIRS), $(eval $(call MAKE_DIRECTORY, $(subdir))))


LIBASSEMBLE = $(MACHTYPE)/libassemble.${LIB_EXT}
LIBPATH     = $(MACHTYPE)/libpath.${LIB_EXT}

libassemble: arch_dir $(LIBASSEMBLE)
libpath: arch_dir $(LIBPATH)

include ../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libassemble.${LIB_EXT}: ${LIBOBJS}  $(LIBCOMMON)
	$(LIB_TOOL) ${LIB_CREATE} $@ $^
	cp $(MACHTYPE)/libassemble.${LIB_EXT} ../lib/$(MACHTYPE)/


$(MACHTYPE)/libpath.${LIB_EXT}: $(PATH_LIBOBJS) $(MACHTYPE)/libassemble.${LIB_EXT} $(LIBCOMMON) Trace.h
	$(LIB_TOOL) ${LIB_CREATE} $@ $(PATH_LIBOBJS) #-lcommon $(LIBDIR) -lassemble
		cp $(MACHTYPE)/libpath.${LIB_EXT} ../lib/$(MACHTYPE)/libpath.$(LIB_EXT)

include Makefile.execs
$(LIBDIR)/libassemble.${LIB_EXT}: $(MACHTYPE)/libassemble.$(LIB_EXT)
		-cp $(MACHTYPE)/libpath.${LIB_EXT} ../lib/$(MACHTYPE)/libassemble.$(LIB_EXT)

