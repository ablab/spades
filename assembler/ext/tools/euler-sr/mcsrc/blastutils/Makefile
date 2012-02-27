include ../common.mak

EXECS = profileBlastAlignments intersectProfiles printMap
.PHONY: arch_dir $(EXECS)

all: $(CURDIR)/make.dep $(EXECS) install

include ../make.rules
include $(CURDIR)/make.dep

profileBlastAlignments: arch_dir $(MACHTYPE)/profileBlastAlignments $(MACHTYPE)/intersectProfiles
printMap: arch_dir $(MACHTYPE)/printMap

LIBCOMMON = $(SRC_DIR)/lib/$(MACHTYPE)/libcommon.${LIB_EXT} 
LIBPB     = $(SRC_DIR)/lib/$(MACHTYPE)/libparseblast.${LIB_EXT} 

$(MACHTYPE)/profileBlastAlignments: ${MACHTYPE}/ProfileBlastAlignments.o $(LIBPB)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lparseblast -lcommon  $(LIBS)

$(MACHTYPE)/intersectProfiles: $(MACHTYPE)/IntersectProfiles.o
	$(CPP) $(CPPOPTS_lib) $< -o $@ 

$(MACHTYPE)/printMap: $(MACHTYPE)/PrintMap.o $(LIBPB)
	$(CPP) $(CPPOPTS_lib) $< -o $@  $(LIBDIR) -lparseblast  -lcommon

$(LIBPB):
	cd ../lib/blast && $(MAKE)
