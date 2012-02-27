include ../../common.mak

HEADERS = DNASequence.h SeqReader.h alignutils.h modify_tuple.h graphalign.h


.PHONY:   libalign arch_dir lib_arch_dir

debug: all
all: arch_dir lib_arch_dir $(CURDIR)/make.dep libalign


LIBOBJS = alignutils.o graphalign.o
ARCH_OBJS = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

include ../../make.rules
include $(CURDIR)/make.dep

libalign:  arch_dir $(CURDIR)/make.dep ./$(MACHTYPE)/libalign.$(LIB_EXT) $(SRC_DIR)/lib/$(MACHTYPE)/libalign.$(LIB_EXT)

$(SRC_DIR)/lib/$(MACHTYPE)/libalign.$(LIB_EXT): $(MACHTYPE)/libalign.$(LIB_EXT)
	-cp $(MACHTYPE)/libalign.$(LIB_EXT) $@

$(MACHTYPE)/libalign.$(LIB_EXT): $(ARCH_OBJS) $(LIBCOMMON)
	$(LIB_TOOL) $(LIB_CREATE) $@  $(ARCH_OBJS)



$(LIBCOMMON):
	cd $(SRC_DIR)/lib && ${MAKE} $(MAKECMDGOALS) common	


