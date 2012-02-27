include ../../common.mak

LIBOBJS    = SRSPairParser.o EmbossAlign.o EmbossAlignment.o


ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

.PHONY: emboss arch_dir lib_arch_dir

all: arch_dir lib_arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libemboss.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libemboss.${LIB_EXT}: ${ARCH_OBJS} ${SRC_DIR}/lib/${MACHTYPE}/libcommon.$(LIB_EXT) $(LIBCOMMON)
	$(LIB_TOOL) $(LIB_CREATE) $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/libemboss.${LIB_EXT} ../$(MACHTYPE)
$(LIBCOMMON):
	cd .. && ${MAKE} common


