include ../../common.mak

LIBOBJS    =  BlastParser.o 


ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

.PHONY: parseblast arch_dir

all: arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libparseblast.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libparseblast.${LIB_EXT}: ${ARCH_OBJS} ${SRC_DIR}/lib/${MACHTYPE}/libcommon.$(LIB_EXT)
	$(LIB_TOOL) $(LIB_CREATE) $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/libparseblast.${LIB_EXT} ../$(MACHTYPE)

$(LIBDIR)/${MACHTYPE}/libcommon.$(LIB_EXT):
	cd .. && ${MAKE} common


