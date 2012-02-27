include ../../common.mak

LIBOBJS    = HashUtils.o 

ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

.PHONY: 
debug: all

all: arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libhash.${LIB_EXT} ${SRC_DIR}/lib/${MACHTYPE}/libhash.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libhash.${LIB_EXT}: ${ARCH_OBJS}
	$(LIB_TOOL) $(LIB_CREATE) $@ ${ARCH_OBJS}

${SRC_DIR}/lib/$(MACHTYPE)/libhash.${LIB_EXT}:  $(MACHTYPE)/libhash.${LIB_EXT}
	-cp $(MACHTYPE)/libhash.${LIB_EXT} ${SRC_DIR}/lib/$(MACHTYPE)/libhash.${LIB_EXT}
