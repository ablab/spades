include ../../common.mak

LIBOBJS    = NewettTree.o

ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

.PHONY: lav arch_dir

all: arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libtree.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libtree.${LIB_EXT}: ${ARCH_OBJS}
	$(LIB_TOOL) $(LIB_CREATE) $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/libtree.${LIB_EXT} ../$(MACHTYPE)

