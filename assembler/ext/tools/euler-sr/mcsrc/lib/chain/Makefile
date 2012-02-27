include ../../common.mak

LIBOBJS    = ChainReader.o ChainUtils.o 

ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

.PHONY: lav arch_dir

all: arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libchain.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libchain.${LIB_EXT}: ${ARCH_OBJS}  $(LIBLAV)
	${LIB_TOOL} ${LIB_CREATE} $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/libchain.${LIB_EXT} ../$(MACHTYPE)

$(LIBLAV):
	cd ../lav && ${MAKE} 


