include ../../common.mak

LIBOBJS    = AxtReader.o \
	     AxtLib.o

ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

.PHONY: axt

all: arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libaxt.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libaxt.${LIB_EXT}: ${ARCH_OBJS}
	$(LIB_TOOL) $(LIB_CREATE) $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/libaxt.$(LIB_EXT) ../$(MACHTYPE)


