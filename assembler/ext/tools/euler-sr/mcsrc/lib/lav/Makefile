include ../../common.mak

LIBOBJS  = LAVReader.o \
	     LAVPrinter.o \
	     LAVTable.o \
	     LAVUtils.o \
	     LAVFile.o \
	     LAVBlock.o


ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))

.PHONY: lav arch_dir lib_arch_dir

all: arch_dir lib_arch_dir $(CURDIR)/make.dep $(MACHTYPE)/liblav.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/liblav.${LIB_EXT}: ${ARCH_OBJS}
	${LIB_TOOL} ${LIB_CREATE} $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/liblav.${LIB_EXT} ../$(MACHTYPE)/

