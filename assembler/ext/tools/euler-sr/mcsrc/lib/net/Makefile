include ../../common.mak


LIB_OBJS    = NetReader.o Net.o 
LIBNETDB_OBJS = NetDB.o

ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIB_OBJS))
DB_ARCH_OBJS = $(addprefix $(MACHTYPE)/,$(LIBNETDB_OBJS))
INCLUDEDIR += $(MYSQL_INC) -I/$(HOME)/sw/include


.PHONY: arch_dir

#ifeq ($(MACHTYPE), x86_64)
all: arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libnet.${LIB_EXT}
#else
#all: arch_dir $(CURDIR)/make.dep ../$(MACHTYPE)/libnet.${LIB_EXT} ../$(MACHTYPE)/libnetdb.${LIB_EXT}
#endif


include ../../make.rules
include $(CURDIR)/make.dep

$(MACHTYPE)/libnet.${LIB_EXT}: ${ARCH_OBJS}  \
	$(SRC_DIR)/lib/${MACHTYPE}/libcommon.$(LIB_EXT) $(LIBBLOCK)
	$(LIB_TOOL) $(LIB_CREATE) $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/libnet.${LIB_EXT} ../$(MACHTYPE)

$(MACHTYPE)/libnetdb.${LIB_EXT}: ${DB_ARCH_OBJS} \
	$(SRC_DIR)/lib/${MACHTYPE}/libcommon.$(LIB_EXT) 
	$(LIB_TOOL) $(LIB_CREATE) $@ ${DB_ARCH_OBJS}
	-cp $(MACHTYPE)/libnetdb.${LIB_EXT} ../$(MACHTYPE)

$(LIBDIR)/${MACHTYPE}/libcommon.$(LIB_EXT):
	cd .. && ${MAKE} common

$(LIBBLOCK):
	cd .. && ${MAKE} blocks
