include ../../common.mak

LIBOBJS    = BlockLib.o \
	     BlockGraph.o \
	     BlockMap.o \
       BlockMapTree.o \
	     BlockPath.o \
	     Block.o 

BLOCKDBOBJS = BlockDB.o

ARCH_OBJS  = $(addprefix $(MACHTYPE)/,$(LIBOBJS))
BLOCKDB_ARCH_OBJS = $(addprefix $(MACHTYPE)/,$(BLOCKDBOBJS))

.PHONY: 


#ifeq ($(MACHTYPE), x86_64)
#all: arch_dir $(CURDIR)/make.dep ../$(MACHTYPE)/libblocks.${LIB_EXT} 
#else
#all: arch_dir $(CURDIR)/make.dep ../$(MACHTYPE)/libblocks.${LIB_EXT} ../$(MACHTYPE)/libblockdb.${LIB_EXT}
#endif
all: arch_dir $(CURDIR)/make.dep $(MACHTYPE)/libblocks.${LIB_EXT}

include ../../make.rules
include $(CURDIR)/make.dep

BlockDB.o: $(MACHTYPE)/BlockDB.o

$(MACHTYPE)/BlockDB.o: BlockDB.cpp BlockDB.h
	$(CPP) $(CPPOPTS) -c BlockDB.cpp -o $(MACHTYPE)/BlockDB.o $(MYSQL_INC) $(INCLUDEDIR)  

$(MACHTYPE)/libblocks.${LIB_EXT}: ${ARCH_OBJS}
	$(LIB_TOOL) $(LIB_CREATE) $@ ${ARCH_OBJS}
	-cp $(MACHTYPE)/libblocks.${LIB_EXT} ../$(MACHTYPE)/

$(MACHTYPE)/libblockdb.${LIB_EXT}: $(BLOCKDB_ARCH_OBJS) 
	$(LIB_TOOL) $(LIB_CREATE) $@ $(BLOCKDB_ARCH_OBJS)
	-cp $(MACHTYPE)/libblockdb.${LIB_EXT} ../(MACHTYPE)/

$(LIBDIR)/${MACHTYPE}/libcommon.$(LIB_EXT):
	cd .. && ${MAKE} common
