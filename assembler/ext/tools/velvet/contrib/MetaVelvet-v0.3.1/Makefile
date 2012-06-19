CC = gcc
CFLAGS = -Wall
DEBUG = -g
LDFLAGS = -lm
OPT = -O3
MAXKMERLENGTH=31
CATEGORIES=2
DEF = -D MAXKMERLENGTH=$(MAXKMERLENGTH) -D CATEGORIES=$(CATEGORIES)

Z_LIB_DIR=../../third-party/zlib-1.2.3
Z_LIB_FILES=$(Z_LIB_DIR)/*.o

# Mac OS users: uncomment the following lines
# Z_LIB_FILES=
# LDFLAGS = -lm -lz
# CFLAGS = -Wall -m64

# Sparc/Solaris users: uncomment the following line
# CFLAGS = -Wall -m64

OBJ = obj/tightString.o obj/run.o obj/splay.o obj/splayTable.o obj/graph.o obj/run2.o obj/fibHeap.o obj/fib.o obj/concatenatedGraph.o obj/passageMarker.o obj/graphStats.o obj/correctedGraph.o obj/dfib.o obj/dfibHeap.o obj/recycleBin.o obj/readSet.o obj/shortReadPairs.o obj/locallyCorrectedGraph.o obj/graphReConstruction.o obj/roadMap.o obj/preGraph.o obj/preGraphConstruction.o obj/concatenatedPreGraph.o obj/readCoherentGraph.o obj/crc.o obj/utility.o obj/kmer.o obj/scaffold.o
OBJDBG = $(subst obj,obj/dbg,$(OBJ))

default : cleanobj zlib obj meta-velveth meta-velvetg

clean :
	-rm obj/*.o obj/dbg/*.o ./meta-velvet* 
	cd $(Z_LIB_DIR) && make clean

cleanobj: 
	-rm obj/*.o obj/dbg/*.o 

zlib : 
	cd $(Z_LIB_DIR); ./configure; make; rm minigzip.o; rm example.o

meta-velveth : obj 
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o meta-velveth obj/tightString.o obj/run.o obj/recycleBin.o obj/splay.o obj/splayTable.o obj/readSet.o obj/crc.o obj/utility.o obj/kmer.o $(Z_LIB_FILES)


meta-velvetg : obj
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o meta-velvetg obj/tightString.o obj/graph.o obj/run2.o obj/fibHeap.o obj/fib.o obj/concatenatedGraph.o obj/passageMarker.o obj/graphStats.o obj/correctedGraph.o obj/dfib.o obj/dfibHeap.o obj/recycleBin.o obj/readSet.o obj/shortReadPairs.o obj/scaffold.o obj/locallyCorrectedGraph.o obj/graphReConstruction.o obj/roadMap.o obj/preGraph.o obj/preGraphConstruction.o obj/concatenatedPreGraph.o obj/readCoherentGraph.o obj/utility.o obj/kmer.o $(Z_LIB_FILES)

debug : cleanobj obj/dbg
	$(CC) $(CFLAGS) $(LDFLAGS) $(DEBUG) -o meta-velveth obj/dbg/tightString.o obj/dbg/run.o obj/dbg/recycleBin.o obj/dbg/splay.o obj/dbg/splayTable.o obj/dbg/readSet.o obj/dbg/crc.o obj/dbg/utility.o obj/dbg/kmer.o $(Z_LIB_FILES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DEBUG) -o meta-velvetg obj/dbg/tightString.o obj/dbg/graph.o obj/dbg/run2.o obj/dbg/fibHeap.o obj/dbg/fib.o obj/dbg/concatenatedGraph.o obj/dbg/passageMarker.o obj/dbg/graphStats.o obj/dbg/correctedGraph.o obj/dbg/dfib.o obj/dbg/dfibHeap.o obj/dbg/recycleBin.o obj/dbg/readSet.o obj/dbg/shortReadPairs.o obj/dbg/scaffold.o obj/dbg/locallyCorrectedGraph.o obj/dbg/graphReConstruction.o obj/dbg/roadMap.o obj/dbg/preGraph.o obj/dbg/preGraphConstruction.o obj/dbg/concatenatedPreGraph.o obj/dbg/readCoherentGraph.o obj/dbg/utility.o obj/dbg/kmer.o $(Z_LIB_FILES)

color : cleanobj obj_de
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o meta-velveth_de obj/tightString.o obj/run.o obj/recycleBin.o obj/splay.o obj/splayTable.o obj/readSet.o obj/crc.o obj/utility.o obj/kmer.o $(Z_LIB_FILES)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o velvetg_de obj/tightString.o obj/graph.o obj/run2.o obj/fibHeap.o obj/fib.o obj/concatenatedGraph.o obj/passageMarker.o obj/graphStats.o obj/correctedGraph.o obj/dfib.o obj/dfibHeap.o obj/recycleBin.o obj/readSet.o obj/shortReadPairs.o obj/scaffold.o obj/locallyCorrectedGraph.o obj/graphReConstruction.o obj/roadMap.o obj/preGraph.o obj/preGraphConstruction.o obj/concatenatedPreGraph.o obj/readCoherentGraph.o obj/utility.o obj/kmer.o $(Z_LIB_FILES)

colordebug : cleanobj obj/dbg_de
	$(CC) $(CFLAGS) $(LDFLAGS) $(DEBUG) -o meta-velveth_de obj/dbg/tightString.o obj/dbg/run.o obj/dbg/recycleBin.o obj/dbg/splay.o obj/dbg/splayTable.o obj/dbg/readSet.o obj/dbg/crc.o obj/dbg/utility.o obj/dbg/kmer.o $(Z_LIB_FILES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DEBUG) -o velvetg_de obj/dbg/tightString.o obj/dbg/graph.o obj/dbg/run2.o obj/dbg/fibHeap.o obj/dbg/fib.o obj/dbg/concatenatedGraph.o obj/dbg/passageMarker.o obj/dbg/graphStats.o obj/dbg/correctedGraph.o obj/dbg/dfib.o obj/dbg/dfibHeap.o obj/dbg/recycleBin.o obj/dbg/readSet.o obj/dbg/shortReadPairs.o obj/dbg/scaffold.o obj/dbg/locallyCorrectedGraph.o obj/dbg/graphReConstruction.o obj/dbg/roadMap.o obj/dbg/preGraph.o obj/dbg/preGraphConstruction.o obj/dbg/concatenatedPreGraph.o obj/dbg/readCoherentGraph.o obj/dbg/utility.o obj/dbg/kmer.o $(Z_LIB_FILES)

objdir:
	mkdir -p obj

obj: zlib objdir $(OBJ)

obj_de: override DEF := $(DEF) -D COLOR
obj_de: zlib cleanobj objdir $(OBJ)

obj/dbgdir: 
	mkdir -p obj/dbg

obj/dbg: zlib obj/dbgdir $(OBJDBG)

obj/dbg_de: override DEF := $(DEF) -D COLOR
obj/dbg_de: zlib cleanobj obj/dbgdir $(OBJDBG)

obj/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPT) $(DEF) -c $? -o $@ 

obj/dbg/%.o: src/%.c
	$(CC) $(CFLAGS) $(DEBUG) $(DEF) -c $? -o $@ 

